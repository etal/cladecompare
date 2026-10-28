"""Core API shared by cladecompare and cladeweb for running the algorithms."""

import contextlib
import logging
import os
import re
import subprocess
import sys
import tempfile
from copy import deepcopy
from io import StringIO
from math import fsum, log10
from os.path import basename

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biofrills import alnutils, consensus

from . import InputError, gtest, hypg, jsd, pairlogo, phospho, pmlscript, urn
from .shared import combined_frequencies

logger = logging.getLogger(__name__)

GAP_THRESH = 0.8
PSEUDO_SIZE = 0.5

# The contrast statistics, by the name that -s and the web form give, each
# with the description that the help text, the web form and the log show.
STRATEGIES = {
    "gtest": ("G-test of amino acid frequencies", gtest),
    "urn": ("ball-in-urn model of consensus residue conservation", urn),
    "phospho": ("urn model for phosphorylatable residues (S/T/Y)", phospho),
    "jsd": ("Jensen-Shannon divergence", jsd),
    "hypg": ("hypergeometric model", hypg),
}
DEFAULT_STRATEGY = "gtest"


def parse_alpha(text):
    """Parse a significance threshold, which must be between 0 and 1.

    Text that is not a number raises ValueError, and a number outside the
    range raises InputError, which is also a ValueError.
    """
    value = float(text)
    if not 0.0 <= value <= 1.0:
        raise InputError(f"must be between 0 and 1: {text!r}")
    return value


# ---- FLOW --------------------------------------------------------------


def process_args(args):
    """Main."""
    if (args.mapgaps or args.hmm) and args.format != "fasta":
        # run_gaps requires FASTA input
        raise InputError("Sequences to align with a profile must be in FASTA format.")

    all_alns = []
    for alnfname in [args.foreground] + args.background:
        if args.hmm:
            logger.info("Aligning %s with HMM profile %s", alnfname, args.hmm)
            aln = hmm_align_and_read(args.hmm, alnfname)
        elif args.mapgaps:
            logger.info("Aligning %s with MAPGAPS profile %s", alnfname, args.mapgaps)
            aln = mapgaps_align_and_read(args.mapgaps, alnfname)
        else:
            aln = read_aln(alnfname, args.format)
        all_alns.append(aln)

    pdb_data = []
    if args.pdb:
        if args.hmm:
            for pdbfname in args.pdb:
                logger.info("Aligning %s with HMM profile %s", pdbfname, args.hmm)
                pdb_rec, pdb_resnums, pdb_inserts = pdb_hmm(args.hmm, pdbfname)
                pdb_data.append((pdbfname, pdb_rec, pdb_resnums, pdb_inserts))
        elif args.mapgaps:
            for pdbfname in args.pdb:
                logger.info(
                    "Aligning %s with MAPGAPS profile %s", pdbfname, args.mapgaps
                )
                pdb_rec, pdb_resnums, pdb_inserts = pdb_mapgaps(args.mapgaps, pdbfname)
                pdb_data.append((pdbfname, pdb_rec, pdb_resnums, pdb_inserts))
        else:
            logger.error("PDB alignment requires a MAPGAPS or HMM profile.")
            # ENH - realign to fg, bg

    description, module = STRATEGIES[args.strategy]
    logger.info("Using %s", description)
    output = args.output or sys.stdout

    if len(all_alns) == 1:
        aln, hits = process_one(all_alns[0], module, args.weight)
        process_output(
            aln,
            None,
            hits,
            args.alpha,
            output,
            args.pattern,
            pdb_data,
            args.pmlout,
        )
    elif len(all_alns) == 2:
        fg_clean, bg_clean, hits = process_pair(
            all_alns[0], all_alns[1], module, args.weight
        )
        process_output(
            fg_clean,
            bg_clean,
            hits,
            args.alpha,
            output,
            args.pattern,
            pdb_data,
            args.pmlout,
        )
    else:
        # Output fnames are based on fg filenames; ignore what's given
        if args.output or args.pattern:
            logger.warning(
                "Comparing more than two alignments: ignoring -o and -p, and "
                "writing <alignment>.out and <alignment>.pttrn for each alignment"
            )
        outfnames_ptnfnames = [
            (basename(alnfname) + ".out", basename(alnfname) + ".pttrn")
            for alnfname in ([args.foreground] + args.background)
        ]
        for idx, fg_aln in enumerate(all_alns):
            # Combine other alns into bg
            _other_alns = all_alns[:idx] + all_alns[idx + 1 :]
            bg_aln = deepcopy(_other_alns[0])
            for otra in _other_alns[1:]:
                bg_aln.extend(deepcopy(otra))
            fg_clean, bg_clean, hits = process_pair(
                deepcopy(fg_aln), bg_aln, module, args.weight
            )
            outfname, ptnfname = outfnames_ptnfnames[idx]
            process_output(
                fg_clean,
                bg_clean,
                hits,
                args.alpha,
                outfname,
                ptnfname,
                pdb_data,
                args.pmlout,
            )
            logger.info("Wrote %s and %s", outfname, ptnfname)


def process_pair(fg_aln, bg_aln, module, do_weight):
    """Calculate a mapping of alignment column positions to "contrast".

    Return a list of tuples:
        (foreground consensus aa, background consensus aa, p-value)
        for each column position.
    """
    fg_aln, bg_aln = clean_alignments(fg_aln, bg_aln)
    if do_weight:
        fg_weights = alnutils.sequence_weights(fg_aln, "none")
        bg_weights = alnutils.sequence_weights(bg_aln, "none")
    else:
        fg_weights = [1 for i in range(len(fg_aln))]
        bg_weights = [1 for i in range(len(bg_aln))]
    fg_size = fsum(fg_weights) if module != urn else len(fg_aln)
    bg_size = fsum(bg_weights)
    # Overall aa freqs for pseudocounts
    aa_freqs = combined_frequencies(fg_aln, fg_weights, bg_aln, bg_weights)
    fg_cons = consensus.consensus(
        fg_aln, weights=fg_weights, trim_ends=False, gap_threshold=GAP_THRESH
    )
    bg_cons = consensus.consensus(
        bg_aln, weights=bg_weights, trim_ends=False, gap_threshold=GAP_THRESH
    )

    hits = []
    for faa, baa, fg_col, bg_col in zip(
        fg_cons,
        bg_cons,
        zip(*fg_aln, strict=True),
        zip(*bg_aln, strict=True),
        strict=True,
    ):
        if faa == "-" or baa == "-":
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.0
        else:
            pvalue = module.compare_cols(
                fg_col,
                faa,
                fg_size,
                fg_weights,
                bg_col,
                baa,
                bg_size,
                bg_weights,
                aa_freqs,
                PSEUDO_SIZE,
            )
        hits.append((faa, baa, pvalue))

    return fg_aln, bg_aln, hits


def process_one(aln, module, do_weight):
    """Calculate a mapping of alignment column positions to "contrast"."""
    if do_weight:
        weights = alnutils.sequence_weights(aln, "none")
    else:
        weights = [1 for i in range(len(aln))]
    aln_size = fsum(weights) if module != urn else len(aln)
    aa_freqs = alnutils.aa_frequencies(aln, weights, gap_chars="-.X")
    cons = consensus.consensus(
        aln, weights=weights, trim_ends=False, gap_threshold=GAP_THRESH
    )
    hits = []
    for cons_aa, col in zip(cons, zip(*aln, strict=True), strict=True):
        if cons_aa == "-":
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.0
        else:
            pvalue = module.compare_one(
                col, cons_aa, aln_size, weights, aa_freqs, PSEUDO_SIZE
            )
        hits.append((cons_aa, "_", pvalue))
    return aln, hits


def process_output(
    fg_aln, bg_aln, hits, alpha, output, pattern, pdb_data, pml_output=None
):
    """Generate the output files from the processed data."""
    with as_handle(output, "w+") as outfile:
        write_pvalues(hits, outfile, alpha)
    tophits = top_hits(hits, alpha)
    if pattern:
        with open(pattern, "w+") as ptnfile:
            write_mcbpps(tophits, ptnfile)
        # XXX hack: don't make pairlogo in single mode
        if bg_aln:
            pairlogo.make_pairlogos(
                fg_aln, bg_aln, tophits, pattern.rsplit(".", 1)[0], 10
            )
    if pdb_data:
        patterns = [t[0] for t in tophits]
        if len(pdb_data) == 1:
            # Single-PBD mode
            pdb_fname, pdb_rec, pdb_resnums, pdb_inserts = pdb_data[0]
            script = pmlscript.build_single(
                pdb_resnums,
                pdb_inserts,
                patterns,
                pdb_fname,
                pdb_rec.annotations["chain"],
            )
            if not pml_output:
                pml_output = pdb_fname + ".pml"
        else:
            # Multi-PBD mode
            pdb_fnames, pdb_recs, pdb_resnumses, pdb_insertses = zip(
                *pdb_data, strict=True
            )
            script = pmlscript.build_multi(
                pdb_resnumses,
                pdb_insertses,
                patterns,
                pdb_fnames,
                [r.annotations["chain"] for r in pdb_recs],
            )
            if not pml_output:
                pml_output = pdb_fnames[0] + "-etc.pml"
        with open(pml_output, "w+") as pmlfile:
            pmlfile.write(script)
        logger.info("Wrote %s", pml_output)


# --- Output ---


def write_pvalues(hits, outfile, alpha):
    """Write p-values & "contrast" stars for each site. (It's noisy.)"""
    for idx, data in enumerate(hits):
        fg_char, bg_char, pvalue = data
        if not (0.0 <= pvalue <= 1.0):
            logger.warning("Out-of-domain p-value at site %s: %s", idx, pvalue)
        stars = "*" * int(-log10(pvalue)) if 0 < pvalue < alpha else ""
        outfile.write(
            "%s (%s) %d : prob=%g\t%s\n" % (fg_char, bg_char, idx + 1, pvalue, stars)
        )


def write_mcbpps(tophits, ptnfile):
    """Write a .pttrn file in the style of mcBPPS."""
    ptnfile.write(
        "1:" + ",".join([("%s%d" % (faa, posn)) for posn, faa, baa in tophits])
    )


# --- Input magic ---


def call_quiet(*args):
    """Safely run a command and get stdout; print stderr if there's an error.

    Like subprocess.check_output, but silent in the normal case where the
    command logs unimportant stuff to stderr. If there is an error, then the
    full error message(s) is shown in the exception message. If the program
    is not found, a RuntimeError naming it is raised as well.
    """
    try:
        proc = subprocess.run(
            args, capture_output=True, encoding="utf-8", errors="replace"
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Subprocess command failed:\n$ %s\n\n%s\n"
            "Is %s installed and on the PATH?" % (" ".join(args), exc, args[0])
        ) from None
    if proc.returncode != 0:
        raise RuntimeError(
            "Subprocess command failed:\n$ %s\n\n%s"
            % (" ".join(args), proc.stdout + proc.stderr)
        )
    return proc.stdout


def hmmalign_a2m(hmm_profile, seq_fname):
    """Align sequences to an HMM profile and return hmmalign's a2m output.

    Every match state of the profile becomes a column, including one that is
    a gap in every sequence, so that sequence sets aligned separately to the
    same profile share the same columns. HMMER 3.1 and later always keep
    such columns and reject the --allcol option; HMMER 3.0 keeps them only
    with --allcol, so the option is given when hmmalign's help lists it.
    """
    allcol = ["--allcol"] if "--allcol" in call_quiet("hmmalign", "-h") else []
    return call_quiet(
        "hmmalign",
        *allcol,
        "--trim",
        "--amino",
        "--outformat",
        "a2m",
        hmm_profile,
        seq_fname,
    )


def hmm_align_and_read(hmm_profile, fasta_fname):
    """Align a FASTA file with HMMER 3 and read the alignment."""
    out = hmmalign_a2m(hmm_profile, fasta_fname)
    # ENH: write to file, then parse incrementally
    records = list(SeqIO.parse(StringIO(out), "fasta"))
    # Remove inserts, i.e. lowercase characters
    for rec in records:
        rec.seq = Seq("".join([c for c in str(rec.seq) if not c.islower()]))
    return MultipleSeqAlignment(records)


def mapgaps_align_and_read(mapgaps_profile, fasta_fname):
    """Align a FASTA file with MAPGAPS and read the CMA alignment."""
    call_quiet("run_gaps", mapgaps_profile, fasta_fname)
    aln = cma_blocks(fasta_fname + "_aln.cma")
    return aln


# A CMA sequence line: the aligned sequence between '{(...)' and '(...)}',
# where flanking residues may appear outside the braces or inside the
# parentheses, e.g.
#   {()SEQ()}*   MENINNI{()SEQ()}*   {(HY)SEQ(AA)}*
# The aligned part is taken verbatim, so that gap characters such as '.' pass
# through as before; the flanks must be residues, so that other text around
# the block is reported as an error.
CMA_SEQ_LINE = re.compile(
    r"[A-Za-z-]*\{\([A-Za-z-]*\)(?P<aligned>[^()]*)"
    r"\([A-Za-z-]*\)\}[A-Za-z-]*\*\Z"
)


def cma_records(cma_fname):
    """Read the records of a CMA file.

    Return a list of (ID, aligned sequence) pairs. The aligned sequence is
    taken verbatim, with lowercase inserts, and without the flanking residues
    that some records carry around it. A sequence line that does not have the
    expected shape raises InputError.
    """
    records = []
    with open(cma_fname) as infile:
        lines = iter(infile)
        for line in lines:
            if line.startswith(">"):
                acc = line.split(None, 1)[0][1:]
                seq_line = next(lines, "").strip()
                match = CMA_SEQ_LINE.match(seq_line)
                if match is None:
                    raise InputError(
                        "Unparseable CMA sequence line for "
                        "record %s in %s: %s" % (acc, cma_fname, seq_line)
                    )
                records.append((acc, match.group("aligned")))
    return records


def cma_blocks(cma_fname):
    """Select the conserved/consensus columns in the alignment.

    This removes inserts relative to the consensus sequence, and the flanking
    residues that some records carry around the aligned sequence. A sequence
    line that does not have the expected shape raises InputError.

    Return a Biopython MultipleSeqAlignment.
    """
    return MultipleSeqAlignment(
        [
            SeqRecord(
                Seq("".join(c for c in aligned if not c.islower())),
                id=acc,
                description="",
            )
            for acc, aligned in cma_records(cma_fname)
        ]
    )


def read_aln(fname, format):
    """Read a sequence alignment."""
    if fname.endswith(".cma") or format == "cma":
        return cma_blocks(fname)
    else:
        aln = AlignIO.read(fname, format)
        for rec in aln:
            rec.description = ""
        return aln


# The program name and version at the start of `muscle -version`, e.g.
# "MUSCLE v3.8.1551 by Robert C. Edgar" or "muscle 5.3.osxarm64 []".
MUSCLE_VERSION_LINE = re.compile(r"muscle +v?(\d+)\.(\d+)", re.IGNORECASE)


def combine_alignments(fg_aln, bg_aln):
    """Align FG and BG to each other so column numbers match.

    Uses MUSCLE 3 for profile-profile alignment, which MUSCLE 5 no longer
    offers; a MUSCLE of version 4 or later raises RuntimeError. If the
    version cannot be read, a warning is logged and MUSCLE is run anyway.
    """
    match = MUSCLE_VERSION_LINE.search(call_quiet("muscle", "-version"))
    if match is None:
        logger.warning(
            "Could not determine the version of MUSCLE; "
            "trying profile-profile alignment anyway"
        )
    elif int(match.group(1)) >= 4:
        raise RuntimeError(
            "Aligning alignments of different widths to each "
            "other requires MUSCLE 3, but the muscle on the "
            "PATH is version %s.%s, which has no "
            "profile-profile alignment mode. Install a "
            "MUSCLE 3 release such as 3.8.1551, or give "
            "alignments of equal width." % match.groups()
        )
    with tempfile.TemporaryDirectory() as tmpdir:
        aseqfname = os.path.join(tmpdir, "fg.fasta")
        bseqfname = os.path.join(tmpdir, "bg.fasta")
        AlignIO.write(fg_aln, aseqfname, "fasta")
        AlignIO.write(bg_aln, bseqfname, "fasta")
        output = call_quiet("muscle", "-profile", "-in1", aseqfname, "-in2", bseqfname)

    full_aln = AlignIO.read(StringIO(output), "fasta")
    full_aln = MultipleSeqAlignment(alnutils.remove_empty_cols(full_aln))
    return full_aln


def clean_alignments(fg_aln, bg_aln):
    """Fix simple issues in the alignments:

    - Remove duplicated sequences and IDs from the background
    - Ensure alignments are the same width (if not, align to each other)
    - Remove all-gap columns from the full alignment
    """
    if not len(fg_aln):
        raise InputError("Foreground set is empty")
    if not len(bg_aln):
        raise InputError("Background set is empty")

    # Remove FG seqs from BG -- by equal sequence IDs, here
    killme = []
    for fg_seq in fg_aln:
        for idx, bg_seq in enumerate(bg_aln):
            if fg_seq.id == bg_seq.id != "consensus":
                if str(fg_seq.seq) != str(bg_seq.seq):
                    logger.warning("Different sequences for %s in fg, bg", fg_seq.id)
                killme.append(idx)
    if killme:
        logger.info(
            "Removing %d duplicated sequence IDs from the background", len(killme)
        )
        for idx in sorted(killme, reverse=True):
            del bg_aln._records[idx]

    if not len(bg_aln):
        raise InputError(
            "Background set is a subset of the foreground. "
            "To fix this, use a background containing sequences "
            "not in the foreground."
        )

    # Remove identical sequences from the FG, then from the BG. Both calls
    # share one set, so a BG sequence identical to an FG sequence is removed
    # from the BG even under a different ID.
    seen = set()

    def purge_duplicates(aln, label):
        out_recs = []
        for rec in aln:
            seq = str(rec.seq)
            if seq not in seen:
                out_recs.append(rec)
                seen.add(seq)
        n_purged = len(aln) - len(out_recs)
        if n_purged:
            logger.warning(
                "Purging %d identical sequences from the %s", n_purged, label
            )
        return out_recs

    fg_aln._records = purge_duplicates(fg_aln, "foreground")
    bg_aln._records = purge_duplicates(bg_aln, "background")
    if not len(bg_aln):
        raise InputError(
            "Every background sequence is identical to a foreground sequence. "
            "To fix this, use a background containing sequences "
            "not in the foreground."
        )

    # Ensure alignments are the same width
    if len(fg_aln[0]) != len(bg_aln[0]):
        logger.warning("Alignments are not of equal width; fixing with MUSCLE.")
        full_aln = combine_alignments(fg_aln, bg_aln)
    else:
        full_aln = deepcopy(fg_aln)
        full_aln.extend(bg_aln)

    # Remove columns that are all gaps in both fg and bg (full_aln)
    seqstrs = [str(rec.seq) for rec in full_aln]
    clean_cols = [
        col for col in zip(*seqstrs, strict=True) if not all(c == "-" for c in col)
    ]
    clean_seqs = ["".join(row) for row in zip(*clean_cols, strict=True)]
    for rec, clean_seq in zip(full_aln, clean_seqs, strict=True):
        rec.seq = Seq(clean_seq)

    # Split the full alignment back into FG and BG sets
    fg_labels = set([seq.id for seq in fg_aln])
    fg_recs = []
    bg_recs = []
    for rec in full_aln:
        if rec.id in fg_labels:
            fg_recs.append(rec)
        else:
            bg_recs.append(rec)
    fg_aln._records = fg_recs
    bg_aln._records = bg_recs
    return fg_aln, bg_aln


def top_hits(hits, alpha, N=50):
    """Take the top (up to N) hits with corrected p-value <= alpha.

    Return a list of triplets, sorted by significance:
        (position, fg_aa, bg_aa)
    """
    hit_quads = [
        (i + 1, faa_baa_pval[0], faa_baa_pval[1], faa_baa_pval[2])
        for i, faa_baa_pval in enumerate(hits)
    ]

    def get_pval(ifbp):
        return ifbp[3]

    # Benjamini-Hochberg multiple-testing FDR correction (BH step-up)
    hit_quads.sort(key=get_pval)
    m = len(hit_quads)  # Num. hypotheses tested
    if m < N:
        N = m
    tophits = [(posn, faa, baa) for posn, faa, baa, pval in hit_quads]
    for k, ifbp in list(zip(range(m, 0, -1), reversed(hit_quads), strict=True))[-N:]:
        if get_pval(ifbp) <= alpha * k / m:
            return tophits[:k]
    return []


# --- PDB alignment magic ---


@contextlib.contextmanager
def read_pdb_seq(pdb_fname):
    """Write the PDB atom sequence(s) to a FASTA file in a temporary directory.

    Yields the name of the FASTA file and a list of the sequences themselves
    (as SeqRecords).

    Context manager: on exit, the directory is removed together with any files
    that an aligner wrote beside the FASTA file.
    """
    pdbseqs = list(SeqIO.parse(pdb_fname, "pdb-atom"))
    with tempfile.TemporaryDirectory() as tmpdir:
        pdbseqfname = os.path.join(tmpdir, "pdbseq")
        SeqIO.write(pdbseqs, pdbseqfname, "fasta")
        yield pdbseqfname, pdbseqs


def choose_best_aligned(aligned):
    """Choose the longest profile match as the "reference" chain.

    Returns a tuple: (sequence ID, sequence string)
    """

    def aligned_len(seq):
        return sum(c.isupper() for c in seq.replace("X", ""))

    if not aligned:
        raise RuntimeError("No PDB sequences were aligned by the profile!")
    elif len(aligned) == 1:
        ref_id, ref_aln = list(aligned.items())[0]
    else:
        ref_id = max(aligned.items(), key=lambda kv: aligned_len(kv[1]))
        ref_aln = aligned[ref_id]
    return ref_id, ref_aln


def get_aln_offset(full, aln):
    aln_trim = aln.replace("-", "").replace(".", "").upper()
    if "X" in aln_trim:
        aln_trim = aln_trim[: aln_trim.index("X")]
    return full.index(aln_trim)


def aln_resnums_inserts(aln, offset):
    """Return two lists: residue numbers for model columns; inserts."""
    aln_resnums = []
    aln_inserts = []
    in_insert = False
    del_len = 0
    curr_ins_start = None
    for i, c in enumerate(aln):
        if c.islower():
            if not in_insert:
                # Start of a new insert region
                curr_ins_start = offset + i + 1 - del_len
                in_insert = True
            continue

        if in_insert:
            # End of the current insert region
            aln_inserts.append((curr_ins_start, offset + i - del_len))
            in_insert = False

        if c.isupper():
            # Match position
            aln_resnums.append((c, offset + i - del_len + 1))
        elif c == "-":
            # Deletion position
            aln_resnums.append(("-", None))
            del_len += 1
        else:
            raise ValueError("Unexpected character '%s'" % c)
    return aln_resnums, aln_inserts


def pdb_ref_resnums(aligned, seqs):
    """Choose the reference chain and map the alignment to its residues.

    Returns a tuple: (SeqRecord,
                      list of aligned residue numbers,
                      list of insert ranges as tuple pairs)
    """
    ref_id, ref_aln = choose_best_aligned(aligned)
    ref_record = SeqIO.to_dict(seqs)[ref_id]
    # Calculate aligned residue numbers & insert ranges
    offset = (
        ref_record.annotations["start"]
        + get_aln_offset(str(ref_record.seq), ref_aln)
        - 1
    )
    resnums, inserts = aln_resnums_inserts(ref_aln, offset)
    return ref_record, resnums, inserts


def pdb_hmm(hmm_profile, pdb_fname):
    """Align a PDB structure to an HMM profile.

    Returns a tuple as pdb_ref_resnums does.
    """
    with read_pdb_seq(pdb_fname) as (seqfname, seqs):
        out = hmmalign_a2m(hmm_profile, seqfname)
    aligned = dict(
        (rec.id, str(rec.seq)) for rec in SeqIO.parse(StringIO(out), "fasta")
    )
    return pdb_ref_resnums(aligned, seqs)


def pdb_mapgaps(mapgaps_profile, pdb_fname):
    """Align a PDB structure to a MAPGAPS profile.

    Returns a tuple as pdb_ref_resnums does.
    """
    with read_pdb_seq(pdb_fname) as (seqfname, seqs):
        call_quiet("run_gaps", mapgaps_profile, seqfname)
        aligned = dict(cma_records(seqfname + "_aln.cma"))
    return pdb_ref_resnums(aligned, seqs)
