"""Core API shared by cladecompare and cladeweb for running the algorithms."""

import contextlib
import logging
import math
import os
import subprocess
import tempfile
from cStringIO import StringIO
from copy import deepcopy
from os.path import basename

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.File import as_handle

from biofrills import alnutils, consensus

from . import pairlogo, pmlscript, urn, gtest, jsd, phospho
from .shared import combined_frequencies


GAP_THRESH = 0.8
PSEUDO_SIZE = 0.5


# ---- FLOW --------------------------------------------------------------

def process_args(args):
    """Main."""
    if args.mapgaps or args.hmm:
        # run_gaps requires FASTA input
        assert args.format == 'fasta', \
                "Input sequence format must be FASTA."

    all_alns = []
    for alnfname in [args.foreground] + args.background:
        if args.hmm:
            logging.info("Aligning %s with HMM profile %s",
                         alnfname, args.hmm)
            aln = hmm_align_and_read(args.hmm, alnfname)
        elif args.mapgaps:
            logging.info("Aligning %s with MAPGAPS profile %s",
                         alnfname, args.mapgaps)
            aln = mapgaps_align_and_read(args.mapgaps, alnfname)
        else:
            aln = read_aln(alnfname, args.format)
        all_alns.append(aln)

    pdb_data = []
    if args.pdb:
        if args.hmm:
            for pdbfname in args.pdb:
                logging.info("Aligning %s with HMM profile %s",
                                pdbfname, args.hmm)
                pdb_rec, pdb_resnums, pdb_inserts = pdb_hmm(args.hmm, pdbfname)
                pdb_data.append((pdbfname, pdb_rec, pdb_resnums, pdb_inserts))
        elif args.mapgaps:
            for pdbfname in args.pdb:
                logging.info("Aligning %s with MAPGAPS profile %s",
                                pdbfname, args.mapgaps)
                pdb_rec, pdb_resnums, pdb_inserts = pdb_mapgaps(args.mapgaps,
                                                                pdbfname)
                pdb_data.append((pdbfname, pdb_rec, pdb_resnums, pdb_inserts))
        else:
            logging.error("PDB alignment requires a MAPGAPS or HMM profile.")
            # ENH - realign to fg, bg
            # aln = read_aln(args.pdb, 'pdb-atom')

    # ENH - put strategies in a dict, look up here
    if args.strategy == 'gtest':
        logging.info("Using G-test of amino acid frequencies")
        module = gtest
    elif args.strategy == 'urn':
        logging.info("Using ball-in-urn statistical model")
        module = urn
    elif args.strategy == 'jsd':
        logging.info("Using Jensen-Shannon divergence")
        module = jsd
    elif args.strategy == 'phospho':
        logging.info("Using urn model for phosphorylatable residues")
        module = phospho
    else:
        raise ValueError("Unknown strategy: %s" % args.strategy)

    if len(all_alns) == 1:
        aln, hits = process_one(all_alns[0], module)
        process_output(aln, None, hits, args.alpha, args.output, args.pattern,
                       pdb_data)
    elif len(all_alns) == 2:
        fg_clean, bg_clean, hits = process_pair(all_alns[0], all_alns[1],
                                                module)
        process_output(fg_clean, bg_clean, hits, args.alpha,
                       args.output, args.pattern,
                       pdb_data)
                       # args.pdb, pdb_rec, pdb_resnums, pdb_inserts)
    else:
        # Output fnames are based on fg filenames; ignore what's given
        outfnames_ptnfnames = [(basename(alnfname) + '.out',
                                basename(alnfname) + '.pttrn')
                               for alnfname in ([args.foreground] +
                                                args.background)]
        for idx, fg_aln in enumerate(all_alns):
            # Combine other alns into bg
            _other_alns = all_alns[:idx] + all_alns[idx+1:]
            bg_aln = deepcopy(_other_alns[0])
            for otra in _other_alns[1:]:
                bg_aln.extend(deepcopy(otra))
            fg_clean, bg_clean, hits = process_pair(deepcopy(fg_aln), bg_aln,
                                                    module)
            outfname, ptnfname = outfnames_ptnfnames[idx]
            process_output(fg_clean, bg_clean, hits, args.alpha,
                           outfname, ptnfname, pdb_data)
                           # args.pdb,
                           # pdb_rec, pdb_resnums, pdb_inserts)
            logging.info("Wrote %s and %s", outfname, ptnfname)


def process_pair(fg_aln, bg_aln, module):
    """Calculate a mapping of alignment column positions to "contrast".

    Return a list of tuples:
        (foreground consensus aa, background consensus aa, p-value)
        for each column position.
    """
    fg_aln, bg_aln = clean_alignments(fg_aln, bg_aln)
    fg_weights = alnutils.sequence_weights(fg_aln, 'none')
                                           # if module != jsd else 'sum1')
    fg_size = sum(fg_weights) if module != urn else len(fg_aln)
    bg_weights = alnutils.sequence_weights(bg_aln, 'none')
                                           # if module != jsd else 'sum1')
    bg_size = sum(bg_weights)
    # Overall aa freqs for pseudocounts
    aa_freqs = combined_frequencies(fg_aln, fg_weights, bg_aln, bg_weights)
    fg_cons = consensus.consensus(fg_aln, weights=fg_weights, trim_ends=False,
                                  gap_threshold=GAP_THRESH)
    bg_cons = consensus.consensus(bg_aln, weights=bg_weights, trim_ends=False,
                                  gap_threshold=GAP_THRESH)

    hits = []
    for faa, baa, fg_col, bg_col in zip(fg_cons, bg_cons,
                                        zip(*fg_aln), zip(*bg_aln)):
        if faa == '-' or baa == '-':
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.
        else:
            pvalue = module.compare_cols(
                fg_col, faa, fg_size, fg_weights,
                bg_col, baa, bg_size, bg_weights,
                aa_freqs, PSEUDO_SIZE)
        hits.append((faa, baa, pvalue))

    return fg_aln, bg_aln, hits


def process_one(aln, module):
    """Calculate a mapping of alignment column positions to "contrast"."""
    weights = alnutils.sequence_weights(aln, 'none')
                                        # if module != jsd else 'sum1')
    aln_size = sum(weights) if module != urn else len(aln)
    aa_freqs = alnutils.aa_frequencies(aln, weights, gap_chars='-.X')
    cons = consensus.consensus(aln, weights=weights, trim_ends=False,
                               gap_threshold=GAP_THRESH)
    hits = []
    for cons_aa, col in zip(cons, zip(*aln)):
        if cons_aa == '-':
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.
        else:
            pvalue = module.compare_one(col, cons_aa, aln_size, weights,
                                        aa_freqs, PSEUDO_SIZE)
        hits.append((cons_aa, '_', pvalue))
    return aln, hits


def process_output(fg_aln, bg_aln, hits, alpha, output, pattern, pdb_data):
    """Generate the output files from the processed data."""
    with as_handle(output, 'w+') as outfile:
        write_pvalues(hits, outfile, alpha)
    tophits = top_hits(hits, alpha)
    if pattern:
        with open(pattern, 'w+') as ptnfile:
            write_mcbpps(tophits, ptnfile)
        # XXX hack: don't make pairlogo in single mode
        if bg_aln:
            pairlogo.make_pairlogos(fg_aln, bg_aln, tophits,
                                    pattern.rsplit('.', 1)[0],
                                    10)
    if pdb_data:
        patterns = [t[0] for t in tophits]
        if len(pdb_data) == 1:
            pdb_fname, pdb_rec, pdb_resnums, pdb_inserts = pdb_data[0]
            script = pmlscript.build_single(pdb_resnums, pdb_inserts,
                                            patterns, pdb_fname,
                                            pdb_rec.annotations['chain'])
            pml_fname = pdb_fname + ".pml"
        else:
            pdb_fnames, pdb_recs, pdb_resnumses, pdb_insertses = zip(*pdb_data)
            # TODO multi-PDB mode
            pml_fname = pdb_fnames[0] + "-etc.pml"
        with open(pml_fname, 'w+') as pmlfile:
            pmlfile.write(script)
        logging.info("Wrote %s", pml_fname)



# --- Output ---

def write_pvalues(hits, outfile, alpha):
    """Write p-values & "contrast" stars for each site. (It's noisy.)"""
    for idx, data in enumerate(hits):
        fg_char, bg_char, pvalue = data
        if not (0.0 <= pvalue <= 1.0):
            logging.warn("Out-of-domain p-value at site %s: %s",
                         idx, pvalue)
        stars = ('*'*int(-math.log10(pvalue)) if 0 < pvalue < alpha else '')
        outfile.write("%s (%s) %d : prob=%g\t%s\n"
                      % (fg_char, bg_char, idx + 1, pvalue, stars))


def write_mcbpps(tophits, ptnfile):
    """Write a .pttrn file in the style of mcBPPS."""
    ptnfile.write("1:" + ','.join([("%s%d" % (faa, posn))
                                   for posn, faa, baa in tophits]))



# --- Input magic ---

def hmm_align_and_read(hmm_profile, fasta_fname):
    """Align a FASTA file with HMMer 3 and read the alignment."""
    out = subprocess.check_output(['hmmalign', '--allcol', '--trim', '--amino',
                                   '--outformat', 'a2m',
                                   hmm_profile, fasta_fname])
    # ENH: write to file, then parse incrementally
    records = list(SeqIO.parse(StringIO(out), 'fasta'))
    # Remove inserts, i.e. lowercase characters
    for rec in records:
        rec.seq._data = ''.join([c for c in str(rec.seq) if not c.islower()])
    return MultipleSeqAlignment(records, generic_protein)


def mapgaps_align_and_read(mapgaps_profile, fasta_fname):
    """Align a FASTA file with MAPGAPS and read the CMA alignment."""
    subprocess.check_call(['run_gaps', mapgaps_profile, fasta_fname])
    aln = cma_blocks(fasta_fname + '_aln.cma')
    return aln


def cma_blocks(cma_fname):
    """Select the conserved/consensus columns in the alignment.

    This removes inserts relative to the consensus sequence.

    Return a Biopython MultipleSeqAlignment.
    """
    records = []
    with open(cma_fname) as infile:
        lines = iter(infile)
        for line in lines:
            if line.startswith('>'):
                acc = line.split(None, 1)[0][1:]
                seq = next(lines).strip()
                seq = ''.join((c for c in seq[3:-4] if not c.islower()))
                records.append(
                        SeqRecord(Seq(seq, generic_protein),
                                  id=acc, description=''))
    return MultipleSeqAlignment(records, generic_protein)


def read_aln(fname, format):
    """Read a sequence alignment."""
    if format:
        assert format.islower()
    if fname.endswith('.cma') or format == 'cma':
        return cma_blocks(fname)
    else:
        aln = AlignIO.read(fname, format)
        # Avoids trouble w/ cogent later on
        for rec in aln:
            rec.description = ''
        return aln


def combine_alignments(fg_aln, bg_aln):
    """Align FG and BG to each other so column numbers match.

    Uses MUSCLE for profile-profile alignment.
    """
    # This would be simpler with NamedTemporaryFile, but Windows doesn't allow
    # multiple open file handles on the same file, so here we are.
    afd, aseqfname = tempfile.mkstemp(text=True)
    os.close(afd)
    bfd, bseqfname = tempfile.mkstemp(text=True)
    os.close(bfd)
    try:
        AlignIO.write(fg_aln, aseqfname, 'fasta')
        AlignIO.write(bg_aln, bseqfname, 'fasta')
        output = subprocess.check_output([
            'muscle', '-profile',
            '-in1', aseqfname,
            '-in2', bseqfname,
        ])
    finally:
        if os.path.exists(aseqfname):
            os.remove(aseqfname)
        if os.path.exists(bseqfname):
            os.remove(bseqfname)

    full_aln = AlignIO.read(StringIO(output), 'fasta')
    full_aln = MultipleSeqAlignment(alnutils.remove_empty_cols(full_aln),
                                    generic_protein)
    # Save a copy
    # ENH: choose a reasonable name
    AlignIO.write(full_aln, '_cc_combined.seq', 'fasta')
    logging.info("Wrote _cc_combined.seq")
    return full_aln


def clean_alignments(fg_aln, bg_aln):
    """Fix simple issues in the alignments:

    - Remove duplicated sequences and IDs from the background
    - Ensure alignments are the same width (if not, align to each other)
    - Remove all-gap columns from the full alignment
    """
    if not len(fg_aln):
        raise ValueError("Foreground set is empty")
    if not len(bg_aln):
        raise ValueError("Background set is empty")

    # Remove FG seqs from BG -- by equal sequence IDs, here
    killme = []
    for fg_seq in fg_aln:
        for idx, bg_seq in enumerate(bg_aln):
            if fg_seq.id == bg_seq.id != 'consensus':
                if str(fg_seq.seq) != str(bg_seq.seq):
                    logging.warn("Different sequences for %s in fg, bg",
                            fg_seq.id)
                killme.append(idx)
    if killme:
        logging.info("Removing %d duplicated sequence IDs from the background",
                    len(killme))
        for idx in sorted(killme, reverse=True):
            del bg_aln._records[idx]

    if not len(bg_aln):
        raise ValueError("Background set is a subset of the foreground. "
                         "To fix this, use a background containing sequences "
                         "not in the foreground.")

    # Remove identical sequences from the FG and BG
    def purge_duplicates(aln, seen=set()):
        out_recs = []
        for rec in aln:
            if str(rec.seq) not in seen:
                out_recs.append(rec)
                seen.add(str(rec.seq))
        diff_fg = len(aln) - len(out_recs)
        if diff_fg:
            logging.warn("Purging %d identical sequences from the alignment",
                    diff_fg)
        return out_recs

    fg_aln._records = purge_duplicates(fg_aln)
    bg_aln._records = purge_duplicates(bg_aln)

    # Ensure alignments are the same width
    if len(fg_aln[0]) != len(bg_aln[0]):
        logging.warn("Alignments are not of equal width; fixing with MUSCLE.")
        full_aln = combine_alignments(fg_aln, bg_aln)
    else:
        full_aln = deepcopy(fg_aln)
        full_aln.extend(bg_aln)

    # Remove columns that are all gaps in both fg and bg (full_aln)
    seqstrs = [str(rec.seq) for rec in full_aln]
    clean_cols = [col for col in zip(*seqstrs)
                  if not all(c == '-' for c in col)]
    clean_seqs = [''.join(row) for row in zip(*clean_cols)]
    for rec, clean_seq in zip(full_aln, clean_seqs):
        rec.seq = Seq(clean_seq, rec.seq.alphabet)

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
    hit_quads = [(i+1, faa_baa_pval[0], faa_baa_pval[1], faa_baa_pval[2])
                 for i, faa_baa_pval in enumerate(hits)]
    get_pval = lambda ifbp: ifbp[3]
    # Benjamini-Hochberg multiple-testing FDR correction (BH step-up)
    hit_quads.sort(key=get_pval)
    m = len(hit_quads)  # Num. hypotheses tested
    if m < N:
        N = m
    tophits = [(posn, faa, baa) for posn, faa, baa, pval in hit_quads]
    for k, ifbp in zip(range(m, 0, -1), reversed(hit_quads))[-N:]:
        # logging.info("BH: a=%s, m=%s, k=%s, p=%s, compare=%s",
        #              alpha, m, k, get_pval(ifbp), alpha * k / m)
        if get_pval(ifbp) <= alpha * k / m:
            return tophits[:k]
    return []



# --- PDB alignment magic ---

# ENH: handle multiple PDBs
@contextlib.contextmanager
def read_pdb_seq(pdb_fname):
    """Return the name of a temporary file containing the PDB atom sequence(s),
    a list of the sequences themselves (as SeqRecords), and a same-length list
    of tuples of (chain ID, chain start resnum, chain end resnum).

    Context manager, so temporary file is automatically removed.
    """
    pdbseqs = list(SeqIO.parse(pdb_fname, 'pdb-atom'))
    try:
        _fd, pdbseqfname = tempfile.mkstemp()
        SeqIO.write(pdbseqs, pdbseqfname, 'fasta')
        yield pdbseqfname, pdbseqs
    finally:
        os.close(_fd)
        if os.path.exists(pdbseqfname):
            os.remove(pdbseqfname)


def choose_best_aligned(aligned):
    """Choose the longest profile match as the "reference" chain.

    Returns a tuple: (sequence ID, sequence string)
    """
    def aligned_len(seq):
        return sum(c.isupper() for c in seq.replace('X', ''))

    if not aligned:
        raise RuntimeError("No PDB sequences were aligned by the profile!")
    elif len(aligned) == 1:
        ref_id, ref_aln = aligned.items()[0]
    else:
        ref_id = max(aligned.iteritems(), key=lambda kv: aligned_len(kv[1]))
        ref_aln = aligned[ref_id]
    return ref_id, ref_aln


def get_aln_offset(full, aln):
    aln_trim = aln.replace('-', '').replace('.', '').upper()
    if 'X' in aln_trim:
        aln_trim = aln_trim[:aln_trim.index('X')]
    return full.index(aln_trim)


def aln_resnums_inserts(record, aln, offset):
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
        elif c == '-':
            # Deletion position
            aln_resnums.append(('-', None))
            del_len += 1
        else:
            raise ValueError("Unexpected character '%s'" % c)
    return aln_resnums, aln_inserts


def pdb_hmm(hmm_profile, pdb_fname):
    """Align a PDB structure to an HMM profile.

    Returns a tuple: (SeqRecord,
                    //chain ID,
                      list of aligned residue numbers,
                      list of insert ranges as tuple pairs)
    """
    with read_pdb_seq(pdb_fname) as (seqfname, seqs):
        out = subprocess.check_output(['hmmalign', '--allcol', '--trim',
                                       '--amino', '--outformat', 'a2m',
                                       hmm_profile, seqfname])
    ref_id, ref_aln = choose_best_aligned(
        dict((rec.id, str(rec.seq))
             for rec in SeqIO.parse(StringIO(out), 'fasta')))
    ref_record = SeqIO.to_dict(seqs)[ref_id]
    # Calculate aligned residue numbers & insert ranges
    offset = (ref_record.annotations['start']
              + get_aln_offset(str(ref_record.seq), ref_aln)
              - 1)
    resnums, inserts = aln_resnums_inserts(ref_record, ref_aln, offset)
    return ref_record, resnums, inserts


def pdb_mapgaps(mapgaps_profile, pdb_fname):
    """Align a PDB structure to a MAPGAPS profile.

    Returns a tuple: (SeqRecord, list of aligned residue numbers)
    """
    from biocma import cma

    with read_pdb_seq(pdb_fname) as (seqfname, seqs):
        subprocess.check_call(['run_gaps', mapgaps_profile, seqfname])
    pdb_cma = cma.read(seqfname + '_aln.cma')
    hits = {}
    head_lengths = {}
    for seq in pdb_cma['sequences']:
        hits[seq['id']] = seq['seq']
        head_lengths[seq['id']] = seq['head_len']
    ref_id, ref_aln = choose_best_aligned(hits)
    ref_record = SeqIO.to_dict(seqs)[ref_id]
    offset = (ref_record.annotations['start']
              + head_lengths[ref_id]
              - 1)
    resnums, inserts = aln_resnums_inserts(ref_record, ref_aln, offset)
    return ref_record, resnums, inserts

