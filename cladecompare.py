#!/usr/bin/env python

"""Compare alignments to identify distinguishing/diagnostic residues.

Examples:

# Basic
cladecompare.py fg.seq bg.seq > fg-v-bg.out

# Save per-site p-values and "pattern" file of significant sites
cladecompare.py fg.seq bg.seq -p fg-v-bg.pttrn -o fg-v-bg.out

# CHAIN style
# Align sequences on the fly by giving the MAPGAPS profile
cladecompare.py fg.fasta bg.fasta -s urn --mapgaps /share/data/PK \\
        -p fg.pttrn -o fg.out

# mcBPPS style
# (Run MAPGAPS on each subfamily beforehand)
# Writes files named: subfam{1,2,3}.fa_aln.cma.{pttrn,out}
cladecompare.py subfam1.fa_aln.cma subfam2.fa_aln.cma subfam3.fa_aln.cma

"""

from __future__ import absolute_import

import logging
import math
import os
import subprocess
import sys
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

from biofrills import alnutils

from cladecompare import pairlogo, urn, gtest, jsd


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


# --- Output ---

def write_report(hits, outfile, alpha):
    """Write p-values & "contrast" stars for each site. (It's noisy.)"""
    for idx, data in enumerate(hits):
        fg_char, bg_char, pvalue = data
        if not (0 <= pvalue <= 1):
            logging.warn("Out-of-domain p-value at site %s: %s",
                         idx, pvalue)
        stars = ('*'*int(-math.log10(pvalue)) if 0 < pvalue < alpha else '')
        outfile.write("%s (%s) %d : prob=%g\t%s\n"
                      % (fg_char, bg_char, idx + 1, pvalue, stars))


def write_mcbpps(tophits, ptnfile):
    """Write a .pttrn file in the style of mcBPPS."""
    ptnfile.write("1:" + ','.join([("%s%d" % (faa, posn))
                                   for posn, faa, baa in tophits]))


# ---- FLOW --------------------------------------------------------------

# ENH:
#   single alignment
#       compare_aln main.cma [--tree foo]
#       --> longest-branch decomposition
#       * stop at some threshold -- clade size 3, 5, 6?
def process_args(args):
    """Main."""
    if args.mapgaps:
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

    if args.strategy == 'urn':
        logging.info("Using ball-in-urn statistical model")
    elif args.strategy == 'gtest':
        logging.info("Using G-test of amino acid frequencies")
    elif args.strategy == 'jsd':
        logging.info("Using Jensen-Shannon divergence")
    elif args.strategy == 'ancestrallrt':
        logging.info("Using maximum-likelihood ancestral character states")
    else:
        raise ValueError("Unknown strategy: %s" % args.strategy)

    if len(all_alns) == 2:
        fg_aln, bg_aln = all_alns
        fg_clean, bg_clean, hits = process_pair(fg_aln, bg_aln,
                                                args.strategy, args.tree)
        process_output(fg_clean, bg_clean, hits, args.alpha,
                       args.output, args.pattern)
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
                                                    args.strategy, args.tree)
            outfname, ptnfname = outfnames_ptnfnames[idx]
            process_output(fg_clean, bg_clean, hits, args.alpha,
                           outfname, ptnfname)
            logging.info("Wrote %s and %s", outfname, ptnfname)


def process_pair(fg_aln, bg_aln, strategy, tree=None):
    """Calculate a mapping of alignment column positions to "contrast".

    Return a list of tuples:
        (foreground consensus aa, background consensus aa, p-value)
        for each column position.
    """
    fg_aln, bg_aln = clean_alignments(fg_aln, bg_aln)
    if strategy == 'urn':
        hits = urn.compare_aln(fg_aln, bg_aln)
    elif strategy == 'gtest':
        hits = gtest.compare_aln(fg_aln, bg_aln)
    elif strategy == 'jsd':
        hits = jsd.compare_aln(fg_aln, bg_aln)
    elif strategy == 'ancestrallrt':
        # hits = ancestrallrt.compare_aln(fg_aln, bg_aln, tree)
        pass
    return fg_aln, bg_aln, hits


def process_output(fg_aln, bg_aln, hits, alpha, output, pattern):
    with as_handle(output, 'w+') as outfile:
        write_report(hits, outfile, alpha)
    if pattern:
        tophits = top_hits(hits, alpha)
        with open(pattern, 'w+') as ptnfile:
            write_mcbpps(tophits, ptnfile)
        pairlogo.make_pairlogos(fg_aln, bg_aln, tophits,
                                pattern.rsplit('.', 1)[0],
                                10)



if __name__ == '__main__':
    import argparse 
    AP = argparse.ArgumentParser(__doc__)
    # Input
    AP.add_argument('foreground',
            help="Foreground sequence (alignment) file.")
    AP.add_argument('background',
            nargs='+',
            help="""Background sequence (alignment) file, or additional
            foreground alignments for each-vs-all comparisons.""")
    AP.add_argument('-f', '--format',
            default='fasta',
            help="Input alignment format (default fasta).")
    AP.add_argument('--hmm',
            help="""Align the foreground and background sequences with this
            HMMer 3.0 profile.""")
    AP.add_argument('--mapgaps',
            help="""Align the foreground and background sequences with this
            MAPGAPS profile.""")
    AP.add_argument('-t', '--tree',
            help="Input tree.")
    # Options
    # add_argument_group?
    AP.add_argument('-a', '--alpha',
            default=0.05, type=float,
            help="Significance threshold for pattern columns.")
    AP.add_argument('-s', '--strategy',
            default='gtest',
            help="""Strategy used to compare alignments:
            'gtest' = G-test of all character frequencies;
            'urn' = ball-in-urn model of consensus residue conservation;
            'jsd' = Jensen-Shannon divergence.
            """)
            # 'ancestrallrt' = likelihood ratio test of ancestral state
            #         likelihoods between the foreground clade and the full tree.
    # Output
    AP.add_argument('-o', '--output',
            default=sys.stdout,
            help="""Write per-column probabilities (standard output) to this
            filename, for use with cladereport.py. (Single-foreground
            comparison only.""")
    AP.add_argument('-p', '--pattern',
            help="""Write an mcBPPS-style pattern file to this filename.
            (Single-foreground comparison only.""")
    AP.add_argument('-q', '--quiet',
            action='store_true',
            help="Don't print status messages, only warnings and errors.")

    args = AP.parse_args()
    if args.quiet:
        logging.basicConfig(level=logging.WARNING,
                format="%(module)s: %(message)s")
    else:
        logging.basicConfig(level=logging.INFO,
                format="%(module)s [@%(lineno)s]: %(message)s")
    process_args(args)

