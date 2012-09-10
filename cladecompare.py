#!/usr/bin/env python

"""Compare alignments to identify distinguishing/diagnostic residues.

Examples:

# CHAIN style
# Align sequences on the fly by giving the MAPGAPS profile
compare_aln.py fg.fasta bg.fasta -s ballinurn --mapgaps /share/data/PK \
        -p fg.pttrn -o fg.noise

# mcBPPS style
# (Run MAPGAPS on each subfamily beforehand)
compare_aln.py subfam1.fa_aln.cma subfam2.fa_aln.cma subfam3.fa_aln.cma
# Writes files named: subfam{1,2,3}.fa_aln.cma.{pttrn,noise}

"""
# Big ENH: process noise (as a script, or build it in here) to show a wrapped
# HTML alignment, where contrast (p-value, #stars) is shown by color, e.g.
# redness

from __future__ import absolute_import

import logging
import math
import subprocess
import sys
from cStringIO import StringIO
from copy import deepcopy
from os.path import basename

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.File import as_handle

from esbglib.sugar import log_config
log_config()

from cladecompare import pairlogo, ballinurn, gtest #, ancestrallrt

# --- Input magic ---

def mapgaps_align_and_read(mapgaps_profile, fasta_fname):
    """Align a FASTA file with MAPGAPS and read the CMA alignment.

    """
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
                acc = line[1:].split(None, 1)[0]
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


def clean_alignments(fg_aln, bg_aln):
    """Fix simple issues in the alignments:

    - Remove duplicated sequences and IDs from the background
    - Ensure alignments are the same width (if not, realign w/ MAFFT)
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
        logging.warn("Alignments are not of equal width; fixing with MAFFT.")
        # Attempt meta-alignment w/ mafft
        AlignIO.write(fg_aln, '_tmp_fg.seq', 'fasta')
        AlignIO.write(bg_aln, '_tmp_bg.seq', 'fasta')
        # Empty/dummy file for MAFFT -- we're really only aligning seeds
        subprocess.check_call('echo > _tmp_mt.seq', shell=True)
        output = subprocess.check_output(['mafft',
            '--globalgenafpair',
            '--maxiterate', '1000',
            '--seed', '_tmp_fg.seq',
            '--seed', '_tmp_bg.seq',
            '_tmp_mt.seq'])
        full_aln = AlignIO.read(StringIO(output), 'fasta')
        # Remove "_seed_" prefixes
        for rec in full_aln:
            rec.id = rec.id[len('_seed_'):]
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


# --- Output ---

def write_noise(hits, outfile, alpha):
    """Write p-values & "contrast" stars for each site. (It's noisy.)"""
    for idx, data in enumerate(hits):
        fg_char, bg_char, pvalue = data
        outfile.write("%s (%s) %d : prob=%1.3g\t%s\n"
                      % (fg_char, bg_char, idx + 1, pvalue,
                         ('*'*int(-math.log(pvalue, 10))
                          if pvalue < alpha else '')))


def write_mcbpps(hits, ptnfile, alpha):
    """Write a .pttrn file in the style of mcBPPS."""
    hit_triples = [(i+1, faa_baa_pval[0], faa_baa_pval[2])
                   for i, faa_baa_pval in enumerate(hits)]
    hit_triples.sort(key=lambda icp: icp[2])
    ptnfile.write("1:" +
                    ','.join([("%s%d" % (aa, posn))
                            for posn, aa, pval in hit_triples
                            # Bonferroni correction (crude)
                            if pval <= alpha #/ len(hits)
                            # Take no more than the top 20 hits
                            ][:20]))

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
        if args.mapgaps:
            aln = mapgaps_align_and_read(args.mapgaps, alnfname)
        else:
            aln = read_aln(alnfname, args.format)
        all_alns.append(aln)

    if args.strategy == 'ballinurn':
        logging.info("Using ball-in-urn statistical model")
    elif args.strategy == 'gtest':
        logging.info("Using G-test of amino acid frequencies")
    elif args.strategy == 'ancestrallrt':
        logging.info("Using maximum-likelihood ancestral character states")
    else:
        raise ValueError("Unknown strategy: %s" % args.strategy)

    if len(all_alns) == 2:
        fg_aln, bg_aln = all_alns
        hits = process_pair(fg_aln, bg_aln, args.strategy, args.tree)
        process_output(hits, args.alpha, args.output, args.pattern)
    else:
        # Output fnames are based on fg filenames; ignore what's given
        outfnames_ptnfnames = [(basename(alnfname) + '.noise',
                                basename(alnfname) + '.pttrn')
                               for alnfname in ([args.foreground] +
                                                args.background)]
        for idx, fg_aln in enumerate(all_alns):
            # Combine other alns into bg
            _other_alns = all_alns[:idx] + all_alns[idx+1:]
            bg_aln = deepcopy(_other_alns[0])
            for otra in _other_alns[1:]:
                bg_aln.extend(deepcopy(otra))
            hits = process_pair(deepcopy(fg_aln), bg_aln,
                                args.strategy, args.tree)
            outfname, ptnfname = outfnames_ptnfnames[idx]
            process_output(hits, args.alpha, outfname, ptnfname)
            logging.info("Wrote %s and %s", outfname, ptnfname)


def process_pair(fg_aln, bg_aln, strategy, tree=None):
    """Calculate a mapping of alignment column positions to "contrast".

    Return a list of tuples:
        (foreground consensus aa, background consensus aa, p-value)
        for each column position.
    """
    fg_aln, bg_aln = clean_alignments(fg_aln, bg_aln)
    if strategy == 'ballinurn':
        # CHAIN-style comparison
        hits = ballinurn.compare_aln(fg_aln, bg_aln)
    elif strategy == 'gtest':
        # Likelihood-based character frequency comparison
        hits = gtest.compare_aln(fg_aln, bg_aln)
    elif strategy == 'ancestrallrt':
        # LRT of ancestral character state likelihoods
        # hits = ancestrallrt.compare_aln(fg_aln, bg_aln, tree)
        pass
    return fg_aln, bg_aln, hits


def process_output(hits, alpha, output, pattern):
    # ENH: Benjamini-Hochberg multiple-testing correction to adjust
    # alpha/p-values for alignment size (#cols), to identify truly significant
    # sites
    with as_handle(output, 'w+') as outfile:
        write_noise(hits, outfile, alpha)
    if pattern:
        with open(pattern, 'w+') as ptnfile:
            write_mcbpps(hits, ptnfile, alpha)



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
    AP.add_argument('-m', '--mapgaps',
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
            # TODO enforce one of: 'gtest', 'ballinurn', 'ancestrallrt'
            default='gtest',
            help="""Strategy used to compare alignments:
            'ballinurn' = ball-in-urn comparison method (like CHAIN);
            'gtest' = G-test of all character frequencies;
            'ancestrallrt' = likelihood ratio test of ancestral state
                    likelihoods between the foreground clade and the full tree.
            """)
    # Output
    AP.add_argument('-p', '--pattern',
            help="""Write an mcBPPS-style pattern file to this filename.
            (Single-foreground comparison only.""")
    AP.add_argument('-o', '--output',
            default=sys.stdout,
            help="""Write per-column probabilities ('noise') to this filename.
            (Single-foreground comparison only.""")

    process_args(AP.parse_args())

