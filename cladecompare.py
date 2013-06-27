#!/usr/bin/env python

"""Compare protein sequence alignments to identify contrasting sites."""

import logging
import sys

from cladecomparelib.core import process_args


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__,
                                 epilog="""Examples:

# Basic
cladecompare.py fg.seq bg.seq > fg-v-bg.out
cladecompare.py fg.cma bg.cma > fg-v-bg.out

# Save per-site p-values and "pattern" file of significant sites
cladecompare.py fg.seq bg.seq -p fg-v-bg.pttrn -o fg-v-bg.out

# Align sequences on the fly by giving the MAPGAPS profile
cladecompare.py fg.fasta bg.fasta -s urn --mapgaps /share/data/PK \\
        -p fg.pttrn -o fg.out

# HMMer/PDB style
cladecompare.py fg.fasta bg.fasta --hmm Pkinase.hmm --pdb 1ATP.pdb

See README or http://github.com/etal/cladecompare for full documentation.
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    # Input
    AP.add_argument('foreground',
            help="Foreground sequence (alignment) file.")
    AP.add_argument('background',
            nargs='*',
            help="""Background sequence (alignment) file, or additional
            foreground alignments for each-vs-all comparisons.""")
    AP.add_argument('-f', '--format',
            default='fasta',
            help="Input alignment format (default: fasta).")
    AP.add_argument('--hmm',
            help="""Align the foreground and background sequences with this
            HMMer 3.0 profile.""")
    AP.add_argument('--mapgaps',
            help="""Align the foreground and background sequences with this
            MAPGAPS profile.""")
    AP.add_argument('--pdb',
            action='append',
            help="""3D structure coordinates in Protein Data Bank format.""")
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
            'phospho' = urn model for phosphorylatable residues (S/T/Y);
            'jsd' = Jensen-Shannon divergence.
            """)
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

