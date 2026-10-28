"""Compare protein sequence alignments to identify contrasting sites."""

import argparse

from cladecomparelib import InputError, configure_logging
from cladecomparelib.core import (
    DEFAULT_STRATEGY,
    STRATEGIES,
    parse_alpha,
    process_args,
)


def main():
    AP = argparse.ArgumentParser(
        description=__doc__,
        epilog="""Examples:

# Basic
cladecompare fg.seq bg.seq > fg-v-bg.out
cladecompare fg.cma bg.cma > fg-v-bg.out

# Save per-site p-values and "pattern" file of significant sites
cladecompare fg.seq bg.seq -p fg-v-bg.pttrn -o fg-v-bg.out

# Align sequences on the fly by giving the MAPGAPS profile
cladecompare fg.fasta bg.fasta -s urn --mapgaps /share/data/PK \\
        -p fg.pttrn -o fg.out

# HMMER/PDB style
cladecompare fg.fasta bg.fasta --hmm Pkinase.hmm --pdb 1ATP.pdb

See README or https://github.com/etal/cladecompare for full documentation.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Input
    AP.add_argument("foreground", help="Foreground sequence (alignment) file.")
    AP.add_argument(
        "background",
        nargs="*",
        help="""Background sequence (alignment) file, or additional
            foreground alignments for each-vs-all comparisons.""",
    )
    AP.add_argument(
        "-f",
        "--format",
        default="fasta",
        help="Input alignment format (default: fasta).",
    )
    AP.add_argument(
        "--hmm",
        help="""Align the foreground and background sequences with this
            HMMER 3 profile.""",
    )
    AP.add_argument(
        "--mapgaps",
        help="""Align the foreground and background sequences with this
            MAPGAPS profile.""",
    )
    AP.add_argument(
        "--pdb",
        action="append",
        help="""3D structure coordinates in Protein Data Bank format.""",
    )
    # Options
    # add_argument_group?
    AP.add_argument(
        "-a",
        "--alpha",
        default=0.005,
        type=alpha,
        help="Significance threshold for pattern columns, between 0 and 1.",
    )
    AP.add_argument(
        "-s",
        "--strategy",
        default=DEFAULT_STRATEGY,
        choices=STRATEGIES,
        help="Strategy used to compare alignments: "
        + "; ".join(f"'{name}' = {desc}" for name, (desc, _) in STRATEGIES.items())
        + ".",
    )
    AP.add_argument(
        "-w",
        "--weight",
        action="store_true",
        help="""Enable sequence weighting for statistical tests.""",
    )
    # Output
    AP.add_argument(
        "-o",
        "--output",
        help="""Write per-column probabilities to this filename, for use
            with cladereport (default: standard output). Ignored when
            comparing more than two alignments, which writes <alignment>.out
            for each alignment.""",
    )
    AP.add_argument(
        "-p",
        "--pattern",
        help="""Write an mcBPPS-style pattern file to this filename. Ignored
            when comparing more than two alignments, which writes
            <alignment>.pttrn for each alignment.""",
    )
    AP.add_argument(
        "--pmlout", help="Write a PyMOL script to this filename (if --pdb given)."
    )
    AP.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Don't print status messages, only warnings and errors.",
    )

    args = AP.parse_args()
    configure_logging(args.quiet)
    try:
        process_args(args)
    except InputError as exc:
        AP.exit(1, f"{AP.prog}: error: {exc}\n")


def alpha(text):
    """Parse the -a value with core.parse_alpha.

    argparse reports the ValueError of a non-numeric value as an invalid
    alpha value, and the message of an ArgumentTypeError unchanged.
    """
    try:
        return parse_alpha(text)
    except InputError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None
