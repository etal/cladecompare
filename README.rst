============
CladeCompare
============

License: BSD

Compare alignments. Identify diagnostic residues between clades.

Usage::

    # Compare two alignments
    cladecompare.py fg_aln.seq bg_aln.seq

    # Align w/ MAPGAPGS and then compare two sequence sets
    cladecompare.py fg_raw.fasta bg_raw.fasta -m ~/db/PK

Modes:

    - Ball-in-urn model, a la CHAIN, for counts of the "consensus" amino acid
      type in FG and BG
    - G-test for goodness-of-fit of amino acid counts in FG vs. BG columns

