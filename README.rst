============
CladeCompare
============

Compare protein sequence alignments. Identify diagnostic residues between the
given "foreground" (FG) and  "background" (BG) clades.

License: BSD

Installation
------------

You can probably use this package in-place, without installing it.

A proper installation looks like::

    python setup.py install

If you have setuptools installed, the dependencies on Biopython_, BioFrills_,
SciPy_ and ReportLab_ will be fetched and installed automatically.

.. _Biopython: http://biopython.org/wiki/Download
.. _biofrills: https://github.com/etal/biofrills
.. _SciPy: http://scipy.org/
.. _ReportLab: http://pypi.python.org/pypi/reportlab


Usage
-----

::

    # Compare two alignments
    cladecompare.py fg_aln.seq bg_aln.seq

This prints the following information for each column in the alignment(s):

- The consensus amino acid types of the foreground and background
- p-value indicating the significance of the contrast in amino acid frequencies
- a little ASCII bar chart indicating contrast, based on the p-value.

P-values are adjusted for number of columns in the alignment with the
Benjamini-Hochberg "step-up" multiple-testing correction (false discovery rate,
FDR).

I redirect the output to a file with the extension ".noise"::

    # Compare two alignments
    cladecompare.py fg_aln.seq bg_aln.seq > fg-v-bg.noise

If you're not using MAPGAPS_, it would make sense to either:

- Create a sequence alignment of all sequences, foreground and background,
  together; then divide the alignment into two FASTA files (e.g. by sequence
  ID).
- Align both the foreground and background sets with HMMer_ using the same
  profile, then use a script (your own) to delete any insert columns.

In case you botch all that, CladeCompare will check then number of columns in
the FG and BG alignments, and if they don't match, will automatically run MAFFT
to align the alignments to each other.

If you are using MAPGAPS, it looks like this::

    # Align sequence sets w/ MAPGAPS, then compare
    cladecompare.py fg_raw.fasta bg_raw.fasta -m ~/mapgaps_profiles/SomeFamily

    # Use pre-aligned sequence sets (MAPGAPS "CMA" format)
    cladecompare.py fg_raw.fasta_aln.cma bg_raw.fasta_aln.cma

The script noise2html.py converts the "\*.noise" files to an HTML file showing
the alignment of the FG and BG consensus sequences, with the FG sequence
colorized to show per-site contrasts (red=significant difference,
blue=non-significant/columns are similar), inserts (yellow) and deletions (gray
gaps)::

    # Visualize the per-site contrasts as a colorized alignment
    noise2html.py fg-v-bg.noise > fg-v-bg.html

When more than 2 sequence sets are given, each set is individually treated as a
foreground and the rest treated as the background for evaluation::

    # Compare several related alignments, e.g. all subfamilies
    cladecompare.py subfam1.seq subfam2.seq subfam3.seq ...

This multi-mode generates and names the "\*.noise" files according to the
corresponding sequence file names. You can visualize these all together::

    # Visualize each subfamily's contrasts together
    noise2html.py subfam2.noise subfam2.noise subfam3.noise ... > somefamily.html

Finally, given the '-p' option, cladecompare.py will write a "pattern" file
listing the alignment column numbers with significant contrasts, in decreasing
order (this can be useful input to other scripts of your own)::

    # Specify where the outputs go
    cladecompare.py fg_aln.seq bg_aln.seq -o fg-v-bg.noise -p fg-v-bg.pttrn

To read the built-in help and detailed options::

    cladecompare.py --help


.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/
.. _HMMer: http://hmmer.janelia.org/

Strategies
----------

Statistical tests ("-s" options) for column comparison:

:gtest:
    (default) G-test for goodness-of-fit of FG amino acid counts vs. those of
    the BG column. BG frequencies include pseudocounts calculated from the
    amino acid frequencies of the full sequence set.
:ballinurn:
    Ball-in-urn model (binomial), a la CHAIN_, for counts of the "consensus"
    amino acid type in FG and BG.
:ancestrallrt:
    (in progress) Likelihood ratio test of ancestral states
:entropy:
    (in progress) Difference in column entropies, i.e. information gain/loss

.. _CHAIN: http://chain.igs.umaryland.edu/

