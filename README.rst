============
CladeCompare
============

Compare protein sequence alignments. Identify diagnostic residues between the
given "foreground" (FG) and  "background" (BG) clades.

If you use this software in a publication, please cite our paper that describes
it:

    Talevich, E. & Kannan, N. (2013) Structural and evolutionary adaptation of
    rhoptry kinases and pseudokinases, a family of coccidian virulence factors.
    *BMC Evolutionary Biology* 13:117
    doi:10.1186/1471-2148-13-117

    Available at: http://www.biomedcentral.com/1471-2148/13/117


Freely distributed under the permissive BSD 2-clause license (see LICENSE).

Installation
------------

CladeCompare requires Python 3.12 or later. It depends on Biopython_,
BioFrills_ (with its ``logo`` extra, which brings in WebLogo), SciPy_,
ReportLab_ and bottle_, and installs three commands: ``cladecompare``,
``cladereport`` and ``cladeweb``.

Until BioFrills 0.4.0 is released on PyPI, install CladeCompare from a copy of
this repository. With uv_, the dependencies are resolved automatically,
BioFrills included, from the commit recorded in ``pyproject.toml``::

    uv tool install .

With pip, install that BioFrills commit first, then CladeCompare::

    pip install "biofrills[logo] @ git+https://github.com/etal/biofrills.git@ff11246a93e14a35723e6f04a3952c88c0cbf29c"
    pip install .

.. _Biopython: https://biopython.org/
.. _biofrills: https://github.com/etal/biofrills
.. _SciPy: https://scipy.org/
.. _ReportLab: https://pypi.org/project/reportlab/
.. _bottle: https://bottlepy.org/
.. _uv: https://docs.astral.sh/uv/

Optionally, CladeCompare can align your sequences for you if you have HMMER3_
(3.0 or later) or MAPGAPS_ installed. If a foreground and a
background alignment differ in width, CladeCompare aligns them to each other
with MUSCLE_, which must be MUSCLE 3, such as 3.8.1551: MUSCLE 5 has no
profile-profile alignment mode, and CladeCompare stops with an error if it
finds MUSCLE 4 or later.

.. _MUSCLE: http://www.drive5.com/muscle/
.. _HMMER3: http://hmmer.org/
.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/

The pixi_ environment described under Development below provides HMMER and
MUSCLE 3.8.1551 from Bioconda, together with CladeCompare itself.

.. _pixi: https://pixi.sh/


Usage
-----

Web interface
~~~~~~~~~~~~~

Launch the command ``cladeweb`` and fill in the form in your web browser.
The form accepts sequences in FASTA or CMA format, and you can upload an HMM
profile to align unaligned FASTA sequence sets. (See below for details about
each field.) Sequences are not weighted unless you check the weighting box,
which corresponds to the ``-w`` option of ``cladecompare``.

If you launched the application from the command line, press Ctrl-C (on
Unix-like systems) to stop the web server application.

Note that only one instance of the server will run on your system at a time; if
you launch ``cladeweb`` twice in a row, another browser tab or window will
open but the server will not restart.


Command line
~~~~~~~~~~~~

The command-line interface ``cladecompare`` provides the same functionality
as the web interface, plus a few more options.  To read the built-in help and
option descriptions::

    cladecompare --help

Two alignments are compared by specifying the foreground and background sets,
in that order, as arguments::

    # Compare two alignments
    cladecompare fg_aln.seq bg_aln.seq

The program prints the following information for each column in the alignment(s):

- The consensus amino acid types of the foreground and background
- p-value indicating the significance of the contrast in amino acid frequencies
- an ASCII bar chart indicating contrast, based on the p-value: one asterisk
  per order of magnitude, shown only where the p-value is below the
  significance threshold alpha (the ``-a`` option, between 0 and 1, default
  0.005).

These per-column p-values are not adjusted for multiple testing, so with
hundreds of columns some may fall below alpha by chance. The multiple-testing
correction applies to the pattern file described below.

Redirect the output to a file with the extension ".out"::

    # Compare two alignments
    cladecompare fg_aln.seq bg_aln.seq > fg-v-bg.out

Or specify the output file name with the ``-o`` option (same effect)::

    cladecompare fg_aln.seq bg_aln.seq -o fg-v-bg.out

If you're not using MAPGAPS_, it would make sense to either:

- Create a sequence alignment of all sequences, foreground and background,
  together; then divide the alignment into two FASTA files (e.g. by sequence
  ID).
- Align both the foreground and background sets with hmmalign (HMMER3_) using
  the same profile, then use a script (your own) to delete any insert columns.

If you do neither, and the foreground and background alignments differ in
width, CladeCompare aligns them to each other with MUSCLE.

To align the FASTA sequences with MAPGAPS on the fly, specify the profile name
(minus the extension) with the ``--mapgaps`` flag::

    # Align sequence sets w/ MAPGAPS, then compare
    cladecompare test/scttlre-domain.fasta test/cdc2-domain.fasta --mapgaps test/CDK_CMGC \
        -o scttlre-v-cdc2.out

Pre-aligned sequences in CMA format (.cma) are also accepted::

    # Use pre-aligned sequence sets (MAPGAPS "CMA" format)
    cladecompare test/scttlre-domain.fasta_aln.cma test/cdc2-domain.fasta_aln.cma \
        -o scttlre-v-cdc2.out

Finally, given the '-p' option, cladecompare will write a "pattern" file
listing the alignment column numbers with significant contrasts, most
significant first (this can be useful input to other scripts of your own), as
well as PDF files of paired sequence logos representing the foreground and
background alignments around each significant site. Only these outputs, and the
residues highlighted in the PyMOL script written with ``--pdb``, are corrected
for multiple testing: a column is significant if it passes the
Benjamini-Hochberg "step-up" procedure at a false discovery rate (FDR) of alpha
(the ``-a`` option), taken over all columns of the alignment. At most the 50
most significant columns are considered::

    # Specify where the outputs go
    cladecompare fg_aln.seq bg_aln.seq -o fg-v-bg.out -p fg-v-bg.pttrn

Outputs
```````

The command ``cladereport`` converts the "\*.out" files to an HTML file showing
the alignment of the FG and BG consensus sequences, with the FG sequence
colorized by the uncorrected p-value of each site (red for p below 0.1, darker
as p decreases; white from 0.1 to 0.5; blue for larger p-values, i.e. similar
columns), inserts (yellow) and deletions (gray gaps). The colors do not depend
on alpha::

    # Visualize the per-site contrasts as a colorized alignment
    cladereport scttlre-v-cdc2.out > scttlre-v-cdc2.html

Single- and multi-profile modes
```````````````````````````````

If a single sequence set is given, the aligned columns are compared to the
overall amino-acid frequencies of the alignment::

    cladecompare subfam1.seq -o subfam1-single.out

When more than 2 sequence sets are given, each set is individually treated as a
foreground and the rest treated as the background for evaluation::

    # Compare several related alignments, e.g. all subfamilies
    cladecompare subfam1.seq subfam2.seq subfam3.seq ...

This multi-mode generates and names the "\*.out" files according to the
corresponding sequence file names; the ``-o`` and ``-p`` options are ignored,
and ``cladecompare`` warns if either is given. You can visualize these all
together::

    # Visualize each subfamily's contrasts together
    cladereport subfam1.out subfam2.out subfam3.out ... > somefamily.html


Strategies
----------

Statistical tests ("-s" options) for column comparison:

:gtest:
    (default) G-test for goodness-of-fit of FG amino acid counts vs. those of
    the BG column. BG frequencies include pseudocounts calculated from the
    amino acid frequencies of the full sequence set.
:urn:
    Ball-in-urn model (binomial), a la CHAIN_, for counts of the "consensus"
    amino acid type in FG and BG.
:phospho:
    Ball-in-urn model, as in "urn", for the combined count of the
    phosphorylatable residues Ser, Thr and Tyr in the FG column, against
    their combined frequency in the FG and BG columns together.
:jsd:
    Jensen-Shannon divergence of column compositions, a la INTREPID_.
:hypg:
    Hypergeometric model: the probability of drawing at least as many
    sequences of the FG consensus amino acid type as the FG column holds,
    when as many sequences as the FG holds are drawn from the FG and BG
    together. With a single alignment, the fraction of the column that
    matches the consensus, which is not a p-value.

.. _CHAIN: http://chain.igs.umaryland.edu/
.. _INTREPID: http://bioinformatics.oxfordjournals.org/content/24/21/2445.full


Development
-----------

Two development workflows are supported from the same ``pyproject.toml``, and
both lock files are committed. Both install CladeCompare in editable mode with
the tools in the ``dev`` dependency group, pytest and ruff. Both also need a
``git`` executable, because BioFrills is installed from its repository until
it is released on PyPI.

uv_ takes the Python packages from PyPI, apart from BioFrills, at the versions
recorded in ``uv.lock``::

    uv sync
    uv run pytest
    uv run ruff check
    uv run ruff format

pixi_ takes Python and most Python packages from conda-forge, and also
installs HMMER and MUSCLE 3.8.1551 from Bioconda, at the versions recorded in
``pixi.lock``. The ``default`` environment uses the latest Python that
conda-forge provides, and ``py312`` the oldest supported one::

    pixi run test
    pixi run -e py312 test
    pixi run test-integration
    pixi run lint            # ruff check
    pixi run fmt             # ruff format

The tests run ``cladecompare``, ``cladereport`` and the ``cladeweb``
application on the inputs in ``tests/data/inputs`` and compare what they
write with the reference outputs in ``tests/data/py27``, apart from the
deliberate differences that ``tests/data/README.md`` lists; unit tests check
the parsing, site-selection and statistics functions against hand
computations. The tests marked ``integration`` also align the inputs with
``hmmalign`` or ``muscle``. Each is skipped unless the program it needs is on
``PATH``; ``pixi run test-integration`` runs them alone and fails instead of
skipping when a program is missing.

The reference outputs come from the Python 2 code base. On macOS, the
``py27-oracle`` environment reproduces the environment that captured them,
and ``pixi run -e py27-oracle capture`` regenerates them, as described in
``tests/data/README.md``.

Code is formatted with ``ruff format``, and the formatting-only commits are
listed in ``.git-blame-ignore-revs``. To skip them in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

