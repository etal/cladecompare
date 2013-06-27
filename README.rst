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

A proper installation looks like::

    python setup.py install

If you have setuptools installed, the dependencies on Biopython_, BioFrills_,
SciPy_ and ReportLab_ will be fetched and installed automatically.

.. _Biopython: http://biopython.org/wiki/Download
.. _biofrills: https://github.com/etal/biofrills
.. _SciPy: http://scipy.org/
.. _ReportLab: http://pypi.python.org/pypi/reportlab

Optionally, CladeCompare can align your sequences for you if you have MUSCLE_,
HMMer3_ or MAPGAPS_ installed.

.. _MUSCLE: http://www.drive5.com/muscle/
.. _HMMer3: http://hmmer.janelia.org/
.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/


If you have the dependencies installed, you can use this package in-place,
without installing it. Just download the source code (git clone, or download
the ZIP file and unpack it) and include the top-level directory in your system
path, or add symbolic links to cladecompare.py, cladereport.py and cladeweb.py
to an existing directory that's in your path (e.g. ``~/bin``).


Testing
~~~~~~~

Finally, if you are on a Unix-like system (i.e. Linux, Mac or Cygwin), you can
verify your installation by running the test suite. Change to the ``test/``
directory and run ``make``::

    cd test
    make

If CladeCompare is installed correctly, the program will run in several modes
and generate output files. View the ``.html`` files in your web browser to see
what happened.


Usage
-----

Web interface
~~~~~~~~~~~~~

Launch the script ``cladeweb.py`` and fill in the form in your web browser.
The form accepts sequences in FASTA or CMA format, and you can upload an HMM
profile to align unaligned FASTA sequence sets. (See below for details about
each field.)

If you launched the application from the command line, press Ctrl-C (on
Unix-like systems) to stop the web server application.

Note that only one instance of the server will run on your system at a time; if
you launch ``cladeweb.py`` twice in a row, another browser tab or window will
open but the server will not restart.


Command line
~~~~~~~~~~~~

The command-line interface ``cladecompare.py`` provides the same functionality
as the web interface, plus a few more options.  To read the built-in help and
option descriptions::

    cladecompare.py --help

Two alignments are compared by specifying the foreground and background sets,
in that order, as arguments::

    # Compare two alignments
    cladecompare.py fg_aln.seq bg_aln.seq

The program prints the following information for each column in the alignment(s):

- The consensus amino acid types of the foreground and background
- p-value indicating the significance of the contrast in amino acid frequencies
- a little ASCII bar chart indicating contrast, based on the p-value.

P-values are adjusted for number of columns in the alignment with the
Benjamini-Hochberg "step-up" multiple-testing correction (false discovery rate,
FDR).

Redirect the output to a file with the extension ".out"::

    # Compare two alignments
    cladecompare.py fg_aln.seq bg_aln.seq > fg-v-bg.out

Or specify the output file name with the ``-o`` option (same effect)::

    cladecompare.py fg_aln.seq bg_aln.seq -o fg-v-bg.out

If you're not using MAPGAPS_, it would make sense to either:

- Create a sequence alignment of all sequences, foreground and background,
  together; then divide the alignment into two FASTA files (e.g. by sequence
  ID).
- Align both the foreground and background sets with hmmalign (HMMer3_) using
  the same profile, then use a script (your own) to delete any insert columns.

In case you botch all that, CladeCompare will check then number of columns in
the FG and BG alignments, and if they don't match, will automatically run MUSCLE
to align the alignments to each other.

To align the FASTA sequences with MAPGAPS on the fly, specify the profile name
(minus the extension) with the ``--mapgaps`` flag::

    # Align sequence sets w/ MAPGAPS, then compare
    cladecompare.py test/scttlre-domain.fasta test/cdc2-domain.fasta --mapgaps test/CDK_CMGC \
        -o scttlre-v-cdc2.out

Pre-aligned sequences in CMA format (.cma) are also accepted::

    # Use pre-aligned sequence sets (MAPGAPS "CMA" format)
    cladecompare.py test/scttlre-domain.fasta_aln.cma test/cdc2-domain.fasta_aln.cma \
        -o scttlre-v-cdc2.out

Finally, given the '-p' option, cladecompare.py will write a "pattern" file
listing the alignment column numbers with significant contrasts, in decreasing
order (this can be useful input to other scripts of your own), as well as PDF
files of paired sequence logos representing the foreground and background
alignments around each significant site::

    # Specify where the outputs go
    cladecompare.py fg_aln.seq bg_aln.seq -o fg-v-bg.out -p fg-v-bg.pttrn

Outputs
```````

The script ``cladereport.py`` converts the "\*.out" files to an HTML file showing
the alignment of the FG and BG consensus sequences, with the FG sequence
colorized to show per-site contrasts (red=significant difference,
blue=non-significant/columns are similar), inserts (yellow) and deletions (gray
gaps)::

    # Visualize the per-site contrasts as a colorized alignment
    cladereport.py scttlre-v-cdc2.out > scttlre-v-cdc2.html

Single- and multi-profile modes
```````````````````````````````

If a single sequence set is given, the aligned columns are compared to the
overall amino-acid frequencies of the alignment::

    cladecompare.py subfam1.seq -o subfam1-single.out

When more than 2 sequence sets are given, each set is individually treated as a
foreground and the rest treated as the background for evaluation::

    # Compare several related alignments, e.g. all subfamilies
    cladecompare.py subfam1.seq subfam2.seq subfam3.seq ...

This multi-mode generates and names the "\*.out" files according to the
corresponding sequence file names. You can visualize these all together::

    # Visualize each subfamily's contrasts together
    cladereport.py subfam2.out subfam2.out subfam3.out ... > somefamily.html


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
:jsd:
    Jensen-Shannon divergence of column compositions, a la INTREPID_.

.. _CHAIN: http://chain.igs.umaryland.edu/
.. _INTREPID: http://bioinformatics.oxfordjournals.org/content/24/21/2445.full

