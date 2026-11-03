============
CladeCompare
============

CladeCompare compares two protein sequence alignments, a *foreground* (FG)
clade and a *background* (BG) of related sequences, and reports for each
alignment column how strongly the amino acid composition of the foreground
differs from that of the background. Its purpose is to identify the residues
that distinguish a protein subfamily from its relatives.

If you use this software in a publication, please cite the paper that
describes it:

    Talevich E, Kannan N (2013) Structural and evolutionary adaptation of
    rhoptry kinases and pseudokinases, a family of coccidian virulence
    factors. *BMC Evolutionary Biology* 13:117.
    doi:`10.1186/1471-2148-13-117 <https://doi.org/10.1186/1471-2148-13-117>`_

CladeCompare is distributed under the BSD 2-clause license (see
``LICENSE``). The font it ships for drawing sequence logos, DejaVu Sans Mono,
is under the Bitstream Vera license (see
``cladecomparelib/data/LICENSE.DejaVu``).

What CladeCompare is for
------------------------

The question
    Given a protein family, which positions of the alignment distinguish
    one clade of the family, the foreground, from related sequences, the
    background? A position where the foreground holds residues that the
    background rarely holds, or rarely holds residues the background
    commonly holds, is a candidate site of the functional specialization
    that set the clade apart. CladeCompare answers the question one column
    at a time and ranks the columns by the strength of the contrast.

Why a contrast
    The reasoning follows Kannan et al. (2007), who compared the families
    of the protein kinase-like superfamily. Ten key residues are conserved
    across its families and likely mediate the core functions of the
    catalytic domain, while "family-specific functions are mediated by
    features that are highly conserved within families, but that are
    divergent between families". It follows that conservation within a
    clade does not by itself separate the two kinds of residue; a contrast
    with the clade's relatives does, since a residue that the background
    conserves as well belongs to what the clades share, and one that it
    does not marks the clade. The same reasoning applies at each level of
    a hierarchy of families and subfamilies, as in the ROPK comparisons
    that the next entry describes.

Where it comes from
    The measurement is that of CHAIN analysis (Neuwald et al. 2003;
    Neuwald 2007), which measures how the composition of a column in a
    subalignment departs from, and thus contrasts with, the composition of
    the same column in a larger alignment of related sequences. CHAIN
    works on PSI-BLAST alignments; CladeCompare applies the same contrast
    to alignments made with HMMER profiles or MAPGAPS (Neuwald 2009),
    offers several statistics for it, among them CHAIN's ball-in-urn model,
    and draws the results as reports, sequence logos and PyMOL scripts.
    CladeCompare was written for Talevich & Kannan (2013), which compared
    each rhoptry kinase (ROPK) subfamily with the other ROPKs, and the ROPK
    family with the protein kinase superfamily, to find sites of
    subfunctionalization and neofunctionalization, and mapped the
    significant sites onto the structures of ROP2, ROP8 and ROP5.

What it needs
    The contrast is taken column by column, so column N of the foreground
    and column N of the background must be the same position of the
    protein. That holds for two parts of one alignment, and for two sets
    aligned separately to the same profile, with ``hmmalign`` (``--hmm``)
    or MAPGAPS (``--mapgaps`` or its CMA files), where column N is match
    state N of the profile. Two alignments of different widths are first
    aligned to each other with MUSCLE 3. Choosing the clades, and building
    the alignments, are done before CladeCompare runs, for instance with
    MAPGAPS, CHAIN's own partitioning, or a phylogeny.

What it gives
    For each column, the value of the chosen statistic, unadjusted for
    multiple testing; the significant columns, selected by the
    Benjamini-Hochberg procedure at a false discovery rate of alpha and
    capped at 50; paired foreground and background sequence logos around
    them; and, given a structure aligned to the same profile, a PyMOL
    script that highlights them. A column that cannot be tested, because
    the consensus of either set is a gap or an unknown residue, has the
    value 1 and still counts among the columns of the selection. Input
    that CladeCompare recognizes as unusable stops the run with an error
    message, which names the file where one is at fault, rather than a
    traceback; a missing external program, or one that fails on its
    input, still ends with a traceback.

What it does not do
    It does not define clades, infer a phylogeny or align sequences other
    than through the programs named above. It tests each column on its
    own, so it does not detect positions that change together. Its values
    are meant to be read as p-values only where the strategy computes
    p-values, which Methods_ records for each strategy, together with the
    calibration measured where that has been examined.


Installation
------------

CladeCompare requires Python 3.12 or later. It depends on Biopython_,
BioFrills_ (with its ``logo`` extra, which brings in WebLogo_), SciPy_,
ReportLab_ and bottle_, and installs three commands:

``cladecompare``
    compares alignments and writes the results to text files;
``cladereport``
    renders one or more result files of ``cladecompare`` as an HTML page;
``cladeweb``
    serves a form in your web browser that runs a comparison and shows the
    resulting page.

Until BioFrills 0.4.0 is released on PyPI, install CladeCompare from a copy
of this repository. Both routes below need ``git``, because BioFrills is
installed from its repository. With uv_, the dependencies are resolved
automatically, BioFrills included, from the commit recorded in
``pyproject.toml``::

    git clone https://github.com/etal/cladecompare.git
    cd cladecompare
    uv tool install .

With pip, install that BioFrills commit first, then CladeCompare::

    pip install "biofrills[logo] @ git+https://github.com/etal/biofrills.git@1b96a800e7258ee4d6f8b58a4e6d4a3162d5788f"
    pip install .

.. _Biopython: https://biopython.org/
.. _BioFrills: https://github.com/etal/biofrills
.. _WebLogo: https://weblogo.threeplusone.com/
.. _SciPy: https://scipy.org/
.. _ReportLab: https://pypi.org/project/reportlab/
.. _bottle: https://bottlepy.org/
.. _uv: https://docs.astral.sh/uv/

External programs
~~~~~~~~~~~~~~~~~

Three external programs are optional. Each is needed only by the options
that run it, and a program that is missing is reported when it is first run,
not when CladeCompare starts.

HMMER_ 3, any release from 3.0
    ``hmmalign`` aligns unaligned sequences, and the chains of a structure,
    to a profile HMM, for ``--hmm`` and ``--pdb``.
MUSCLE_ 3, such as 3.8.1551
    aligns a foreground and a background alignment of different widths to
    each other. MUSCLE 5 has no profile-profile alignment mode, and
    CladeCompare stops with an error if it finds MUSCLE 4 or later. Bioconda_
    packages MUSCLE 3.8.1551, the version the tests use, and the `MUSCLE 3
    download page`_ offers binaries of 3.8.31.
MAPGAPS_
    ``run_gaps`` aligns sequences to a hierarchy of profiles, for
    ``--mapgaps``. Its distribution holds binaries compiled for Linux on
    x86-64.

The pixi_ environment described under Development below provides HMMER 3.4
and MUSCLE 3.8.1551 from Bioconda, together with CladeCompare itself.

.. _HMMER: https://hmmer.org/
.. _MUSCLE: https://drive5.com/muscle/
.. _MUSCLE 3 download page: https://drive5.com/muscle/downloads_v3.htm
.. _MAPGAPS: https://mapgaps.igs.umaryland.edu/
.. _Bioconda: https://bioconda.github.io/
.. _pixi: https://pixi.prefix.dev/


Usage
-----

The examples below use the two alignments of protein kinase domains in this
repository's ``test`` directory, which come from Talevich, Mirza & Kannan
(2011): ``scttlre-domain``, 14 sequences of the apicomplexan subfamily of
cyclin-dependent kinases whose αC-helix consensus is SCTTLRE and which
includes Pfcrk-5 (PFF0750w), and ``cdc2-domain``, 37 sequences of the CDC2
subfamily, with which that paper contrasted it using CHAIN. Where an
example names an input file that the repository does not hold, such as
``fg.fasta`` or ``subfam1.aln``, the name stands for a file of your own.
``cladecompare --help`` lists every option.

Comparing two alignments
~~~~~~~~~~~~~~~~~~~~~~~~

Give the foreground and the background alignments, in that order::

    cladecompare test/scttlre-domain.fasta_aln.cma test/cdc2-domain.fasta_aln.cma \
        -o scttlre-v-cdc2.out

Without ``-o``, the result is written to standard output. A file whose name
ends in ``.cma``, or any file when ``-f cma`` is given, is read as a CMA
alignment, the format MAPGAPS writes; any other file is read in the format
that ``-f`` names, using Biopython's names for alignment formats (``fasta``
by default, or for instance ``clustal`` or ``stockholm``). The two
alignments should share their columns, as they do when both are aligned to
the same profile. If they differ in width, CladeCompare aligns them to each
other with MUSCLE, and the column numbers then refer to that combined
alignment, which is not written out.

Each line of the result, the ``.out`` file, describes one column::

    V (Y) 1 : prob=1.83534e-10       *********
    R (T) 2 : prob=2.67812e-05       ****
    ...
    K (K) 94 : prob=0.190823

The fields are the consensus residue of the foreground column, that of the
background column in parentheses, the column number, and the value of the
statistic, which for most strategies is a p-value. The value of urn and hypg
refers to the residue that the pattern file names, which can differ from the
consensus, and that of phospho to Ser, Thr and Tyr together, as `The residue
tested`_ explains. A consensus of ``-`` marks a gap column, as `Consensus
residues and untested columns`_ defines it, which is not tested and is given
the value 1. The asterisks, after a tab, show one per order of magnitude of
the value where it is below the significance threshold alpha (``-a``,
between 0 and 1, default 0.005). These values are not adjusted for multiple
testing, so among hundreds of columns some fall below alpha by chance; the
correction applies to the significant sites described below. `Preparing the
alignments`_ describes which columns the numbers refer to.

``-s`` chooses the statistic, as described under Strategies_, and ``-w``
weights the sequences, so that a group of closely related sequences counts
for less than the same number of unrelated ones; Methods_ describes the
weights each strategy uses. With ``-q``, only warnings and errors are
printed.

Pattern files and sequence logos
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``-p``, ``cladecompare`` also writes a pattern file listing the
significant columns, and PDF files of paired sequence logos around them::

    cladecompare test/scttlre-domain.fasta_aln.cma test/cdc2-domain.fasta_aln.cma \
        -o scttlre-v-cdc2.out -p scttlre-v-cdc2.pttrn

The pattern file holds one line in a format modeled on the pattern files of
mcBPPS (Neuwald 2010): ``1:`` and the significant columns, most
significant first, each written as the residue that names the site followed
by the column number::

    1:L135,V1,L158,M123,Y154,L47,G9,I153,...

For gtest and jsd, which compare the whole composition of the columns, the
residue named is the foreground consensus; for urn, hypg and phospho it is
the residue that `The residue tested`_ describes. The significant columns are
those that pass the Benjamini-Hochberg procedure at a false discovery rate of
alpha, taken over all columns of the alignment, and at most 50 columns are
selected, even where the procedure would accept more, as `Significant sites`_
under Methods explains. Each logo file shows the foreground above the
background for a group of significant columns that lie at most 10 columns
apart, with up to 10 columns on either side, and is named after the pattern
file and the columns it highlights, as in
``scttlre-v-cdc2-31,32,34,37,40,43,47,52,53.pdf``.

Aligning sequences to a profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unaligned sequences in FASTA format can be aligned when the comparison
runs. With ``--hmm``, each set of sequences is aligned to a profile HMM with
``hmmalign``, and only the columns of the profile's match states are kept,
so the two sets share their columns, and column N of the output is match
state N of the profile, even where no sequence of either set occupies it::

    cladecompare test/scttlre-domain.fasta test/cdc2-domain.fasta \
        --hmm test/CDK.hmm -o scttlre-v-cdc2.out

With ``--mapgaps``, each set is aligned with MAPGAPS' ``run_gaps``, which
writes its alignment beside the input file, named after it with
``_aln.cma`` added; the profile is named without its file extension, and,
as with ``--hmm``, column N of the output is match state N::

    cladecompare fg.fasta bg.fasta --mapgaps path/to/profile -o fg-v-bg.out

Without these options, align the foreground and background sets together
and then divide the alignment into two files, or align each to the same
profile.

Mapping sites onto a structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``--pdb`` and a profile given with ``--hmm`` or ``--mapgaps``,
``cladecompare`` aligns a structure in PDB format to the same profile and
writes a PyMOL_ script that draws the aligned chain as a cartoon and the
residues of the significant columns as sticks. For instance, with the
structure 1HCK from the Protein Data Bank::

    curl -O https://files.rcsb.org/download/1HCK.pdb
    cladecompare test/scttlre-domain.fasta test/cdc2-domain.fasta \
        --hmm test/CDK.hmm --pdb 1HCK.pdb -o scttlre-v-cdc2.out

The script is written to the structure's file name with ``.pml`` added
(here ``1HCK.pdb.pml``), or to the file that ``--pmlout`` names. The file
name of the structure must give its PDB ID, as in ``1abc.pdb``,
``1abcA.pdb``, ``1abc_A.pdb`` or ``pdb1abc.ent``. Without a profile, the
comparison runs and an error is logged, but no script is written. A
significant column that is a deletion in the structure has no residue to
draw; it is named in a warning and left out of the script. Two
limitations apply in this release: a file holding more than one protein
chain, such as a kinase in complex with a cyclin, stops the run with an
error, so keep only the chain of interest in the file; and ``--pdb`` given
more than once stops the run with an error, after the ``.out`` file and,
with ``-p``, the pattern file and any logos are written.

.. _PyMOL: https://pymol.org/

One alignment
~~~~~~~~~~~~~

Given a single alignment, ``cladecompare`` compares each column with the
amino acid composition of the whole alignment::

    cladecompare test/cdc2-domain.fasta_aln.cma -o cdc2.out

The background consensus in the ``.out`` file is then ``_``. With ``-p``,
the pattern file is written, but no logos are drawn.

More than two alignments
~~~~~~~~~~~~~~~~~~~~~~~~

Given three or more alignments, ``cladecompare`` takes each in turn as the
foreground and all the others together as the background::

    cladecompare subfam1.aln subfam2.aln subfam3.aln

All the alignments must have the same width. For each alignment, the result,
the pattern file and the logos are written to the current directory, named
after the alignment file, as in ``subfam1.aln.out``, ``subfam1.aln.pttrn``
and ``subfam1.aln-12,15.pdf``. ``-o`` and ``-p`` are ignored, with a warning.
With ``--pdb``, every comparison writes its PyMOL script to the same file,
the one named after the structure or by ``--pmlout``, so that the file holds
only the script of the last alignment on the command line.

Viewing the results
~~~~~~~~~~~~~~~~~~~

``cladereport`` renders a ``.out`` file as an HTML page::

    cladereport scttlre-v-cdc2.out -p scttlre-v-cdc2.pttrn -o scttlre-v-cdc2.html

The page shows the foreground and background consensus sequences, with each
foreground residue colored by the unadjusted value of its column: red below
0.1, darker as the value decreases; white from 0.1 to 0.5; blue above 0.5,
for columns that are similar in both sets. Columns whose foreground
consensus is a gap are gray, and those whose background consensus is a gap
are yellow. The colors do not depend on alpha. With ``-p``, an asterisk
above the sequences marks each column of the pattern file. Without ``-o``,
the page is written to standard output.

Given several ``.out`` files, such as those written for `More than two
alignments`_, ``cladereport`` shows the foreground consensus of each as one
row of the page, labeled with the file name without ``.out``::

    cladereport subfam1.aln.out subfam2.aln.out subfam3.aln.out -o family.html

The files must number the same columns. The files of one such run do, and
so do the files of separate runs whose sets were aligned to the same profile
with ``--hmm`` or ``--mapgaps``, since column N of each is match state N;
otherwise ``cladereport`` stops with an error. The background consensus,
which differs between those files, is not shown, although each row colors
yellow the columns where its own background consensus is a gap. With
``-p``, the asterisks of the one pattern file given are drawn above all the
rows. Additional file 2 of Talevich & Kannan (2013) is a page of this kind,
with a row for each of 47 result files.

Web interface
~~~~~~~~~~~~~

``cladeweb`` serves a form at http://localhost:8080/cladecompare, on this
computer only, and opens that address in your web browser once the server
is listening. ``--port`` chooses another port, and ``--no-browser`` leaves
the browser closed. If another ``cladeweb`` already serves that port, the
command stops with an address-in-use error and opens no browser; the earlier
server keeps running. Press Ctrl-C to stop the server.

The form takes one or two sequence files: with one, the comparison is that
of a single alignment, described above. A file whose name ends in ``.cma``
is read as CMA, and any other as FASTA. If an HMM profile, whose file name
must end in ``.hmm``, is given as well, the sequences are aligned to it with
``hmmalign``, as with ``--hmm``, and must then be in FASTA format. The form
also chooses the strategy, alpha and, with its checkbox, the weighting of
``-w``. The submission name becomes the title of the result, which is the
page that ``cladereport`` would write for the ``.out`` and pattern files.
The web interface does not offer the logos, structures or more than two
alignments; use ``cladecompare`` for these.


Strategies
----------

The statistic that compares the columns is chosen with ``-s``. Methods_
below gives the formula each one computes and the sources it follows:

:gtest:
    (default) G-test of homogeneity between the amino acid compositions of
    the FG and BG columns. With a single alignment, G-test of goodness of fit
    of each column's composition against that of the whole alignment.
:urn:
    Ball-in-urn model (binomial), after CHAIN (Neuwald et al. 2003), for
    the count in the FG column of the amino acid type with the smallest tail
    probability, against the frequency of that type in the FG and BG
    together.
:phospho:
    Ball-in-urn model, as in "urn", for the combined count of the
    phosphorylatable residues Ser, Thr and Tyr in the FG column, against
    their combined frequency in the FG and BG columns together.
:jsd:
    Jensen-Shannon divergence between the FG and BG column compositions,
    after INTREPID (Sankararaman & Sjölander 2008). The value reported is a
    score that is smaller for more divergent columns, not a p-value.
:hypg:
    Hypergeometric model: the probability of drawing at least as many
    sequences of the amino acid type with the smallest tail probability as
    the FG column holds, when as many sequences as the FG holds are drawn
    from the FG and BG together. hypg compares two alignments; given one,
    ``cladecompare`` stops with an error.

gtest, urn, phospho, and hypg report p-values. The jsd score is not a
p-value, but in this release it still receives asterisks where it is below
alpha and passes through the same selection of significant sites. For jsd,
the pattern file lists the columns with the smallest scores, without the
meaning of a false discovery rate.


Methods
-------

This section records what each strategy computes, as implemented, and the
published methods it follows or adapts. CladeCompare belongs to a set of
projects by the same author, with the libraries BioFrills_ and biocma_ and the
program fammer_, that implement the central ideas of Andrew Neuwald's CHAIN
and MAPGAPS programs with variations, as `What CladeCompare is for`_
describes.

.. _biocma: https://github.com/etal/biocma
.. _fammer: https://github.com/etal/fammer

Preparing the alignments
~~~~~~~~~~~~~~~~~~~~~~~~

Before comparing a foreground with a background, CladeCompare removes from
the background every sequence whose ID also occurs in the foreground, removes
repeated identical sequences (a background sequence identical to a foreground
sequence is removed from the background), aligns the two sets to each other
with MUSCLE if they differ in width, and, unless the sets were aligned to a
profile with ``--hmm`` or ``--mapgaps``, drops the columns in which every
sequence of both sets has the gap character ``-``. The column numbers of the
``.out`` and pattern files, and of the pair logos, refer to this cleaned
alignment. With a profile, every match state keeps its column, so the PyMOL
script maps each column to the residue that the structure's own alignment to
the profile places in the same match state; a structure aligned to a
different number of match states stops the run with an error before any
output is written.

Counts, weights and pseudocounts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For one alignment column, ``n_a`` is the count of amino acid type ``a`` in the
foreground column and ``m_a`` its count in the background column. The
characters ``-``, ``.`` and ``X`` are not counted. Without ``-w``, each
sequence counts as 1, and ``N_F`` and ``N_B`` are the numbers of foreground and
background sequences. With ``-w``, each sequence counts as its sequence
weight, and ``N_F`` and ``N_B`` are the sums of the weights; urn and phospho
weight the sequences in their own way, which their sections below describe.

The sequence weights are the position-based weights of Henikoff & Henikoff
(1994), computed separately for each set by BioFrills
(``alnutils.sequence_weights`` with the scaling ``none``), which scales them
to sum to an estimate of the effective number of independent sequences in
the set: the number ``k`` of sequences with residues drawn at random from the
20 amino acids whose expected number of distinct types in a column,
``20 (1 - (19/20)^k)``, equals the mean number of distinct residue types
(gaps not counted) over the columns that contribute to the weights, but at
most the number of sequences. Identical sequences thus count as one. The
analysis published in Talevich & Kannan (2013) always weighted the sequences,
with an earlier scaling whose sum fell as the sequences became more diverse;
since October 2013, CladeCompare weights them only when ``-w`` is given.

``f_a`` is the frequency of type ``a`` over the whole alignment (in paired
mode, foreground and background together), weighted in the same way. Several
strategies add pseudocounts of total size ``π = 0.5``, distributed over the
amino acid types in proportion to ``f_a``.

Consensus residues and untested columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The consensus residue of a column in one set is the type ``a`` that maximizes
the relative entropy term ``q_a ln(q_a / b_a)``, where ``q_a`` is the
unweighted frequency of ``a`` in the column and ``b_a`` the weighted frequency
of ``a`` over all columns of the same set; ties go to the alphabetically first
type. A column whose weighted fraction of gap characters, ``-`` or ``.``, is
0.8 or more has a gap as its consensus, except where such columns form a run
at the start or the end of the alignment: those receive a residue consensus
if they hold any aligned residue, however many gaps accompany it. Where an
alignment contains uppercase residues, its lowercase residues are insertions:
a column of lowercase residues and gaps alone has a gap as its consensus,
while a column that holds lowercase residues beside aligned characters, which
are uppercase residues and other characters such as ``*``, is read,
consensus, counts and sequence weights alike, as if it were uppercase. An
alignment without any uppercase residue marks nothing by case and is read as
if it were uppercase.
A column whose foreground or background consensus is a gap, or X, an unknown
residue that the counts skip, is not tested and is reported with p = 1.
The site of a column is named in the pattern file after its foreground
consensus residue for gtest and jsd, and after the residue that the next
section describes for urn, hypg and phospho.

The residue tested
~~~~~~~~~~~~~~~~~~

urn and hypg each test one residue type ``a`` of the foreground column: is
its count ``n_a`` larger than the frequency of ``a`` in the two sets together
makes likely? Each computes the tail probability ``P_a`` that its section
below gives, for each of the ``K`` residue types that the foreground column
holds, and reports ::

    P = min(1, K · min_a P_a)

for the type ``a`` with the smallest ``P_a``, which names the site in the
pattern file; a tie goes to the foreground consensus, then to the
alphabetically first type. The type is chosen from the same column that is
then tested, so the smallest ``P_a`` alone is not a p-value; multiplying by
``K`` is the Bonferroni correction for that choice. A column that holds a
single residue type has ``K = 1``. With a single alignment, urn does the same
over the types of each column.

Until November 2026, urn and hypg tested the foreground consensus without
correction. The consensus is also chosen from the column, so the selection
was the same, and the values it gave were too small to be read as p-values,
the more so the more background sequences there were. With columns drawn
from the compositions of the example pair's own columns, 11 foreground
sequences against 400 in the background, the consensus test of hypg gave
7.1% of values below 0.05, 1.5% below 0.01 and 0.16% below 0.001, and that
of urn 6.7%, 1.3% and 0.13%; with the correction, hypg gives 3.5%, 0.64% and
0.06%, and urn 3.2%, 0.55% and 0.04%. On the example pair itself, 726,000
columns from 3,000 random divisions of the 47 pooled sequences into 11
foreground and 36 background sequences gave, before and after, 4.9% and 1.9%
below 0.05, 0.82% and 0.28% below 0.01, and 0.06% and 0.02% below 0.001 for
hypg; urn, whose pooled frequency includes the foreground count, gave 2.9%
and 0.68%, 0.25% and 0.07%, and 0.007% and 0.002%. The values were smallest
with a single alignment, where the frequency of the type is that of the whole
alignment: for columns of 37 sequences drawn from the amino acid frequencies
of ``cdc2-domain``, urn gave 46% of values below 0.05, 9.5% below 0.01 and
0.9% below 0.001, and gives 2.9%, 0.6% and 0.09%.

Of the columns that the Benjamini-Hochberg procedure accepts on the example
pair at the default alpha, 91 for hypg and 62 for urn, the type tested is the
foreground consensus in all but one, column 180 under hypg: the foreground
holds 6 A and 5 C, its consensus is A, and the background holds 24 A and no
C, so hypg tests C, with ``P = 6.0e-04``, where A gave 0.86.

phospho tests Ser, Thr and Tyr together, a set fixed in advance, so it needs
no correction; it names the site after the most frequent of the three in the
foreground column, the alphabetically first on a tie, or after the
consensus where the column holds none of them.

G-test (gtest)
~~~~~~~~~~~~~~

In paired mode, the foreground and background columns form a 2 × K table of
counts over the ``K`` residue types that occur in either column, and the
G-test of homogeneity (Dunning 1993, section 5.3) tests whether the two sets
draw the column's residues from one distribution. With ``F = Σ n_a`` and
``B = Σ m_a`` the residue totals of the two columns, gaps excluded (unlike
``N_F`` and ``N_B``), and with all sums over the ``K`` types::

    E_a = F (n_a + m_a) / (F + B)        (expected foreground count)
    D_a = B (n_a + m_a) / (F + B)        (expected background count)
    G   = 2 Σ [n_a ln(n_a / E_a) + m_a ln(m_a / D_a)]
    p   = P(χ²_(K-1) ≥ G)

where a term whose count is zero is zero. A column with fewer than two residue
types gives ``p = 1``. No pseudocounts are added, so each p-value depends on
its column alone, and gaps do not enter the table. The upper tail of the
chi-squared distribution is computed by SciPy (``scipy.stats.chi2.sf``). On
the example alignments, 726,000 columns from 3,000 random divisions of the 47
pooled ``scttlre-domain`` and ``cdc2-domain`` sequences into 11 foreground and
36 background sequences, unweighted, gave 4.7% below 0.05, 0.66% below 0.01
and 0.04% below 0.001. Under ``-w`` the counts are sequence weights, whose
sums estimate the effective numbers of sequences, 2.1 and 7.7 for that pair,
so the tables are small and the test seldom rejects: with ``-w``, no column of
the example pair falls below 0.005, the smallest p-value being 0.017, and none
did under the earlier test with the same weights.

Until October 2026, gtest instead tested the foreground counts for goodness
of fit to the background column's composition, including pseudocounts, on 19
degrees of freedom::

    E_a = N_F (m_a + π f_a) / (N_B + π)
    G   = 2 Σ n_a ln(n_a / E_a)        (sum over the types with n_a > 0)

This was the test of Talevich & Kannan (2013), which cites Dunning (1993),
whose test compares two samples. It took the background composition as known,
although it is estimated from the background column, and its expected values
summed to ``N_F`` times the background's share of residues rather than to the
foreground's residues. On the same random divisions it gave 7.2% below 0.05,
3.6% below 0.01 and 1.4% below 0.001. A foreground residue that the background
lacks was compared with its pseudocount alone: in column 95 of the example
pair, where 8 of 11 foreground residues are W and the background has no W,
the earlier test gave 8.6e-21 and the test of homogeneity gives 7.1e-05. In
column 160, all V in the foreground against I, V and L in the background, the
earlier test gave 0.076 and the test of homogeneity gives 1.6e-05.

With a single alignment, each column is compared with the composition of the
whole alignment, by a G-test of goodness of fit whose expected values are the
alignment-wide frequencies scaled to the column's observed total::

    O_a = n_a + π f_a
    E_a = f_a Σ_b O_b
    G   = 2 Σ O_a ln(O_a / E_a)
    p   = P(χ²₁₉ ≥ G)

Observed and expected values have the same total, the residues of the column
plus ``π``, so a column whose composition equals the alignment's gives
``G = 0``. For columns drawn from the alignment-wide composition, these
p-values follow their nominal distribution approximately: of 200,000 such
columns of 37 unweighted sequences, drawn from the frequencies of
``cdc2-domain.fasta_aln.cma``, 6.4% fell below 0.05, 1.0% below 0.01 and
0.06% below 0.001. Until October 2026, the expected values were ``f_a``
themselves, which sum to 1, so ``G`` grew with the number of sequences, and
nearly every column was reported as significant.

Ball-in-urn model (urn)
~~~~~~~~~~~~~~~~~~~~~~~

CHAIN (Neuwald et al. 2003) measures the selective constraint at position
``j`` of a subalignment as the probability of drawing, with replacement, at
least as many conserved residues from the corresponding column of a
superalignment as the subalignment holds: the cumulative binomial
probability ::

    P_j = Σ_{i=c..N} C(N, i) p_j^i (1 - p_j)^(N - i)

where ``c`` and ``N`` are the numbers of conserved residues and of residues in
column ``j`` of the subalignment, and ``p_j`` is the frequency of the
conserved residue in column ``j`` of the superalignment. CHAIN reports
``K_j = 1 / P_j``, the expected number of trials needed to observe this
event; CladeCompare reports ``P_j``.

In paired mode, CladeCompare takes the foreground as the subalignment and the
foreground and background columns together as the superalignment, and for
each residue type ``a`` of the foreground column computes ::

    p_a = (m_a + π f_a + n_a) / (N_B + N + π)
    P_a = Σ_{i=⌈n_a⌉..N} C(N, i) p_a^i (1 - p_a)^(N - i)

with ``N`` the number of foreground sequences, from which `The residue
tested`_ takes the reported value. With a single alignment, ``p_a = f_a``
and ``N`` is the number of sequences.

Under ``-w``, urn weights the sequences as CHAIN does. The foreground,
CHAIN's query family, is not weighted, because CHAIN takes its sequences to
be independent, so ``n_a`` remains a count. The background is weighted with
CHAIN's robust weights: the position-based weights of Henikoff & Henikoff,
normalized to a mean of one-half, with weights greater than one truncated to
one (BioFrills' scaling ``andy``). ``m_a``, ``N_B`` and the background
consensus use these weights, and the frequencies ``f_a`` pool the unweighted
foreground with the weighted background. With a single alignment, which
takes the part of the query family, ``-w`` changes nothing.

CladeCompare departs from CHAIN in these respects:

- **Pseudocounts.** CHAIN adds one pseudocount for each amino acid type;
  CladeCompare adds pseudocounts of total size 0.5, in proportion to ``f_a``.
- **Conserved residues.** CHAIN counts residues biochemically similar to the
  conserved residue as conserved, and takes the conserved residue from its
  query sequence, chosen before the column is examined; CladeCompare counts a
  single residue type, chosen from the foreground column with the correction
  that `The residue tested`_ describes.
- **Superalignment.** CHAIN draws from a separate superalignment;
  CladeCompare pools the foreground and background columns.

Phosphorylatable residues (phospho)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ball-in-urn model is applied to the combined count of Ser, Thr and Tyr
instead of a single residue type, with the pseudocounts of urn. With
``t_F = Σ n_a``, the observed count, and ``t_B = Σ (m_a + π f_a)``, both
summed over ``a`` in {S, T, Y}::

    p = (t_B + t_F) / (N_B + N_F + π)
    P = Σ_{i=⌈t_F⌉..N} C(N, i) p^i (1 - p)^(N - i),   N = N_F

A foreground column without Ser, Thr or Tyr thus has ``P = 1``. With a single
alignment, ``t = Σ n_a`` over {S, T, Y}, ``p = f_S + f_T + f_Y``, and ``N`` is
the number of sequences.

Until October 2026, phospho also added pseudocounts to the foreground count,
computed the sum over ``N + 1`` trials, and divided by ``N_B + N_F + 2``, a
constant left from when the pseudocounts had a total size of 1. A column
without Ser, Thr or Tyr thus received a p-value below 1, and a column of Ser,
Thr and Tyr alone received ``p^(N+1)`` rather than ``p^N``. Without ``-w``,
no p-value of the current form is smaller than the earlier one on the
example alignments.

Under ``-w``, phospho weights the sequences as urn does, following CHAIN: the
foreground is not weighted, so ``N_F`` is the number of foreground sequences,
and the background is weighted with the scaling ``andy``, which ``t_B`` and
``N_B`` carry; ``f_a`` pools the unweighted foreground with the weighted
background. With a single alignment, ``-w`` changes nothing.

Jensen-Shannon divergence (jsd)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In paired mode, ``P`` is the foreground composition ``n_a / Σ n`` and ``Q`` the
background composition including pseudocounts, ``(m_a + π f_a)``, normalized
to sum to 1. With ``M = (P + Q) / 2`` and ``H`` the Shannon entropy in bits::

    JSD   = H(M) - (H(P) + H(Q)) / 2
    score = (1 - JSD)²

With a single alignment, ``Q`` is the uniform distribution over the 20 amino
acid types. The score lies between 0 and 1 and is smaller for more divergent
columns. It is not a p-value, and it does not account for the number of
sequences, so a small foreground set appears divergent partly through
sampling variation. INTREPID (Sankararaman & Sjölander 2008), which this
strategy follows, uses JSD between a column and a background distribution as
a conservation score.

Hypergeometric model (hypg)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In paired mode, for each residue type ``a`` of the foreground column, with
every count and size rounded up to an integer, ``k = ⌈n_a⌉``,
``S = ⌈m_a + n_a⌉`` and ``M = ⌈N_F⌉ + ⌈N_B⌉``, ``P_a`` is the upper tail
``P(X ≥ k)`` of a hypergeometric variable ``X`` counting the residues of
type ``a`` among ``⌈N_F⌉`` sequences drawn without replacement from a
population of ``M`` that holds ``S`` of them: the one-sided Fisher exact
test of the 2 × 2 table of the two sets against ``a`` and the other types.
`The residue tested`_ takes the reported value from these. No
pseudocounts are added. A hypergeometric test needs two sets to draw from,
so hypg has no single-alignment mode. Until November 2026, hypg given one
alignment reported ``⌈c⌉`` divided by the number of sequences, the fraction
of the column that matches the consensus, which is not a p-value; the
Benjamini-Hochberg selection then ranked the least conserved columns first.

Significant sites
~~~~~~~~~~~~~~~~~

The ``.out`` file holds each column's value as computed above, unadjusted.
The pattern file, the pair logos and the residues highlighted in the PyMOL
script come from the Benjamini-Hochberg step-up procedure (Benjamini &
Hochberg 1995) at a false discovery rate of alpha, over all columns of the
cleaned alignment, including the untested columns with p = 1, among them
any match state of a profile that neither set occupies. The step-up search
runs over all columns, and at most 50 of the columns it accepts, the most
significant, are selected, so a selection of fewer than 50 means that the
procedure accepted fewer. On the example pair at the default alpha, without
``-w``, gtest, urn and hypg accept 102, 62 and 91 columns, and each of
their pattern files lists 50. Versions before 0.3.0 confined the search to
the 50 most significant columns, so where the procedure accepted more than
50, they could select fewer than 50, or none: with ``-s hypg -w -a 0.245``,
the example pair accepted 54 columns under the hypg of that time, of which
the earlier search selected 11.

Alignments from MAPGAPS
~~~~~~~~~~~~~~~~~~~~~~~

MAPGAPS (Neuwald 2009) aligns sequences to a hierarchy of profiles and writes
the aligned sequences in CMA format. ``--mapgaps`` runs its ``run_gaps``
program and reads the result; CladeCompare also reads CMA files written
earlier. Lowercase residues in a CMA alignment are insertions relative to the
profile and are not compared.

References
~~~~~~~~~~

- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. *Journal of the Royal
  Statistical Society, Series B* 57(1):289-300.
  doi:10.1111/j.2517-6161.1995.tb02031.x
- Dunning T (1993) Accurate methods for the statistics of surprise and
  coincidence. *Computational Linguistics* 19(1):61-74.
- Henikoff S, Henikoff JG (1994) Position-based sequence weights. *Journal
  of Molecular Biology* 243(4):574-578. doi:10.1016/0022-2836(94)90032-9
- Kannan N, Taylor SS, Zhai Y, Venter JC, Manning G (2007) Structural and
  functional diversity of the microbial kinome. *PLoS Biology* 5(3):e17.
  doi:10.1371/journal.pbio.0050017
- Neuwald AF (2007) The CHAIN program: forging evolutionary links to
  underlying mechanisms. *Trends in Biochemical Sciences* 32(11):487-493.
  doi:10.1016/j.tibs.2007.08.009
- Neuwald AF (2009) Rapid detection, classification and accurate alignment of
  up to a million or more related protein sequences. *Bioinformatics*
  25(15):1869-1875. doi:10.1093/bioinformatics/btp342 (MAPGAPS)
- Neuwald AF (2010) Bayesian classification of residues associated with
  protein functional divergence: Arf and Arf-like GTPases. *Biology Direct*
  5:66. doi:10.1186/1745-6150-5-66 (mcBPPS)
- Neuwald AF, Kannan N, Poleksic A, Hata N, Liu JS (2003) Ran's C-terminal,
  basic patch, and nucleotide exchange mechanisms in light of a canonical
  structure for Rab, Rho, Ras, and Ran GTPases. *Genome Research*
  13(4):673-692. doi:10.1101/gr.862303 (the CHAIN methods)
- Sankararaman S, Sjölander K (2008) INTREPID: INformation-theoretic TREe
  traversal for Protein functional site IDentification. *Bioinformatics*
  24(21):2445-2452. doi:10.1093/bioinformatics/btn474
- Talevich E, Kannan N (2013) Structural and evolutionary adaptation of
  rhoptry kinases and pseudokinases, a family of coccidian virulence factors.
  *BMC Evolutionary Biology* 13:117. doi:10.1186/1471-2148-13-117
- Talevich E, Mirza A, Kannan N (2011) Structural and evolutionary
  divergence of eukaryotic protein kinases in Apicomplexa. *BMC Evolutionary
  Biology* 11:321. doi:10.1186/1471-2148-11-321


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
``tests/data/README.md``. Where a deliberate change moves a whole output
file, as a change to the sequence weights does, the expected output is the
current one, which ``tests/data/replay.py`` writes to ``tests/data/current``;
review such a change against the explanation in ``tests/data/README.md``.

Code is formatted with ``ruff format``, and the formatting-only commits are
listed in ``.git-blame-ignore-revs``. To skip them in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

