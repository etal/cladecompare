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

    pip install "biofrills[logo] @ git+https://github.com/etal/biofrills.git@4aa8b57a5192b79707dc11ed8cdaf14b548fc257"
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
.. _MAPGAPS: https://mapgaps.igs.umaryland.edu/

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

``cladeweb`` serves the form at http://localhost:8080/cladecompare, on this
computer only, and opens that address in your web browser once the server is
listening. ``--port`` chooses another port, and ``--no-browser`` leaves the
browser closed. If another ``cladeweb`` already serves that port, the command
stops with an address-in-use error and opens no browser; the earlier server
keeps running.

If you launched the application from the command line, press Ctrl-C (on
Unix-like systems) to stop the web server application.


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

Statistical tests ("-s" options) for column comparison. The section Methods
below gives the formula each one computes and the sources it follows:

:gtest:
    (default) G-test for goodness-of-fit of FG amino acid counts vs. those of
    the BG column. BG frequencies include pseudocounts calculated from the
    amino acid frequencies of the full sequence set.
:urn:
    Ball-in-urn model (binomial), after CHAIN (Neuwald et al. 2003), for
    counts of the "consensus" amino acid type in FG and BG.
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
    sequences of the FG consensus amino acid type as the FG column holds,
    when as many sequences as the FG holds are drawn from the FG and BG
    together. With a single alignment, the fraction of the column that
    matches the consensus, which is not a p-value.


Methods
-------

This section records what each strategy computes, as implemented, and the
published methods it follows or adapts. CladeCompare belongs to a set of
projects by the same author, with the libraries BioFrills_ and biocma_ and the
program fammer_, that implement the central ideas of Andrew Neuwald's CHAIN
(Neuwald et al. 2003; Neuwald 2007) and MAPGAPS (Neuwald 2009) programs, with
variations. CHAIN builds on PSI-BLAST alignments; CladeCompare aligns
sequences to profiles with HMMER, and also accepts alignments from MAPGAPS.

.. _biocma: https://github.com/etal/biocma
.. _fammer: https://github.com/etal/fammer

Preparing the alignments
~~~~~~~~~~~~~~~~~~~~~~~~

Before comparing a foreground with a background, CladeCompare removes from
the background every sequence whose ID also occurs in the foreground, removes
repeated identical sequences (a background sequence identical to a foreground
sequence is removed from the background), aligns the two sets to each other
with MUSCLE if they differ in width, and drops the columns in which every
sequence of both sets has the gap character ``-``. The column numbers of the
``.out`` and pattern files, and of the pair logos, refer to this cleaned
alignment; the PyMOL script maps them to residues through the structure's own
alignment to the profile.

Counts, weights and pseudocounts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For one alignment column, ``n_a`` is the count of amino acid type ``a`` in the
foreground column and ``m_a`` its count in the background column. The
characters ``-``, ``.`` and ``X`` are not counted. Without ``-w``, each
sequence counts as 1, and ``N_F`` and ``N_B`` are the numbers of foreground and
background sequences. With ``-w``, each sequence counts as its sequence
weight, and ``N_F`` and ``N_B`` are the sums of the weights.

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
type. A column whose weighted fraction of ``-`` characters is 0.8 or more has
a gap as its consensus, except where such columns form a run at the start or
the end of the alignment: those receive a residue consensus however many gaps
they hold. A column with no uppercase residue at all has a gap as its
consensus. A column whose foreground or background consensus is a gap is not
tested and is reported with p = 1. Below, ``c*`` is the foreground consensus
residue and ``c = n_c*`` its count in the foreground column.

G-test (gtest)
~~~~~~~~~~~~~~

In paired mode, the foreground counts are compared with the composition of
the background column, including pseudocounts, by a G-test (Dunning 1993)::

    E_a = N_F (m_a + π f_a) / (N_B + π)
    G   = 2 Σ n_a ln(n_a / E_a)        (sum over the types with n_a > 0)
    p   = P(χ²₁₉ ≥ G)

This is the test of Talevich & Kannan (2013). The upper tail of the
chi-squared distribution with 19 degrees of freedom is computed by BioFrills
(``stats.chisq.chisqprob``, adapted from Cephes). The expected values sum to
``N_F`` times the share of the background column, pseudocounts included,
that is not gaps, while the observed counts sum to the residues of the
foreground column, so in columns with gaps the two totals differ. The
background composition is estimated from the background column but enters the
test as a known distribution, so with small sequence sets these p-values tend
to be smaller than a test that accounted for that uncertainty would give.

With a single alignment, each column is compared with the composition of the
whole alignment::

    O_a = n_a + π f_a
    G   = 2 Σ O_a ln(O_a / f_a)
    p   = P(χ²₁₉ ≥ G)

Here the expected values ``f_a`` sum to 1, while the observed values sum to
the number of residues in the column plus ``π``, so ``G`` grows with the
number of sequences, and these values are not calibrated p-values.

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

In paired mode, CladeCompare takes the foreground as the subalignment, the
foreground consensus residue as the conserved residue, and the foreground and
background columns together as the superalignment::

    p_j = (m_c* + π f_c* + c) / (N_B + N + π)
    P   = Σ_{i=k..N} C(N, i) p_j^i (1 - p_j)^(N - i)

with ``k = max(1, ⌈c⌉)`` and ``N`` the number of foreground sequences. With a
single alignment, ``p_j = f_c*``, ``k = ⌈c⌉`` and ``N`` is the number of
sequences. Under ``-w``, ``N`` remains the number of sequences while ``c``
and ``m_c*`` are sums of weights, so the foreground count is understated
relative to the number of trials.

CladeCompare departs from CHAIN in these respects:

- **Weights.** CHAIN does not weight the query family, whose sequences it
  selects from distinct phyla or kingdoms and treats as independent, and it
  weights the superalignment with position-based weights normalized to a mean
  of one-half, with weights greater than one truncated to one. CladeCompare
  either weights neither set or weights both as described above.
- **Pseudocounts.** CHAIN adds one pseudocount for each amino acid type;
  CladeCompare adds pseudocounts of total size 0.5, in proportion to ``f_a``.
- **Conserved residues.** CHAIN counts residues biochemically similar to the
  conserved residue as conserved; CladeCompare counts only the consensus
  residue type.
- **Superalignment.** CHAIN draws from a separate superalignment;
  CladeCompare pools the foreground and background columns.

Phosphorylatable residues (phospho)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ball-in-urn model is applied to the combined count of Ser, Thr and Tyr
instead of the consensus residue. With ``t_F = Σ (n_a + π f_a)`` and
``t_B = Σ (m_a + π f_a)``, both summed over ``a`` in {S, T, Y}::

    p = (t_B + t_F) / (N_B + N_F + 2)
    P = Σ_{i=⌈t_F⌉..N+1} C(N + 1, i) p^i (1 - p)^(N + 1 - i),   N = ⌈N_F⌉

With a single alignment, ``t = Σ (n_a + π f_a)`` over {S, T, Y},
``p = f_S + f_T + f_Y``, ``N`` is the sum of the weights rounded up, and the
sum runs from ``⌈t⌉`` to ``N + 1`` in the same way.

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

In paired mode, with every count and size rounded up to an integer,
``k = max(1, ⌈c⌉)``, ``K = ⌈m_c* + c⌉`` and ``M = ⌈N_F⌉ + ⌈N_B⌉``, the
p-value is the upper tail ``P(X ≥ k)`` of a hypergeometric variable ``X``
counting the consensus residues among ``⌈N_F⌉`` sequences drawn without
replacement from a population of ``M`` that holds ``K`` of them. No
pseudocounts are added. With a single alignment, the value is ``⌈c⌉`` divided
by the number of sequences, the fraction of the column that matches the
consensus, which is not a p-value.

Significant sites
~~~~~~~~~~~~~~~~~

The ``.out`` file holds each column's value as computed above, unadjusted.
The pattern file, the pair logos and the residues highlighted in the PyMOL
script come from the Benjamini-Hochberg step-up procedure (Benjamini &
Hochberg 1995) at a false discovery rate of alpha, over all columns of the
cleaned alignment, including the untested columns with p = 1, and limited
to the 50 most significant columns.

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
- Neuwald AF (2007) The CHAIN program: forging evolutionary links to
  underlying mechanisms. *Trends in Biochemical Sciences* 32(11):487-493.
  doi:10.1016/j.tibs.2007.08.009
- Neuwald AF (2009) Rapid detection, classification and accurate alignment of
  up to a million or more related protein sequences. *Bioinformatics*
  25(15):1869-1875. doi:10.1093/bioinformatics/btp342 (MAPGAPS)
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

