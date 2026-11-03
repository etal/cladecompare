# CladeCompare — agent instructions

Two filenames reach this one file: OMP and Codex read it as `AGENTS.md`, while
Claude Code reads it as `CLAUDE.md`, a symlink recorded in git. Where symlinks do
not survive checkout — Git for Windows outside Developer Mode, or a DrvFS mount
such as `/mnt/c` — `CLAUDE.md` arrives as a stub naming this file, so read
`AGENTS.md` directly rather than pasting content into the stub.

## What this is

CladeCompare compares two protein sequence alignments — a *foreground* clade
against a *background* — and reports, for each alignment column, how strongly the
amino-acid composition of the two differs. The intended use is finding the
residues that distinguish one subfamily from its relatives. It accompanies
Talevich & Kannan (2013), *BMC Evolutionary Biology* 13:117,
doi:10.1186/1471-2148-13-117, and is BSD 2-clause.

Three console commands share one library of about 2,200 lines under `cladecomparelib/`, and each is the `main()` function of a module there that keeps the name of the script it replaced:

| command | module | role |
|---|---|---|
| `cladecompare` | `cladecomparelib/cladecompare.py` | command line; the full option surface |
| `cladeweb` | `cladecomparelib/cladeweb.py` | local browser form, served on localhost (`--port`, default 8080) |
| `cladereport` | `cladecomparelib/cladereport.py` | renders an existing `.out` file to HTML, no statistics |

The module names matter: the log format prints `%(module)s`, and the web reference cases record log lines from `cladeweb`.

## Five contrast statistics behind one flag

`-s/--strategy` selects one statistics module per run:

- `gtest` — G-test of column amino-acid frequencies
- `urn` — CHAIN's ball-in-urn statistic
- `jsd` — Jensen–Shannon divergence
- `hypg` — hypergeometric
- `phospho` — Ser/Thr/Tyr phosphorylation sites specifically

`STRATEGIES` in `cladecomparelib/core.py` registers each strategy once, mapping its name to a one-line description and its module, and `DEFAULT_STRATEGY` names the default.
The `-s` choices and help text, the radio buttons of the cladeweb form, the module that `process_args` and the web application run, and the log line that names the strategy all come from it, so adding a sixth means writing the module under `cladecomparelib/` and adding one entry to `STRATEGIES`, together with the module's import in `core.py`.
Before the registry, the strategies were listed in four places, and `hypg` reached the command line but never the web form.
One rule remains outside the registry: under `-w`, a strategy's module may name the `sequence_weights` scalings of the foreground and the background in `WEIGHT_SCALINGS`, with `None` for unweighted, and `urn` and `phospho` do so to follow CHAIN, which weights the background only; every other strategy weights both sets with the `none` scaling (`DEFAULT_WEIGHT_SCALINGS` in `core.py`).
`process_one` gives a single alignment the foreground's scaling, and every size is the sum of the weights.
A module without `compare_one` compares two alignments only, as hypg does: `process_one` raises `InputError` for it, before any output is written.

## `.out` values are raw; Benjamini–Hochberg only selects the pattern

The `.out` file holds each column's unadjusted value of the statistic: a p-value for gtest, urn, phospho and hypg, and a score that is not a p-value for jsd, which nonetheless receives asterisks and passes through the same selection of sites.
`write_pvalues` in `cladecomparelib/core.py` writes the value as the statistic returned it, and its asterisks compare that raw value with `-a/--alpha`, one asterisk per order of magnitude.
`cladereport` and the cladeweb result page color each column from the same raw value (`p2class` in `cladecomparelib/report.py`), at fixed cut points that do not depend on alpha.

The multiple-testing correction exists only in `top_hits` in `core.py`: a Benjamini–Hochberg step-up procedure at false discovery rate alpha, with m equal to the number of columns in `hits`, including the columns that `process_pair` and `process_one` skip and assign p = 1, those whose consensus is a gap or X (`NON_RESIDUES` in `cladecomparelib/shared.py`, the characters that no count includes).
In paired mode that count is taken after `clean_alignments` has dropped the columns that are gaps, `-` or `.`, in every sequence of both sets, except when both sets were aligned to a profile with `--hmm` or `--mapgaps`, where every match state keeps its column; the column numbers in `.out` and `.pttrn` refer to that cleaned alignment.
Its result is the set of significant sites, which becomes the `.pttrn` file when `-p` is given, names the pair logos, and lists the residues highlighted in the PyMOL script.
`make_residues` in `cladecomparelib/pmlscript.py` draws the structure's residue at each site's column, from the first list that `aln_resnums_inserts` returns, which holds one entry per match state of the profile, and logs and leaves out a site whose column is a deletion in the structure.
Site N is entry N - 1 of that list because, under a profile, column N is match state N: `process_args` and the cladeweb application pass `profile_aligned` to `process_pair`, so that `clean_alignments` keeps a match state that no sequence of either set occupies.
Such a column has a gap consensus and p = 1 and changes no other column's value, but it counts in m, so a profile much wider than the sequences can select fewer sites than it would without those columns.
`process_output` raises `InputError` before writing anything when a structure is aligned to a different number of match states than the compared alignment has columns, as with HMMER 3.0 without `--allcol`.
A `.cma` file read directly, without `--mapgaps`, is not marked as profile-aligned, so its all-gap columns are still dropped: the example pair has 7 such columns of 249, and its paired `.out` files have 242 columns.
The step-up search runs over all m ranks, and the 50 most significant of the columns it accepts become the sites, so at most 50 sites are selected and a selection of fewer than 50 means that the procedure accepted fewer; `top_hits` takes the cap as its `N` parameter, which neither the command line nor the web form exposes.
Before November 2026 the search examined only the 50 most significant ranks, so where the procedure accepted more than 50 columns it could select fewer than 50, or none; no expected output of the tests moved when that was fixed, as "Recorded defects" in `tests/data/README.md` explains.

The consequence for a local change: altering how many columns enter `hits`, for instance with a filter that drops columns before scoring, does not rescale any value in `.out`, but it changes m and therefore every threshold `alpha * k / m` (k being a column's rank by p-value).
A gap-handling rule that only changes which columns receive p = 1 leaves m unchanged but changes the ranking.
Either can move which columns reach the pattern, the logos and the PyMOL script.
The raw values are not purely local either: except in paired gtest, pseudocounts come from amino-acid frequencies over the whole alignment (`combined_frequencies` in `cladecomparelib/shared.py` in paired mode, `alnutils.aa_frequencies` in single mode), so a change to which sequences or residues enter those frequencies moves the p-value of every column.
Judge such a change on the whole output, pattern included, never on the column you were looking at.

## The program runs on Python 3; `tests/data/py27/` is the check

The library and the three commands run on Python 3.12 with current Biopython; they no longer run on Python 2.
The port removed `Bio.Alphabet` and `Bio._py3k`, which Biopython dropped in 1.78 and 1.77, and otherwise changed only what Python 3 and current Biopython require, so it moves no number by design.

`tests/data/py27/` holds the outputs of the unmodified code, captured under Python 2.7 and Biopython 1.76 from commit `e265b1f`, for the command line, `cladereport.py` and the web application.
`tests/data/README.md` records the environment and how to regenerate it, the inputs, the defects the reference preserves on purpose, and how to read a mismatch.
A change should reproduce those outputs, or differ from them only by a documented, deliberate change; the statistics are exactly the kind of output where a wrong answer looks entirely plausible.

With HMMER 3.0 or 3.4 and MUSCLE 3.8.1551, under Python 3.12 with the ported biofrills, the ported code reproduces the `.out`, `.pttrn`, `.html` and PyMOL files of every command-line case and the response of every successful web submission byte for byte, both with Biopython 1.83, SciPy 1.11.4 and bottle 0.12.25 and with Biopython 1.88, SciPy 1.18.1 and bottle 0.13.4, apart from the differences below, which come from SciPy, from the later fixes to phospho, hypg and gtest, from the removal of single-alignment hypg, from the scale of the sequence weights in later biofrills, from CHAIN's weighting in urn, from the removal of `_cc_combined.seq`, from the fix to the residues of the PyMOL script, from later changes to two failing web cases, or from the form of logs and tracebacks.
Under both SciPy versions, 13 small p-values of `pair.hypg` differ from the reference in the sixth significant digit; the values from both newer releases agree with exact rational arithmetic, and the reference values, computed with SciPy 1.2.1, do not.
hypg now computes the tail with `hypergeom.sf`, which stays accurate for small tails, rather than as `1 - hypergeom.cdf`; on these inputs the two forms print the same values.
`pair.hypg.w` differs throughout: SciPy 1.11.4 and later return `nan` for the non-integral count of successes that weights produce, and hypg now rounds that count up, as it rounds the other counts; with SciPy 1.18.1 its p-values are finite and its pattern is empty.
Six weighted cases differ from the reference: four since the biofrills pin moved to `4aa8b57`, whose `none` scaling makes the weights of an alignment sum to an estimate of its effective number of independent sequences, and the two `urn` cases because `urn` now weights the sequences as CHAIN does; `tests/data/current/` holds the expected outputs of those cases, and `tests/data/README.md` explains both changes.
The two single-alignment hypg cases, `single.hypg.a.2` and `single.hypg.w.a.2`, now stop with `InputError` and write nothing, because hypg compares two alignments only.
`single.jsd.w.a.2` still reproduces the reference, because single-alignment jsd does not depend on the scale of the weights.
The `-s phospho` cases complete where the reference records a failure, because of the earlier fix to the phospho sample size, and their expected outputs are regenerated, because phospho has since taken urn's pseudocounts: the foreground count has none, the tail runs over N trials rather than N + 1, and a column without Ser, Thr or Tyr gives 1; on the example alignments, no p-value of the unweighted phospho cases is smaller than at the last commit at which phospho ran, which the two legacy cases record.
Single-alignment gtest now scales its expected values, the alignment-wide amino-acid frequencies, to each column's observed total, where it used the frequencies themselves, so that G no longer grows with the number of sequences; `single.gtest.a.001`, `single.gtest.w.a.001` and `hmm.single.gtest` move, and their expected outputs are regenerated.
Paired gtest is now the G-test of homogeneity on the 2 × K table of the two columns' counts, on K - 1 degrees of freedom and without pseudocounts, where it tested the foreground for goodness of fit to the pseudocounted background composition on 19 degrees of freedom, which treated the background composition as known; every paired gtest case that completes moves, and their expected outputs are regenerated.
The four gtest web submissions are checked against the reports of the command-line cases that they equal in the reference.
`unequal.gtest` no longer writes `_cc_combined.seq` into the working directory, nor logs that it did.
The PyMOL script of `pdb.1HCK.gtest` draws the residue at each site's own column, where the reference drew the residue of the following column, and its expected output is regenerated.
`web.pair.hypg` returns status 200 with the report that the command line writes for the same alignments, because the web form now takes its strategies from the same registry as the command line; the reference records status 500.
`web.pdb.1HCK.gtest` still returns status 500, with `InputError` in place of `AssertionError`.
Standard error and logs differ in the line numbers that the log format prints, in the wording of the duplicate-sequence message and of the strategy line for `urn` and `phospho`, and in the temporary file names that the web logs show, and tracebacks take Python 3's own form; `tests/data/README.md` lists these and the few differences that come from the bottle and Biopython versions.
The run with Biopython 1.83 installed biofrills without its dependencies, because biofrills declares a higher floor (see "Packaging and the two development workflows").

## External aligners: any HMMER 3, and MUSCLE 3 only

`hmmalign_a2m` in `cladecomparelib/core.py` is the one place that runs `hmmalign`, for the `--hmm` sequences and for `--pdb` structures alike.
Its output must hold one column per match state of the profile, including match states that no sequence occupies, because the foreground, the background and each structure are aligned to the profile separately and are then compared column by column.
HMMER 3.1 and later always produce that, and reject the `--allcol` option; HMMER 3.0 produces it only with `--allcol`, so the option is passed when `hmmalign -h` lists it.
On the two unaligned FASTA inputs, `scttlre-domain.fasta` and `cdc2-domain.fasta`, HMMER 3.1b2, 3.2.1, 3.3.2 and 3.4 without `--allcol` give a2m output byte-identical to HMMER 3.0 with it, and the reference cases that use hmmalign, including the `--pdb` cases, reproduce under each of them.
HMMER 3.0 without the option drops all-gap match columns: none of 289 for `cdc2-domain.fasta`, one of 289 for `scttlre-domain.fasta`, which then reaches MUSCLE as an alignment of unequal width, and 7 of 289 for chain A of `1HCK.pdb`, which would shift the residues in the PyMOL script if the run continued.
The structure would then have 282 entries, while the alignment that MUSCLE makes of the two sets keeps at least the 289 columns of `cdc2-domain.fasta`, so `process_output` would be expected to stop the run with `InputError`; the column counts were measured with HMMER 3.0, but the whole run has not been.

`combine_alignments` aligns two alignments of unequal width to each other with `muscle -profile`, a mode MUSCLE 5 no longer has, so it reads `muscle -version` first and stops with an error for MUSCLE 4 or later; where the version cannot be read, it warns and runs MUSCLE anyway.
Replacing MUSCLE 3 with another profile-profile aligner would change the output of every unequal-width run.
No missing program is checked for at the start of a run.
`call_quiet` names a program that is not on the PATH when it is first run: for `hmmalign` and `muscle` that is the probe preceding each alignment, and for MAPGAPS it is the `run_gaps` call itself.

## Letter case: uppercase is aligned, lowercase is inserted

CladeCompare follows the convention of HMMER's a2m output and of MAPGAPS and CHAIN's CMA files: an uppercase residue occupies a match state of the profile, a lowercase residue is an insertion relative to the profile, and `-` is a deletion.
In a2m, `.` may pad an insert column, but HMMER 3.4 omits it, so the rows of its a2m output differ in length until the lowercase residues are removed.
biofrills counts `.` as a gap wherever it counts `-`, in the gap threshold of the consensus and in `remove_empty_cols`, as `count_col` in `cladecomparelib/shared.py` and the alignment-wide frequencies already did: on the example pair written with `.` for every gap, every strategy in both modes writes the output of the pair written with `-`, where before biofrills `95c29b6` the 7 columns of `.` alone survived `clean_alignments` and raised m from 242 to 249.
Three input paths in `cladecomparelib` apply the convention to an alignment made with a profile.
For `--hmm`, `hmm_align_and_read` in `cladecomparelib/core.py` drops the lowercase residues of hmmalign's rows, and for a CMA file `cma_blocks` drops those of the aligned part of each sequence line, so in both cases only match columns are compared.
For `--pdb`, the lowercase residues of the structure's own alignment are kept, because they are residues of the chain.
`aln_resnums_inserts` numbers each uppercase residue and records each run of lowercase residues as an insert, `get_aln_offset` locates the whole run of residues, inserts included, in the chain sequence, and `choose_best_aligned` ranks the chains by their uppercase residues, although for more than one aligned chain it currently fails with `KeyError` (`tests/data/README.md`).

An alignment read in any other format goes through `read_aln` and Biopython's `AlignIO`, and the convention applies to it only where its letter case can mark anything; biofrills' consensus then applies it, as described below.
An alignment without any uppercase residue, such as the lowercase output of some aligners, marks no inserts: biofrills' consensus reads such an alignment without regard to case, and `uppercase_if_unmarked` in `cladecomparelib/core.py` uppercases it, so that the consensus residue and the residues counted in its column share one case.
It does so in `read_aln` for each file, so that a lowercase sequence identical to an uppercase one in the other set is purged as a duplicate and the multi-alignment mode can combine a lowercase file with uppercase ones, and again in `process_pair` for each set after `clean_alignments`, which can remove the only uppercase rows of a set.
Without it, the consensus would be uppercase while the counts and the pseudocounts held lowercase types.
Measured on the example pair: with the foreground file lowercased, single-alignment urn, phospho and hypg raised `KeyError` and jsd wrote values above 1; with both files lowercased, paired jsd raised `ZeroDivisionError`; and every paired result, apart from gtest on two lowercase files, differed from that of the uppercase files.
With it, every strategy in both modes, weighted or not, writes the uppercase files' output byte for byte.
`uppercase_if_unmarked` decides with `alnutils.case_marks_inserts`, the test that the consensus itself applies.

An alignment that holds any uppercase residue keeps its letter case, and biofrills' consensus reads it column by column.
A column holding only lowercase residues and gaps is an insert column: it has a gap consensus and is reported with p = 1, but it is still numbered and still counted among the columns of the Benjamini–Hochberg selection.
Its residues enter the alignment-wide frequencies, and with them the pseudocounts, as types of their own.
A column that holds lowercase residues beside aligned characters, which biofrills takes to be any non-gap character that is not lowercase, such as an uppercase residue, `X` or `*`, is read by the consensus as uppercase; `uppercase_mixed_columns` in `cladecomparelib/core.py` finds such columns with the predicates that the consensus applies, `alnutils.case_marks_inserts` and `alnutils.has_aligned_residue`, and uppercases them, in `process_one` and in `process_pair` for each set, with a warning that names the set's mixed columns, before the sequence weights, the counts and the pseudocounts are computed, so that these and the pair logos read the column as the consensus does.
Before `uppercase_mixed_columns`, the statistics counted `c` and `C` as two types, so that a column of only Cys could test as different between the sets, and single-alignment urn and hypg raised `KeyError` where the consensus type occurred nowhere in the alignment in uppercase; such a column now gives the value of the same column in uppercase.
Two sequences of such an input that differ only in their inserts are not purged as duplicates.
The a2m output of hmmalign cannot be read as a FASTA alignment through `AlignIO`, because its rows differ in length, but its default Stockholm output, read with `-f stockholm`, can, and Biopython reads its `.` padding as `-`: `test/cdc2-domain.fasta` aligned to `CDK.hmm` then gives 451 columns, 162 of them insert columns, beside the 289 match states.
Such an alignment is compared with its insert columns in place, which changes the column numbers, the selection threshold and the pseudocounts of every column.

## biofrills is a dependency, pinned to a git commit until 0.4.0

CladeCompare calls seven names from three pure-Python [biofrills](https://github.com/etal/biofrills) modules, in `core.py`, `shared.py` and `pairlogo.py`: `alnutils.sequence_weights`, `alnutils.aa_frequencies`, `alnutils.remove_empty_cols`, `alnutils.case_marks_inserts`, `alnutils.has_aligned_residue`, `consensus.consensus` and `logoutils.aln2logodata`.
The latest release on PyPI, 0.3.1 from 2013, cannot be installed on a current Python, because its `setup.py` always builds a Cython extension that no longer compiles.
The biofrills development branch has since removed that extension and been ported to Python 3: it builds as a pure-Python `py3-none-any` wheel for Python 3.12 and later, and all seven names import from it on Python 3.12, 3.14 and 3.15.
The dependency is therefore `biofrills[logo]`, where the `logo` extra supplies WebLogo, which `pairlogo` reaches through `logoutils`.
Until 0.4.0 is published, `[tool.uv.sources]` takes biofrills from its GitHub repository at a pinned commit, now `1b96a80`, so that `uv.lock` and `pixi.lock` record exactly which biofrills a result came from; moving to a newer commit is a deliberate change of that pin.
A path source to a checkout beside this one was the alternative, but the lock then records only the path, and the checkout may hold unpushed or uncommitted changes.
The pinned commit's version, `0.4.0.dev0`, sorts below `0.4` under PEP 440, so the floor is declared as `>=0.4.0.dev0`: with `>=0.4`, no installer that checks the floor could install the built wheel before 0.4.0 exists.
Raise the floor to `0.4.0` with that release.
uv installs the commit named in `[tool.uv.sources]` without checking it against the floor, and `uv.lock` does not record the floor, so the floor does not check the pin.
pip does not read that table; `README.rst` gives the pip route, which installs the pinned commit first.
Vendoring the biofrills modules was the alternative while the extension blocked installation, and it no longer offers anything.

Moving from the biofrills commit behind the reference (`36684bb`) to the ported one (`ff11246`, version `0.4.0.dev0`) leaves the values that biofrills returns to CladeCompare identical or within a few units in the last place.
This was measured on the committed alignments as the program passes them to biofrills, with its weights and gap threshold, under Python 3.12, 3.14 and 3.15: the `scttlre-domain` and `cdc2-domain` pair in both its `.cma` and `.aln.fasta` forms after `clean_alignments`, and the single alignments `cdc2-domain` (both forms), `cdc2-a`, `cdc2-b` and `cdc2-trimmed`.
Consensus sequences, sequence weights, the output of `remove_empty_cols`, the pair-logo letter counts and the combined amino-acid frequencies of a pair are identical.
Weighted amino-acid frequencies of a single alignment differ by at most 2 units in the last place; the counts code is unchanged apart from `iteritems()` becoming `items()`, so the likely cause is the order in which the normalizing sum is accumulated, together with the compensated `sum()` of Python 3.12, though this was not isolated.
`chisqprob` from `biofrills.stats.chisq`, which gtest called until October 2026, differs in 10 of 4,001 sampled values at 19 degrees of freedom, by at most 2 units in the last place; each difference traces to `math.exp`, which rounds differently on the reference's x86_64 Python 2.7 and on arm64 Python 3.
The entropies and weights that scale the letters of a pair logo, which WebLogo computes (version 3.6 in the reference environment, 3.9 here), differ by at most 8 units in the last place, and the reference does not record the logos' contents.
The `.out` file prints p-values to six significant digits, so differences of this size should change a printed value only where it lies on a rounding boundary.
On the reference cases none did: with `ff11246`, the code of that time reproduced every `.out` file of the reference that does not call SciPy's `hypergeom` exactly.
Later biofrills releases can move numbers here, for instance through changes to consensus gap decisions or to frequency summation, so judge a biofrills upgrade against the reference outputs like any other change.
Moving the pin from `ff11246` to `4aa8b57` moved seven weighted outputs and no unweighted one: the weights kept their ratios and changed only their sum, and CladeCompare of that time, run against biofrills `ff11246` with weights rescaled to the new sums, wrote the new outputs byte for byte; it still does for the two of those cases that remain on the `none` scaling and do not run gtest, while `single.gtest.w.a.001` and `pair.gtest.w` have since moved with the G-test, the two `urn` cases with CHAIN's weighting, and `single.hypg.w.a.2` no longer runs.
Moving the pin from `4aa8b57` to `54c6e70` moved no reference output.
The one change in that range that reaches CladeCompare, biofrills `95c29b6`, applies the HMMER and CHAIN case convention in `consensus` and `remove_empty_cols`: both count `.` as a gap, and `consensus` reads an alignment without uppercase residues regardless of case, which `uppercase_if_unmarked` matches (see "Letter case").
`sequence_weights`, `aa_frequencies` and `aln2logodata` are unchanged and still count `a` and `A` as different types; the other commits in the range change pairwise identity and command-line tools and remove `biofrills.stats`, none of which CladeCompare calls.
`consensus` counts `X` as a residue and can return it, where CladeCompare's counts and frequencies skip it, so `process_one` and `process_pair` give a column whose consensus is `X` p = 1; a biofrills change to that rule should be judged against those two functions.
Moving the pin from `54c6e70` to `1b96a80` moved no reference output.
Four of the seven biofrills commits in that range reach the functions that CladeCompare calls.
Before the biofrills commits `192c370` and `eb8a856`, `remove_empty_cols`, which `clean_alignments` and `combine_alignments` call, raised `ValueError` whenever it dropped a column from a record carrying per-letter annotations, such as the `#=GR` lines of a Stockholm alignment; it now keeps those annotations and the feature locations at the kept columns, neither of which CladeCompare writes into any output.
A pair of equal-width Stockholm files with such lines and an all-gap column therefore passes `clean_alignments` where it failed before.
A later step can still fail on a Stockholm pair with `#=GC` lines, because `clean_alignments` leaves the column annotations of those lines at the original width, and the pair logos slice the alignments together with those annotations.
The biofrills commit `5581325` raises `ValueError` where `aa_counts` asserted the length of its weights, which CladeCompare's calls satisfy because its weights come from `sequence_weights` on the same alignment, and for an alignment without rows in `blocks`, which CladeCompare does not call.
The biofrills commit `1b96a80` makes `case_marks_inserts` and `has_aligned_residue` public without changing them, and `uppercase_if_unmarked` and `uppercase_mixed_columns` call them in place of the copies they held before.
The other three change pairwise identity, `seq-pick` and `alnutils.to_graph`, none of which CladeCompare calls.

gtest takes the upper tail of the chi-squared distribution from `scipy.stats.chi2.sf`.
Until October 2026 it called `chisqprob` from biofrills, adapted from Cephes to avoid a SciPy dependency, whose convergence tolerance `MACHEP = 1e-7` leaves a relative error of up to 3.6e-7 at 1 degree of freedom and 4.1e-8 at 19, enough to change the sixth significant digit of a printed p-value and the ranking that the Benjamini–Hochberg selection depends on.
On the reference cases the change moved one value, with the paired test of that time: column 153 of `cdc2-b.aln.fasta.out` in `multi.gtest` printed 0.528721, the correctly rounded tail, where `chisqprob` gave 0.528722.

## Packaging and the two development workflows

`pyproject.toml` is the only build configuration, with setuptools as the backend.
The floors of biopython, bottle, reportlab and scipy are the versions Ubuntu 26.04 ships; biopython's floor is also the one that every biofrills commit CladeCompare has pinned declares, although CladeCompare's own code needs no Biopython newer than 1.83.
The pair-logo font, DejaVu Sans Mono Bold 2.33, ships as package data in `cladecomparelib/data/` together with its license, which is the reason the license expression names `Bitstream-Vera`.
`pairlogo` loads that copy unconditionally; it consults no system font directory and has no fallback font, so every installation draws the logos with the same font and line height.

uv installs the Python packages from PyPI as `uv.lock` records them, apart from biofrills, which comes from GitHub.
pixi installs them from conda-forge as `pixi.lock` records them, apart from biofrills, WebLogo and the three small packages that WebLogo brings in, which come from GitHub and PyPI, and adds HMMER 3.4 and MUSCLE 3.8.1551 from Bioconda; its `test-integration` task fails rather than skips when either program is missing.
The `py27-oracle` pixi environment, for osx-64 only, holds the versions that captured `tests/data/py27/`; `tests/data/README.md` describes how it regenerates them.

## `tests/` compares the program with the reference; `test/` holds example inputs

`tests/data/` holds the reference outputs described above: every combination of mode, strategy and weighting on the example alignments, and further command-line and web cases.
The example alignments are real alignments from Talevich, Mirza & Kannan (2011), *BMC Evolutionary Biology* 11:321.
`test/` holds them under the names that the examples in `README.rst` use, and `tests/data/make_inputs.py` copies them from there into `tests/data/inputs`.

`uv run pytest` or `pixi run test` runs the suite under `tests/`.
`tests/test_cli_cases.py` replays every recorded command-line case that the current command line can run with the installed `cladecompare` and `cladereport` commands, in a scratch directory holding copies of `tests/data/inputs`; the two legacy phospho cases, whose `-nw` option no longer exists, are kept as a record and not replayed.
`tests/test_web.py` submits every recorded form to the cladeweb WSGI application.
Each command or submission must give the recorded return code or status, each written file and successful response must equal the reference byte for byte, and every p-value written to a `.out` file must be finite and within [0, 1].
Where the code departs from the reference by design, the test module declares the difference — a substituted value, a file no longer written, a return code, a different exception, an expected output regenerated into `tests/data/current/` by `tests/data/replay.py`, or a test of its own — and `tests/data/README.md` explains it.
`tests/test_core.py`, `tests/test_strategies.py`, `tests/test_report.py` and `tests/test_pmlscript.py` test the parsing, site-selection and statistics functions on small constructed inputs, against hand computations and, for the site selection, against `scipy.stats.false_discovery_control`.
A known defect that a test exposes is marked as a strict expected failure, so that fixing it fails the suite until the marker is removed.
The cases that run hmmalign or MUSCLE, and a `--mapgaps` case that replaces MAPGAPS with a stand-in `run_gaps` built on hmmalign, are marked `integration` and are skipped unless the programs are on the PATH; `pixi run test-integration` runs them alone and fails instead of skipping.

## Conventions

- The code predates this project's current Python baseline; it is 2013-era style, ported to Python 3 and formatted with `ruff format`, but not otherwise modernized.
  Modernize deliberately and in separate commits from behavior changes, so a numerical difference is never buried in a reformat.
- `ruff check` and `ruff format --check` are clean; `pixi run lint` runs the first, and `pixi run fmt` applies the formatting.
  `[tool.ruff]` in `pyproject.toml` records why each ignored rule is ignored.
- List a formatting-only commit in `.git-blame-ignore-revs`.
- The trailing spaces inside the PyMOL templates of `pmlscript` are output that `1HCK.pdb.pml` records, so W291 is ignored there rather than fixed.
- Every `zip()` states `strict=`: `True` where the operands are equal in length by construction, `False` only where truncation is the current behavior, with a comment saying why.
- Input that CladeCompare cannot use raises `InputError`, a subclass of `ValueError` defined in `cladecomparelib/__init__.py`.
  `cladecompare.main` catches it, printing the message and exiting with status 1, and two parsers of the significance threshold catch the `InputError` of `core.parse_alpha`: the `-a` parser of `cladecompare` turns it into an argparse usage error, and `cladeweb` returns its own error message.
  `cladeweb` also returns the message of the `InputError` that `process_one` raises for a strategy that compares two alignments only.
  Otherwise `cladeweb` and `cladereport` let it propagate, and a missing program or MUSCLE 4 or later still raises `RuntimeError` with a traceback.
- Every module that logs uses `logger = logging.getLogger(__name__)`, and the three commands configure logging through `configure_logging` in `cladecomparelib/__init__.py`.
- The files under `tests/data/py27/` are generated but committed deliberately, as the reference; regenerate them with `tests/data/capture_py27.py`, for instance through `pixi run -e py27-oracle capture`, rather than editing them.
- The files under `tests/data/current/` are generated from the current code and committed as the expected output of the cases that depart from `tests/data/py27/` by design; regenerate them with `tests/data/replay.py`, review the change against `tests/data/README.md`, and never edit them by hand.
