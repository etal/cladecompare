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
| `cladeweb` | `cladecomparelib/cladeweb.py` | local browser form; one server instance per machine |
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
The `-s` choices and help text, the radio buttons of the cladeweb form, the module that `process_args` and the web application run, and the log line that names the strategy all come from it, so adding a sixth means writing the module under `cladecomparelib/` and adding one entry, with that file's module import.
Before the registry, the strategies were listed in four places, and `hypg` reached the command line but never the web form.
One rule remains outside the registry: `process_pair` and `process_one` give `urn` the number of sequences as its sample size, where every other strategy gets the sum of the sequence weights.

## `.out` p-values are raw; Benjamini–Hochberg only selects the pattern

The `.out` file holds each column's unadjusted p-value.
`write_pvalues` in `cladecomparelib/core.py` writes the value as the statistic returned it, and its asterisks compare that raw value with `-a/--alpha`, one asterisk per order of magnitude.
`cladereport` and the cladeweb result page color each column from the same raw value (`p2class` in `cladecomparelib/report.py`), at fixed cut points that do not depend on alpha.

The multiple-testing correction exists only in `top_hits` in `core.py`: a Benjamini–Hochberg step-up procedure at false discovery rate alpha, with m equal to the number of columns in `hits`, including the columns that `process_pair` and `process_one` skip as gaps and assign p = 1.
In paired mode that count is taken after `clean_alignments` has dropped the columns that are gaps in every sequence of both sets; the column numbers in `.out` and `.pttrn` refer to that cleaned alignment.
Its result is the set of significant sites, which becomes the `.pttrn` file when `-p` is given, names the pair logos, and lists the residues highlighted in the PyMOL script.
The step-up search examines only the 50 most significant columns, so at most 50 sites are selected, and where the full procedure would accept more than 50, the result can be smaller than 50.

The consequence for a local change: altering how many columns enter `hits`, for instance with a filter that drops columns before scoring, does not rescale any value in `.out`, but it changes m and therefore every threshold `alpha * k / m` (k being a column's rank by p-value).
A gap-handling rule that only changes which columns receive p = 1 leaves m unchanged but changes the ranking.
Either can move which columns reach the pattern, the logos and the PyMOL script.
The raw values are not purely local either: pseudocounts come from amino-acid frequencies over the whole alignment (`combined_frequencies` in `cladecomparelib/shared.py` in paired mode, `alnutils.aa_frequencies` in single mode), so a change to which sequences or residues enter those frequencies moves the p-value of every column.
Judge such a change on the whole output, pattern included, never on the column you were looking at.

## The program runs on Python 3; `tests/data/py27/` is the check

The library and the three commands run on Python 3.12 with current Biopython; they no longer run on Python 2.
The port removed `Bio.Alphabet` and `Bio._py3k`, which Biopython dropped in 1.78 and 1.77, and otherwise changed only what Python 3 and current Biopython require, so it moves no number by design.

`tests/data/py27/` holds the outputs of the unmodified code, captured under Python 2.7 and Biopython 1.76 from commit `e265b1f`, for the command line, `cladereport.py` and the web application.
`tests/data/README.md` records the environment and how to regenerate it, the inputs, the defects the reference preserves on purpose, and how to read a mismatch.
A change should reproduce those outputs, or differ from them only by a documented, deliberate change; the statistics are exactly the kind of output where a wrong answer looks entirely plausible.

With HMMER 3.0 or 3.4 and MUSCLE 3.8.1551, under Python 3.12 with the ported biofrills, the ported code reproduces the `.out`, `.pttrn`, `.html` and PyMOL files of every command-line case and the response of every successful web submission byte for byte, both with Biopython 1.83, SciPy 1.11.4 and bottle 0.12.25 and with Biopython 1.88, SciPy 1.18.1 and bottle 0.13.4, apart from the differences below, which come from SciPy, from the later fixes to phospho and hypg, from the removal of `_cc_combined.seq`, or from the form of logs and tracebacks.
Under both SciPy versions, 13 small p-values of `pair.hypg` differ from the reference in the sixth significant digit; the values from both newer releases agree with exact rational arithmetic, and the reference values, computed with SciPy 1.2.1, do not.
hypg now computes the tail with `hypergeom.sf`, which stays accurate for small tails, rather than as `1 - hypergeom.cdf`; on these inputs the two forms print the same values.
`pair.hypg.w` differs throughout: SciPy 1.11.4 and later return `nan` for the non-integral count of successes that weights produce, and hypg now rounds that count up, as it rounds the other counts; with SciPy 1.18.1 its p-values are finite and its pattern is empty.
The `-s phospho` cases now complete where the reference records a failure, because of the earlier fix to the phospho sample size, and their unweighted outputs equal the legacy reference cases.
`unequal.gtest` no longer writes `_cc_combined.seq` into the working directory, nor logs that it did.
Standard error and logs differ in the line numbers that the log format prints, and tracebacks take Python 3's own form; `tests/data/README.md` lists these and the few differences that come from the bottle and Biopython versions.
biofrills `ff11246` declares `biopython>=1.86`, so the run with Biopython 1.83 installed biofrills without its dependencies; CladeCompare's own code needs no Biopython newer than 1.83.

## External aligners: any HMMER 3, and MUSCLE 3 only

`hmmalign_a2m` in `cladecomparelib/core.py` is the one place that runs `hmmalign`, for the `--hmm` sequences and for `--pdb` structures alike.
Its output must hold one column per match state of the profile, including match states that no sequence occupies, because the foreground, the background and each structure are aligned to the profile separately and are then compared column by column.
HMMER 3.1 and later always produce that, and reject the `--allcol` option; HMMER 3.0 produces it only with `--allcol`, so the option is passed when `hmmalign -h` lists it.
On the two unaligned FASTA inputs, `scttlre-domain.fasta` and `cdc2-domain.fasta`, HMMER 3.1b2, 3.2.1, 3.3.2 and 3.4 without `--allcol` give a2m output byte-identical to HMMER 3.0 with it, and the reference cases that use hmmalign, including the `--pdb` cases, reproduce under each of them.
HMMER 3.0 without the option drops all-gap match columns: one of 289 for `scttlre-domain.fasta`, which then reaches MUSCLE as an alignment of unequal width, and 7 of 289 for chain A of `1HCK.pdb`, which shifts the residues in the PyMOL script without any error.

`combine_alignments` aligns two alignments of unequal width to each other with `muscle -profile`, a mode MUSCLE 5 no longer has, so it reads `muscle -version` first and stops with an error for MUSCLE 4 or later; where the version cannot be read, it warns and runs MUSCLE anyway.
Replacing MUSCLE 3 with another profile-profile aligner would change the output of every unequal-width run.
No missing program is checked for at the start of a run.
`call_quiet` names a program that is not on the PATH when it is first run: for `hmmalign` and `muscle` that is the probe preceding each alignment, and for MAPGAPS it is the `run_gaps` call itself.

## biofrills is a dependency, pinned to a git commit until 0.4.0

CladeCompare calls six names from four pure-Python [biofrills](https://github.com/etal/biofrills) modules, in `core.py`, `shared.py`, `gtest.py` and `pairlogo.py`: `alnutils.sequence_weights`, `alnutils.aa_frequencies`, `alnutils.remove_empty_cols`, `consensus.consensus`, `logoutils.aln2logodata` and `stats.chisq.chisqprob`.
The latest release on PyPI, 0.3.1 from 2013, cannot be installed on a current Python, because its `setup.py` always builds a Cython extension that no longer compiles.
The biofrills development branch has since removed that extension and been ported to Python 3: it builds as a pure-Python `py3-none-any` wheel for Python 3.12 and later, and all six names import from it on Python 3.12, 3.14 and 3.15.
The dependency is therefore `biofrills[logo]`, where the `logo` extra supplies WebLogo, which `pairlogo` reaches through `logoutils`.
Until 0.4.0 is published, `[tool.uv.sources]` takes biofrills from its GitHub repository at a pinned commit, starting from `ff11246`, the commit measured below, so that `uv.lock` and `pixi.lock` record exactly which biofrills a result came from; moving to a newer commit is a deliberate change of that pin.
A path source to a checkout beside this one was the alternative, but the lock then records only the path, and the checkout may hold unpushed or uncommitted changes.
The pinned commit's version, `0.4.0.dev0`, sorts below `0.4` under PEP 440, so the floor is declared as `>=0.4.0.dev0`: with `>=0.4`, no installer that checks the floor could install the built wheel before 0.4.0 exists.
Raise the floor to `0.4.0` with that release.
uv installs the commit named in `[tool.uv.sources]` without checking it against the floor, and `uv.lock` does not record the floor, so the floor does not check the pin.
pip does not read that table; `README.rst` gives the pip route, which installs the pinned commit first.
Vendoring the four modules was the alternative while the extension blocked installation, and it no longer offers anything.

Moving from the biofrills commit behind the reference (`36684bb`) to the ported one (`ff11246`, version `0.4.0.dev0`) leaves the values that biofrills returns to CladeCompare identical or within a few units in the last place.
This was measured on the committed alignments as the program passes them to biofrills, with its weights and gap threshold, under Python 3.12, 3.14 and 3.15: the `scttlre-domain` and `cdc2-domain` pair in both its `.cma` and `.aln.fasta` forms after `clean_alignments`, and the single alignments `cdc2-domain` (both forms), `cdc2-a`, `cdc2-b` and `cdc2-trimmed`.
Consensus sequences, sequence weights, the output of `remove_empty_cols`, the pair-logo letter counts and the combined amino-acid frequencies of a pair are identical.
Weighted amino-acid frequencies of a single alignment differ by at most 2 units in the last place; the counts code is unchanged apart from `iteritems()` becoming `items()`, so the likely cause is the order in which the normalizing sum is accumulated, together with the compensated `sum()` of Python 3.12, though this was not isolated.
`chisqprob` differs in 10 of 4,001 sampled values at 19 degrees of freedom, by at most 2 units in the last place; each difference traces to `math.exp`, which rounds differently on the reference's x86_64 Python 2.7 and on arm64 Python 3.
The entropies and weights that scale the letters of a pair logo, which WebLogo computes (version 3.6 in the reference environment, 3.9 here), differ by at most 8 units in the last place, and the reference does not record the logos' contents.
The `.out` file prints p-values to six significant digits, so differences of this size should change a printed value only where it lies on a rounding boundary.
On the reference cases none does: with `ff11246`, every `.out` file of the reference that does not call SciPy's `hypergeom` is reproduced exactly.
Later biofrills releases can move numbers here, for instance through changes to consensus gap decisions or to frequency summation, so judge a biofrills upgrade against the reference outputs like any other change.

biofrills implements `chisqprob` itself, adapted from Cephes to avoid a SciPy dependency, and SciPy removed its own `chisqprob` in 1.0, so `scipy.stats.chi2.sf` is the nearest SciPy function.
Do not substitute one for the other in passing: at the 19 degrees of freedom the G-test uses, they differ by up to 4.0e-8 relative, which reflects the convergence tolerance `MACHEP = 1e-7` in biofrills, and which can change the sixth significant digit of a printed p-value and the ranking that the Benjamini–Hochberg selection depends on.

## Packaging and the two development workflows

`pyproject.toml` is the only build configuration, with setuptools as the backend.
The floors of biopython, bottle, reportlab and scipy are the versions Ubuntu 26.04 ships; biopython's floor is also the one that biofrills `ff11246` declares, although CladeCompare's own code needs no Biopython newer than 1.83.
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
`tests/test_cli_cases.py` replays every recorded command-line case that the current command line can run with the installed `cladecompare` and `cladereport` commands, in a scratch directory holding copies of `tests/data/inputs`; the two legacy phospho cases, whose `-nw` option no longer exists, serve as the expected output of the unweighted phospho cases.
`tests/test_web.py` submits every recorded form to the cladeweb WSGI application.
Each command or submission must give the recorded return code or status, each written file and successful response must equal the reference byte for byte, and every p-value written to a `.out` file must be finite and within [0, 1].
Where the code departs from the reference by design, the test module declares the difference — a substituted value, a file no longer written, a different exception, or a test of its own — and `tests/data/README.md` explains it.
`tests/test_core.py`, `tests/test_strategies.py`, `tests/test_report.py` and `tests/test_pmlscript.py` test the parsing, site-selection and statistics functions on small constructed inputs, against hand computations and, for the site selection, against `scipy.stats.false_discovery_control`.
A known defect that a test exposes is marked as a strict expected failure, so that fixing it fails the suite until the marker is removed.
The cases that run hmmalign or MUSCLE, and a `--mapgaps` case that replaces MAPGAPS with a stand-in `run_gaps` built on hmmalign, are marked `integration` and are skipped unless the programs are on the PATH; `pixi run test-integration` runs them alone and fails instead of skipping.

## Conventions

- The code predates this project's current Python baseline; it is 2013-era style, ported to Python 3 and formatted with `ruff format`, but not otherwise modernized.
  Modernize deliberately and in separate commits from behavior changes, so a numerical difference is never buried in a reformat.
- `ruff check` and `ruff format --check` are clean (`pixi run lint`, `pixi run fmt`), and `[tool.ruff]` in `pyproject.toml` records why each ignored rule is ignored.
  List a formatting-only commit in `.git-blame-ignore-revs`.
  The trailing spaces inside the PyMOL templates of `pmlscript` are output that `1HCK.pdb.pml` records, so W291 is ignored there rather than fixed.
  Every `zip()` states `strict=`: `True` where the operands are equal in length by construction, `False` only where truncation is the current behavior, with a comment saying why.
- The files under `tests/data/py27/` are generated but committed deliberately, as the reference; regenerate them with `tests/data/capture_py27.py`, for instance through `pixi run -e py27-oracle capture`, rather than editing them.
- Maintainer planning for this repository lives in a local beads store that is excluded from git and is not present in a clone.
  If you have `bd`, run `bd prime --mcp`, which renders each memory as a one-line summary because the full text printed by plain `bd prime` can exceed an agent's tool-output limits; `bd memories <keyword>` searches the full text, and `bd recall <key>` prints one memory in full.
  If you do not have `bd`, nothing here depends on it.
