# CladeCompare — agent instructions

Two filenames reach this one file: OMP and Codex read it as `AGENTS.md`, while
Claude Code reads it as `CLAUDE.md`, a symlink recorded in git. Where symlinks do
not survive checkout — Git for Windows outside Developer Mode, or a DrvFS mount
such as `/mnt/c` — `CLAUDE.md` arrives as a stub naming this file, so read
`AGENTS.md` directly rather than pasting content into the stub.

## What this is

CladeCompare compares two protein sequence alignments — a *foreground* clade against a *background* — and reports, for each alignment column, how strongly the amino-acid composition of the two differs, in order to find the residues that distinguish one subfamily from its relatives.
`README.rst`, "What CladeCompare is for", states the purpose, its lineage in CHAIN and the reference to Talevich & Kannan (2013); the license is BSD 2-clause.

Three console commands share one library of about 2,750 lines under `cladecompare/`, and each is the `main()` function of a module there that keeps the name of the script it replaced:

| command | module | role |
|---|---|---|
| `cladecompare` | `cladecompare/cladecompare.py` | command line; the full option surface |
| `cladeweb` | `cladecompare/cladeweb.py` | local browser form, served on localhost (`--port`, default 8080) |
| `cladereport` | `cladecompare/cladereport.py` | renders one or more existing `.out` files to HTML, no statistics |

The module names matter: the log format prints `%(module)s`, and the web reference cases record log lines from `cladeweb`.
The import package carries the distribution's name; from 2013 until 0.3.0 it was `cladecomparelib`, because the scripts that `setup.py` installed into `bin/`, such as `cladecompare.py`, would have been imported in place of a package named `cladecompare`, and the console commands that replaced them cannot be imported.
The tracebacks recorded in `tests/data/py27/` therefore name `cladecomparelib`.

## Five contrast statistics behind one flag

`-s/--strategy` selects one statistics module per run:

- `gtest` — G-test of column amino-acid frequencies
- `urn` — CHAIN's ball-in-urn statistic
- `jsd` — Jensen–Shannon divergence
- `hypg` — hypergeometric
- `phospho` — Ser/Thr/Tyr phosphorylation sites specifically

`STRATEGIES` in `cladecompare/core.py` registers each strategy once, mapping its name to a one-line description and its module, and `DEFAULT_STRATEGY` names the default.
The `-s` choices and help text, the radio buttons of the cladeweb form, the module that `process_args` and the web application run, and the log line that names the strategy all come from it, so adding a sixth means writing the module under `cladecompare/` and adding one entry to `STRATEGIES`, together with the module's import in `core.py`.
One rule remains outside the registry: under `-w`, a strategy's module may name the `sequence_weights` scalings of the foreground and the background in `WEIGHT_SCALINGS`, with `None` for unweighted, and `urn` and `phospho` do so to follow CHAIN, which weights the background only; every other strategy weights both sets with the `none` scaling (`DEFAULT_WEIGHT_SCALINGS` in `core.py`).
`process_one` gives a single alignment the foreground's scaling, and every size is the sum of the weights.
A module without `compare_one` compares two alignments only, as hypg does: `process_one` raises `InputError` for it, before any output is written.
`compare_cols` and `compare_one` return the residue that names the site together with the value: gtest and jsd, which compare whole compositions, return the foreground consensus they are given; urn and hypg return the residue type of the foreground column with the smallest tail, and as the value that tail multiplied by the number of types in the column, at most 1 (`most_significant_type` in `cladecompare/shared.py`), because the type is chosen from the column it is tested on; phospho returns the most frequent of Ser, Thr and Tyr.
Each entry of `hits` is therefore (foreground consensus, background consensus, value, residue): the `.out` file prints the consensus, which `cladereport` draws as the consensus sequence, and the `.pttrn` file names each site after the residue.

## `.out` values are raw; Benjamini–Hochberg only selects the pattern

The `.out` file holds each column's unadjusted value of the statistic: a p-value for gtest, urn, phospho and hypg, and a score that is not a p-value for jsd, which nonetheless receives asterisks and passes through the same selection of sites.
`write_pvalues` in `cladecompare/core.py` writes the value as the statistic returned it, and its asterisks compare that raw value with `-a/--alpha`, one asterisk per order of magnitude.
`cladereport` and the cladeweb result page color each column from the same raw value (`p2class` in `cladecompare/report.py`), at fixed cut points that do not depend on alpha.
`report.parse_out` reads a `.out` file for both, and raises `InputError` for a file without lines, a line that `write_pvalues` would not write, or a value that is not finite; a finite value outside [0, 1], of which `write_pvalues` warns, is colored as the nearest bound.
A statistic that returns `nan` writes it to the `.out` file on purpose (`most_significant_type` in `cladecompare/shared.py`), so that a numerical failure is not hidden; `cladereport` then names the file and line instead of drawing a report with that column colored as if it had a value.
A pattern file marks the columns that carry its sites' numbers, and a site that is not among the `.out` file's columns raises `InputError`.

The multiple-testing correction exists only in `top_hits` in `core.py`: a Benjamini–Hochberg step-up procedure at false discovery rate alpha, with m equal to the number of columns in `hits`, including the columns that `process_pair` and `process_one` skip and assign p = 1, those whose consensus is a gap or X (`NON_RESIDUES` in `cladecompare/shared.py`, the characters that no count includes).
In paired mode that count is taken after `clean_alignments` has removed the insert columns of lowercase residues and gaps ("Letter case" below) and dropped the columns that are gaps, `-` or `.`, in every sequence of both sets, except when both sets were aligned to a profile with `--hmm` or `--mapgaps`, where every match state keeps its column; the column numbers in `.out` and `.pttrn` refer to that cleaned alignment.
Its result is the set of significant sites, which becomes the `.pttrn` file when `-p` is given, names the pair logos, and lists the residues highlighted in the PyMOL script.
`make_residues` in `cladecompare/pmlscript.py` draws the structure's residue at each site's column, from the first list that `aln_resnums_inserts` returns, which holds one entry per match state of the profile, and logs and leaves out a site whose column is a deletion in the structure.
Site N is entry N - 1 of that list because, under a profile, column N is match state N: `process_args` and the cladeweb application pass `profile_aligned` to `process_pair`, so that `clean_alignments` keeps a match state that no sequence of either set occupies.
Such a column has a gap consensus and p = 1 and changes no other column's value, but it counts in m, so a profile much wider than the sequences can select fewer sites than it would without those columns.
`process_output` raises `InputError` before writing anything when a structure is aligned to a different number of match states than the compared alignment has columns, as with HMMER 3.0 without `--allcol`.
A `.cma` file read directly, without `--mapgaps`, is not marked as profile-aligned, so its all-gap columns are still dropped: the example pair has 7 such columns of 249, and its paired `.out` files have 242 columns.
The step-up search runs over all m ranks, and the 50 most significant of the columns it accepts become the sites, so at most 50 sites are selected and a selection of fewer than 50 means that the procedure accepted fewer; `top_hits` takes the cap as its `N` parameter, which neither the command line nor the web form exposes.

The consequence for a local change: altering how many columns enter `hits`, for instance with a filter that drops columns before scoring, does not rescale any value in `.out`, but it changes m and therefore every threshold `alpha * k / m` (k being a column's rank by p-value).
A gap-handling rule that only changes which columns receive p = 1 leaves m unchanged but changes the ranking.
Either can move which columns reach the pattern, the logos and the PyMOL script.
The raw values are not purely local either: except in paired gtest, pseudocounts come from amino-acid frequencies over the whole alignment (`combined_frequencies` in `cladecompare/shared.py` in paired mode, `alnutils.aa_frequencies` in single mode), so a change to which sequences or residues enter those frequencies moves the p-value of every column.
Judge such a change on the whole output, pattern included, never on the column you were looking at.

## The program runs on Python 3; `tests/data/py27/` is the check

The library and the three commands run on Python 3.12 with current Biopython; they no longer run on Python 2.
The port removed `Bio.Alphabet` and `Bio._py3k`, which Biopython dropped in 1.78 and 1.77, and otherwise changed only what Python 3 and current Biopython require, so it moves no number by design.

`tests/data/py27/` holds the outputs of the unmodified code, captured under Python 2.7 and Biopython 1.76 from commit `e265b1f`, for the command line, `cladereport.py` and the web application.
A change should reproduce those outputs, or differ from them only by a documented, deliberate change; the statistics are exactly the kind of output where a wrong answer looks entirely plausible.
`tests/data/README.md` records the environment and how to regenerate it, the inputs, the cases that differ from the reference and why ("Regenerated outputs in `current/`", "Recorded defects"), the environments under which the rest reproduces byte for byte, and how to read a mismatch ("Interpreting a mismatch").
A new deliberate difference is recorded there and declared in the test module, never left as an unexplained mismatch.

## External aligners: any HMMER 3, and MUSCLE 3 only

`hmmalign_a2m` in `cladecompare/core.py` is the one place that runs `hmmalign`, for the `--hmm` sequences and for `--pdb` structures alike.
Its output must hold one column per match state of the profile, including match states that no sequence occupies, because the foreground, the background and each structure are aligned to the profile separately and are then compared column by column.
HMMER 3.1 and later always produce that, and reject the `--allcol` option; HMMER 3.0 produces it only with `--allcol`, so the option is passed when `hmmalign -h` lists it.
HMMER 3.0 without the option drops all-gap match columns, which leaves alignments of unequal width or a structure with too few entries; `tests/data/README.md`, "HMMER version", gives the measurements.

`combine_alignments` aligns two alignments of unequal width to each other with `muscle -profile`, a mode MUSCLE 5 no longer has, so it reads `muscle -version` first and stops with an error for MUSCLE 4 or later; where the version cannot be read, it warns and runs MUSCLE anyway.
Replacing MUSCLE 3 with another profile-profile aligner would change the output of every unequal-width run.
No missing program is checked for at the start of a run.
`call_quiet` names a program that is not on the PATH when it is first run: for `hmmalign` and `muscle` that is the probe preceding each alignment, and for MAPGAPS it is the `run_gaps` call itself.

## Letter case: uppercase is aligned, lowercase is inserted

CladeCompare follows the convention of HMMER's a2m output and of MAPGAPS and CHAIN's CMA files: an uppercase residue occupies a match state of the profile, a lowercase residue is an insertion relative to the profile, and `-` is a deletion.
In a2m, `.` may pad an insert column, but HMMER 3.4 omits it, so the rows of its a2m output differ in length until the lowercase residues are removed.
biofrills counts `.` as a gap wherever it counts `-`, in the gap threshold of the consensus and in `remove_empty_cols`, as `count_col` in `cladecompare/shared.py` and the alignment-wide frequencies already did, so an alignment written with `.` for every gap gives the output of the same alignment written with `-`.
Three input paths in the library apply the convention to an alignment made with a profile.
For `--hmm`, `hmm_align_and_read` in `cladecompare/core.py` drops the lowercase residues of hmmalign's rows, and for a CMA file `cma_blocks` drops those of the aligned part of each sequence line, so in both cases only match columns are compared.
For `--pdb`, the lowercase residues of the structure's own alignment are kept, because they are residues of the chain.
`aln_resnums_inserts` numbers each uppercase residue and records each run of lowercase residues as an insert, `get_aln_offset` locates the whole run of residues, inserts included, in the chain sequence, and `choose_chain` takes the chain that `--chain` names or, by default, of the structure's chains the one with the most uppercase residues other than X, and of chains with equally many the one whose chain ID sorts first; only that chain is located and drawn.

An alignment read in any other format goes through `read_aln` and Biopython's `AlignIO`, and `match_columns` in `cladecompare/core.py` applies the convention to it before anything counts, weights or numbers a column: `process_one` calls it for the alignment, `clean_alignments` for the two sets before it removes duplicates, and the multi-alignment mode of `process_args` for all the files before any background joins two of them.
`read_aln` drops the column annotations, such as the `#=GC` lines of a Stockholm file, which nothing reads: once `clean_alignments` has removed a column, an annotation that keeps the input width makes Biopython refuse a slice that reaches past the last column, as a pair logo window can.
An alignment without any uppercase residue, such as the lowercase output of some aligners, marks no inserts: biofrills' consensus reads such an alignment without regard to case, and `uppercase_if_unmarked` uppercases it, for each file or set separately, so that the case of one file is never read through that of another.

In an alignment that holds any uppercase residue, a column of lowercase residues and gaps, holding no aligned character, is an insert column, and `match_columns` removes it, as `hmm_align_and_read` and `cma_blocks` remove lowercase residues.
It is therefore not numbered and not counted among the columns of the Benjamini–Hochberg selection, its residues do not enter the alignment-wide frequencies or the pseudocounts, and two sequences that differ only in their inserts are purged as duplicates.
A column that holds lowercase residues beside aligned characters, which biofrills takes to be any non-gap character that is not lowercase, such as an uppercase residue, `X` or `*`, is read by the consensus as uppercase, and `match_columns` uppercases it, with a warning that names the set's mixed columns, so that the sequence weights, the counts, the pseudocounts and the pair logos read the column as the consensus does.
Such a column gives the value of the same column in uppercase, where the statistics would otherwise count `c` and `C` as two types.
Both decisions use the predicates that the consensus applies, `alnutils.case_marks_inserts` and `alnutils.has_aligned_residue`, and no lowercase residue remains after `match_columns`.

Alignments of equal width are taken to share their columns, as `clean_alignments` already assumes when it compares them without MUSCLE.
A column that is an insert column of one of them is removed from all, with a warning that the alignments are read as parts of one, so that two parts of one alignment, where an insert column of one part holds only gaps in the other, keep the same match columns.
Where any such column holds an aligned character in another alignment, the alignments do not share their columns, and each loses only its own insert columns, as alignments of different widths do.
The rule cannot tell two parts of one alignment from two alignments made separately whose widths are equal by chance and whose insert columns all fall on columns that the other leaves empty; it reads them as parts of one.
Two alignments made separately with one profile, such as the default Stockholm output of hmmalign for each set, read with `-f stockholm`, thereby both arrive at the profile's match states and are compared without MUSCLE.
In the multi-alignment mode, alignments that still differ in width raise `InputError`, because each background joins the rows of the other alignments, and `combine_alignments` aligns two sets only.
On the example pair, that Stockholm output gives the output of `--hmm` byte for byte (`tests/data/README.md`, "Letter case and `.` gaps").
The exception is a match state that no sequence of either set occupies: without the profile it is an all-gap column, which `clean_alignments` drops where `--hmm` keeps it, so that the later column numbers and m differ from those of `--hmm`.
An alignment whose letter case no longer marks its inserts, as in a subset that `esl-alimanip --seq-k` writes, is compared with its insert columns in place, because its `#=GC RF` line, which marks them, is not read.

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

Judge a biofrills upgrade against the reference outputs like any other change, because it can move numbers, for instance through changes to consensus gap decisions or to frequency summation; `tests/data/README.md`, "biofrills version", records what each pin moved and what it left unchanged.
Two behaviors of biofrills are standing assumptions of CladeCompare's code:

- `sequence_weights`, `aa_frequencies` and `aln2logodata` count `a` and `A` as different types, which `match_columns` accounts for.
- `consensus` counts `X` as a residue and can return it, where CladeCompare's counts and frequencies skip it, so `process_one` and `process_pair` give a column whose consensus is `X` p = 1; a biofrills change to that rule should be judged against those two functions.

gtest takes the upper tail of the chi-squared distribution from `scipy.stats.chi2.sf`, not from biofrills' `chisqprob`, whose convergence tolerance is large enough to change the sixth significant digit of a printed p-value and with it the ranking that the Benjamini–Hochberg selection depends on.

## Packaging and the two development workflows

`pyproject.toml` is the only build configuration, with setuptools as the backend.
The floors of biopython, bottle, reportlab and scipy are the versions Ubuntu 26.04 ships; biopython's floor is also the one that every biofrills commit CladeCompare has pinned declares, although CladeCompare's own code needs no Biopython newer than 1.83.
The pair-logo font, DejaVu Sans Mono Bold 2.33, ships as package data in `cladecompare/data/` together with its license, which is the reason the license expression names `Bitstream-Vera`.
`pairlogo` loads that copy unconditionally; it consults no system font directory and has no fallback font, so every installation draws the logos with the same font and line height.

uv installs the Python packages from PyPI as `uv.lock` records them, apart from biofrills, which comes from GitHub.
pixi installs them from conda-forge as `pixi.lock` records them, apart from biofrills, WebLogo and the three small packages that WebLogo brings in, which come from GitHub and PyPI, and adds HMMER 3.4 and MUSCLE 3.8.1551 from Bioconda; its `test-integration` task fails rather than skips when either program is missing.
The `py27-oracle` pixi environment, for osx-64 only, holds the versions that captured `tests/data/py27/`; `tests/data/README.md` describes how it regenerates them.

## `tests/` compares the program with the reference; `test/` holds example inputs

`tests/data/` holds the reference outputs: every combination of mode, strategy and weighting on the example alignments, and further command-line and web cases.
The example alignments are real alignments from Talevich, Mirza & Kannan (2011), *BMC Evolutionary Biology* 11:321.
`test/` holds them under the names that the examples in `README.rst` use, and `tests/data/make_inputs.py` copies them from there into `tests/data/inputs`.

`uv run pytest` or `pixi run test` runs the suite under `tests/`.
`tests/test_cli_cases.py` replays every recorded command-line case that the current command line can run with the installed `cladecompare` and `cladereport` commands, in a scratch directory holding copies of `tests/data/inputs`; the two legacy phospho cases, whose `-nw` option no longer exists, are kept as a record and not replayed.
`tests/test_web.py` submits every recorded form to the cladeweb WSGI application.
Each command or submission must give the recorded return code or status, each written file and successful response must equal the reference byte for byte, and every p-value written to a `.out` file must be finite and within [0, 1].
Where the code departs from the reference by design, the test module declares the difference — a file no longer written, a return code, a different exception, an expected output regenerated into `tests/data/current/` by `tests/data/replay.py`, or a test of its own — and `tests/data/README.md` explains it.
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
- Input that CladeCompare cannot use raises `InputError`, a subclass of `ValueError` defined in `cladecompare/__init__.py`, with a message that names the file and, where it can, the line.
  Each reader of an input file opens it inside `reading`, defined beside `InputError`, which turns an `OSError`, bytes that are not UTF-8 text and, given the form the file should have, a parser's `ValueError` into `InputError`; `core.check_readable` only opens a file that an external program will read, so that a missing or unreadable file is reported before the program runs, while a file the program cannot parse still ends with the program's `RuntimeError`.
  The `main` functions of `cladecompare/cladecompare.py` and `cladecompare/cladereport.py` catch it, printing `<command>: error: <message>` and exiting with status 1; the `-a` parser of the `cladecompare` command turns the `InputError` of `core.parse_alpha` into an argparse usage error.
  `cladeweb` catches it around the whole submission and returns `Error: <message>` with status 200, naming each upload by its file name, as bottle normalizes it, rather than by its temporary path.
  A missing program or MUSCLE 4 or later still raises `RuntimeError`, which, like any other exception, ends a command with a traceback and gives status 500 in `cladeweb`.
- Every module that logs uses `logger = logging.getLogger(__name__)`, and the three commands configure logging through `configure_logging` in `cladecompare/__init__.py`.
- The files under `tests/data/py27/` are generated but committed deliberately, as the reference; regenerate them with `tests/data/capture_py27.py`, for instance through `pixi run -e py27-oracle capture`, rather than editing them.
- The files under `tests/data/current/` are generated from the current code and committed as the expected output of the cases that depart from `tests/data/py27/` by design; regenerate them with `tests/data/replay.py`, review the change against `tests/data/README.md`, and never edit them by hand.
