# Reference outputs from the Python 2 code base

The files under `py27/` record what CladeCompare computed on Python 2.7, before the port to Python 3.
They are the reference against which the port is checked: a ported program should reproduce these outputs, or the difference should be explained by a documented, deliberate change.

## Producing versions

`py27/versions.json` is written by the capture script and is the authoritative record.
At capture time it read:

| Component | Version |
|---|---|
| Python | 2.7.15 (CPython, conda-forge build) |
| Platform | Darwin 25.6.0, x86_64 under Rosetta 2 on Apple silicon |
| Biopython | 1.76 |
| NumPy | 1.16.5 |
| SciPy | 1.2.1 |
| ReportLab | 3.5.42 |
| bottle | 0.12.18 |
| WebLogo | 3.6.0 (from PyPI; provides the `weblogolib` and `corebio` modules) |
| Pillow | 6.2.1 |
| HMMER | 3.0 (`hmmalign`) |
| MUSCLE | 3.8.1551 |
| cladecompare | commit `e265b1f`, the last commit to change code before this reference was captured, extracted from git |
| cladecompare, legacy cases | commit `44c81cc`, the parent of the commit that broke `-s phospho` |
| biofrills | commit `36684bb`, its last Python 2 commit |

The capture never runs the working tree.
It extracts each commit above from git into a scratch directory and runs that, so the reference can be regenerated after the source has been ported.

The unmodified code needs HMMER 3.0, because it calls `hmmalign --allcol`, and the later bioconda releases checked (3.1b2, 3.2.1, 3.3.2 and 3.4) all lack that option.
MUSCLE 3.8 is needed for the same reason: MUSCLE 5.3 rejects the `-profile` option that the code uses.
The current code does not need these particular releases; the external-aligner section of `AGENTS.md` states what it requires.

## Regenerating

conda-forge publishes no Python 2.7 build for osx-arm64, but its osx-64 build runs under Rosetta, and bioconda has osx-64 builds of both tools.
The `py27-oracle` environment in `pyproject.toml` holds the versions of the table above, and pixi installs it as `pixi.lock` records it.
It differs from the environment of the original capture in two respects: that environment installed WebLogo 3.6.0 from PyPI with pip, while pixi takes the same release from bioconda, because pixi resolves PyPI packages with uv, which cannot use Python 2.7; and it names setuptools explicitly, whose `pkg_resources` the capture uses to read the package versions and which the original environment already contained.

The reference was captured on macOS only, and the environment is defined for osx-64 only.
The capture extracts commits with `git archive`, so this repository must contain `e265b1f` and `44c81cc`, and the biofrills repository `36684bb`; a clone of depth 1 of either contains none of them.
Then, with the biofrills repository beside this one (or named by `$BIOFRILLS_REPO`):

```sh
python3 tests/data/make_inputs.py          # optional: rewrites inputs/
pixi run -e py27-oracle capture
```

On the host of the original capture, this rewrites every file under `py27/` byte for byte, `versions.json` included.

`make_inputs.py` needs Python 3.6 or later.
It copies the test files from `test/`, derives the aligned FASTA inputs from them, and downloads two structures from the RCSB Protein Data Bank; rerunning it reproduces the committed files byte for byte unless the RCSB has since remediated those entries.
`capture_py27.py` puts the interpreter's own `bin/` first on `PATH`, so that `hmmalign` and `muscle` come from the same environment, checks both tool versions, and then deletes and rewrites `py27/`.
A capture that fails after that point leaves `py27/` deleted or incomplete; `git checkout -- tests/data/py27` restores the committed files, and `git clean -fd tests/data/py27` removes any that the capture had added.
It passes `PYTHONHASHSEED` through to every command and records it in `versions.json`; the committed reference was captured with it unset.
Two consecutive captures produce identical files, and captures with `PYTHONHASHSEED` set to 1, 2 and 3 differ from the committed reference only in that record.

## Inputs

| File | Source | Used for |
|---|---|---|
| `scttlre-domain.fasta_aln.cma`, `cdc2-domain.fasta_aln.cma` | `test/`, unchanged | paired and single-alignment cases; a web submission of CMA files; `scttlre-domain.fasta_aln.cma` is also one of the three alignments of the three-alignment mode |
| `scttlre-domain.fasta`, `cdc2-domain.fasta` | `test/`, unchanged (unaligned) | `--hmm` and `--pdb` cases |
| `CDK.hmm` | `test/`, unchanged (HMMER 3 profile `CDK_CMGC`) | `--hmm` and `--pdb` cases |
| `scttlre-domain.aln.fasta`, `cdc2-domain.aln.fasta` | the two CMA files, reduced to their consensus columns | web submissions; `scttlre-domain.aln.fasta` is also the foreground of the unequal-width case |
| `cdc2-a.aln.fasta`, `cdc2-b.aln.fasta` | the first 18 and the remaining 19 records of `cdc2-domain.aln.fasta` | three-alignment mode |
| `cdc2-trimmed.aln.fasta` | `cdc2-domain.aln.fasta` without columns 121 to 130 | alignments of unequal width |
| `1HCK.pdb` | RCSB entry 1HCK (CDK2 with ATP, one chain), unchanged | `--pdb` |
| `1FIN-AB-CA.pdb` | RCSB entry 1FIN reduced to its HEADER line and the C-alpha atoms of chains A (CDK2) and B (cyclin A) | `--pdb` with more than one chain |

`1HCK.pdb` is kept complete, with its ATP `HETATM`, `SEQRES` and `CONECT` records, because the ported code reads it with a later Biopython release than the reference did, and that reader should be checked on a full entry; reduced to its C-alpha atoms, it would give the same outputs here.

The `.aln.fasta` files apply the rule of `cladecomparelib.core.cma_blocks`: keep the aligned part of each CMA sequence line and drop its lowercase insert residues.
In these files the aligned part always lies between a leading `{()` and a trailing `()}*`.
They hold the same alignments that CladeCompare computes from the CMA files.

## Layout of `py27/`

### `cli/<case>/`

Each command-line case runs one or more commands in a scratch directory that holds copies of `inputs/`.
`run.json` records the source tree that ran the case (`source`, a key of the `commits` record in `versions.json`) and, for each command, the arguments (`argv`, with the script named relative to the source tree), the return code, standard output and standard error.
A failing command ends the case, and its failure is part of the reference.
`written_files` lists the text files the commands wrote, which sit beside `run.json`; `pdfs` lists the pair-logo PDFs by name only, since they embed their creation time; `leftover_temp_files` counts the files left in the case's temporary directory.

Case names have the form `<mode>.<strategy>[.w][.a<alpha>]`: `pair` compares `scttlre-domain` (foreground) with `cdc2-domain` (background), `single` analyzes `cdc2-domain` alone, `.w` adds `-w`, and `.a<alpha>` gives `-a`.
Every combination of mode, strategy and weighting is recorded, with the `-a` value that `test/Makefile` used for that mode and strategy; that file was an unasserted smoke test of eight configurations, since replaced by the test suite.
Each successful run is followed by `cladereport.py` on its output, as the Makefile did.
The eight Makefile targets are marked in `makefile_target`:

| Makefile target | Case |
|---|---|
| `scttlre-v-cdc2-gtest` | `pair.gtest.w` |
| `scttlre-v-cdc2-urn` | `pair.urn` |
| `scttlre-v-cdc2-jsd` | `pair.jsd.a.2` |
| `scttlre-v-cdc2-hypg` | `pair.hypg` |
| `cdc2-gtest` | `single.gtest.w.a.001` |
| `cdc2-urn` | `single.urn` |
| `cdc2-jsd` | `single.jsd.a.2` |
| `cdc2-hypg` | `single.hypg.a.2` |

Their `.out`, `.pttrn` and `.html` files were byte-identical to those of an independent `make` in `test/` with the same interpreter and packages, apart from the output file names, which appear in the HTML title and heading.

The remaining cases:

| Case | What it runs |
|---|---|
| `multi.gtest` | three alignments, each in turn the foreground against the other two |
| `unequal.gtest` | aligned FASTA inputs of different widths, which CladeCompare aligns to each other with MUSCLE |
| `hmm.pair.gtest`, `hmm.single.gtest` | unaligned sequences aligned with `hmmalign` to `CDK.hmm` |
| `pdb.1HCK.gtest` | as `hmm.pair.gtest`, plus `--pdb 1HCK.pdb`, which writes the PyMOL script `1HCK.pdb.pml` |
| `pdb.1FIN-AB-CA.gtest` | the same with a two-chain structure |
| `legacy.pair.phospho.nw`, `legacy.single.phospho.nw` | `-s phospho -nw -a 0.005` at the legacy commit (`source` is `cladecompare_legacy`; the current code has no `-nw` option); that commit's `cladereport.py` has no `-o` option, so these cases end after `cladecompare.py` |

### `web/<case>/`

Each web case submits the cladeweb form to its WSGI application in the capture process, with the fields and uploads recorded in `run.json`.
`run.json` also records the response status, the traceback the application logged (`error`), its log output in the format `cladeweb.py` configures when run as a server (`log`), what it printed to standard output (`stdout`, one line per temporary file it removed), the files written to the working directory (`written_files`, empty in every captured case), and the number of files left in the temporary directory (`leftover_temp_files`); `response.html` holds the response body.
The form offered gtest, urn, jsd and phospho; `web.pair.hypg` submits the one strategy the command line accepted and the form did not offer.

### Normalization

Standard output, standard error, web response bodies, logs and tracebacks have the scratch working directory, the temporary directory, the extracted cladecompare and biofrills trees and the Python installation prefix replaced by `<workdir>`, `<tmp>`, `<cladecompare>` (or `<cladecompare-legacy>`), `<biofrills>` and `<python-prefix>`.
The random part of a temporary file name is replaced by `tmpXXXXXX`.

## Regenerated outputs in `current/`

`current/<case>/` holds the expected output of each command-line case whose output now differs from `py27/` by a deliberate change that a substitution of a few values cannot express.
It holds the text files the case writes and a `run.json` with only the two fields that can move, `written_files` and `pdfs`, which mean what they mean in `py27/`; the commands, return codes and the count of leftover temporary files remain those of the reference.
`tests/test_cli_cases.py` lists these cases in `REGENERATED`, and fails if that list and the directories disagree or if a regenerated case no longer differs from its reference.

The files are generated from the code under test, so each regeneration is a claim to be reviewed against the explanation below, never a way to make a failing case pass.
With the dependencies of the lock files installed, `replay.py` replays every case that needs no external program and rewrites `current/` for those whose output differs from the reference:

```sh
uv run python tests/data/replay.py
```

Its output is deterministic, so a rerun on unchanged code leaves `current/` unchanged.
It exits with an error while `PYTHONPATH` is set, because the commands it runs would inherit that path, and the files must record the biofrills that the lock files pin rather than a checkout on the path.
The tests inherit `PYTHONPATH`, so they import a checkout on that path.

- **Scale of the sequence weights** (`pair.gtest.w`, `pair.hypg.w`, `pair.jsd.w.a.2`, `single.gtest.w.a.001`, `single.hypg.w.a.2`).
  biofrills `72282be` changed the `none` scaling of `sequence_weights`, which `-w` uses, so that the weights of an alignment sum to an estimate of its effective number of independent sequences; the ratios between the weights did not change.
  The sums moved from 5.004 to 2.117 for the 11 `scttlre-domain` sequences of the paired cases, from 6.146 to 7.734 for the 36 `cdc2-domain` sequences of the paired cases, and from 6.317 to 7.734 for the 37 `cdc2-domain` sequences of the single-alignment cases; the sequence that the paired cases remove from the background is identical to another, so it no longer changes the sum.
  CladeCompare uses the weights as counts and their sum as the sample size, while the pseudocounts have a fixed total of 0.5, so the change of scale moves every weighted value that depends on it.
  The change is exactly this rescaling: the code at biofrills `ff11246`, with each weight vector multiplied to give the new sum, writes every file of these cases byte for byte.
  The pattern of `pair.gtest.w` falls from 9 sites to none, so its six pair logos are no longer written; `pair.jsd.w.a.2` selects 32 sites instead of 29, in 7 logos instead of 8; `single.gtest.w.a.001` selects 50 sites instead of 41; the pattern of `single.hypg.w.a.2` remains empty.
  For `pair.hypg.w`, see "SciPy version" below.
  `single.jsd.w.a.2` is the one weighted case that still reproduces the reference: single-alignment jsd normalizes the weighted counts of a column and adds no pseudocounts, so the scale cancels.
  The weighted phospho cases, whose reference records a failure, have tests of their own.
- **CHAIN's weighting in urn** (`pair.urn.w`, `single.urn.w`).
  Under `-w`, urn now weights the sequences as CHAIN does (Neuwald et al. 2003, Methods): the foreground, CHAIN's query family, is not weighted, and the background is weighted with CHAIN's robust weights, which the `andy` scaling of `sequence_weights` gives from biofrills `205e1d1` on: normalized to a mean of one-half, with weights greater than one truncated to one.
  The reference instead weighted both sets with the `none` scaling while taking the number of foreground sequences as the number of trials, so that the foreground count was understated against it.
  These outputs no longer depend on the `none` scaling.
  `pair.urn.w` had no p-value below 0.06 and an empty pattern; its smallest p-value is now 2.08e-05, and it selects 24 sites, in 8 pair logos.
  `single.urn.w` analyzes CHAIN's query family alone, which is not weighted, so it writes the `.out` and `.pttrn` files that the reference records for `single.urn`, and its pattern has 50 sites instead of 3.

Every unweighted case, including those that run hmmalign and MUSCLE, reproduces the reference as before.

## Recorded defects

The reference records the original behavior, defects included.
A port that fixes one of them differs from the reference by design, and the test should state the expected difference.

- `-s phospho` fails with `TypeError` in every case captured from `e265b1f`, paired and single, with and without `-w`, and in the web form: the alignment size reaches `range()` as a float.
  It ran only without weights, before the unweighted size also became a floating-point sum; the legacy cases record that behavior.
  The code has since been fixed to round the size up to an integer, so it differs from these cases by design: the unweighted cases now match `legacy.pair.phospho.nw` and `legacy.single.phospho.nw`, and the weighted and web cases complete where the reference records the failure.
- `pdb.1FIN-AB-CA.gtest` fails with `KeyError` in `choose_best_aligned`, which indexes the aligned chains with a (key, value) pair whenever more than one chain aligns to the profile, so every multi-chain structure fails.
- `multi.gtest` writes its three `.out` and `.pttrn` pairs, but `cladereport.py` on more than one `.out` file fails with `ValueError` in `report.table_multi`.
- `web.pair.hypg` returns status 500 with `KeyError`: the web application's strategy table lacks hypg.
  The form and the web application now take their strategies from the same registry as the command line, so this case differs from the reference by design: it returns status 200, and its report equals `cli/pair.hypg/pair.hypg.html`, with the same SciPy substitutions, apart from the title.
- `web.pdb.1HCK.gtest` returns status 500 with `AssertionError` from `pmlscript.pdbid_from_fname`: the uploaded structure is saved under a temporary name without the `.pdb` suffix that the function asserts.
  The function now raises `InputError` for a file name it cannot read, so that the check also runs under `python -O`, and this case fails with `InputError` instead, still with status 500.
- Every paired web submission leaves files in the temporary directory (`leftover_temp_files`, between 2 and 17 per case), as does `web.single.phospho`, which fails in the reference; the successful single-alignment submissions leave none.
  For the successful submissions, those files are exactly the pair-logo PDFs named in the log, which are written beside the temporary pattern file and never served.
  The failing submissions also leave the uploads and whatever other temporary files the request had already created, because the exception ends the request before the cleanup that `stdout` records for the successful ones; in the reference, `web.pair.hypg` fails before the output files are opened, so it leaves only its two uploads.
  Each request now writes its uploads and outputs, the pair-logo PDFs included, into one temporary directory that is removed when the request ends, whether it succeeds, returns an error message or raises, so this count differs from the reference by design: no submission leaves a file of its own, and the `Cleaned up` lines in `stdout` are gone.
  The log names those files by fixed names inside that directory, such as `seqfile1.seq` and `cladecompare-<sites>.pdf`, where this reference has `tmpXXXXXX` names.
  The files that bottle 0.13 spools (see "Python 3, bottle and Biopython" below) are bottle's own, outside that directory.
- `single.hypg` cases record the fraction of the column that matches the consensus, not a hypergeometric probability.
  A column whose consensus is X raises `KeyError`, as in the reference code; no recorded case has such a column.
- `unequal.gtest` writes `_cc_combined.seq` into the working directory.
  The code no longer writes that file, so this case differs from the reference by design: `_cc_combined.seq` is absent from the written files, and the log line `Wrote _cc_combined.seq` is gone.

## Interpreting a mismatch

- **SciPy version.**
  The reference computed hypg as `1 - hypergeom.cdf(...)`, which loses precision for small probabilities, and SciPy releases differ there: under SciPy 1.11.4 and 1.18.1 alike, 13 of the small p-values in `pair.hypg` differ from this reference in the sixth significant digit.
  The newer values are the accurate ones: exact rational arithmetic for the three hypergeometric tails behind those 13 sites gives 6.88977e-10, 5.74147e-11 and 2.12434e-09, as both newer releases do, where this reference has 6.8898e-10, 5.74155e-11 and 2.12435e-09.
  The code now computes the tail directly with `hypergeom.sf`: on the 242 columns of `pair.hypg` under SciPy 1.18.1, its largest relative error against exact arithmetic is 6e-16, against 2e-7 for `1 - cdf`, and at the six significant digits of the `.out` file both print the newer values, which `tests/test_cli_cases.py` substitutes for the recorded ones.
  The values recorded for `pair.hypg.w` come from SciPy 1.2.1, which accepted the non-integral number of successes that weighted counts supply; SciPy 1.11.4 and 1.18.1 return `nan` for it.
  In 85 of the 242 columns that number is smaller than the foreground count rounded up, so those recorded values are probabilities of drawing more successes than the population holds.
  The code now rounds the number of successes up, as it rounds the other counts, so `pair.hypg.w` differs from this reference by design, and its expected output is regenerated (see "Regenerated outputs in `current/`" above): its p-values are finite, and its pattern, which records 16 sites, is empty.
  The effective foreground and background sizes, 2.117 and 7.734 sequences, each round up separately, to 3 and 8, so 3 sequences are drawn from 11, and across the 242 columns only 24 distinct p-values occur; the smallest, 1 / C(11, 3) or about 6.1e-3, is attained by 11 columns and already exceeds the default alpha of 0.005, so the Benjamini-Hochberg procedure selects no site at any rank.
  (`single.hypg.w.a.2` does not call `hypergeom`; see the recorded defects.)
- **Consensus ties.**
  The consensus residue of a tied column comes from Python 2 dictionary order in biofrills.
  Every case was also captured with `PYTHONHASHSEED` set to 1, 2 and 3, with identical outputs each time, so tie-breaking does not affect these inputs; a different consensus letter on other data may originate in biofrills' own tie-breaking rule rather than in CladeCompare.
  The ported biofrills breaks ties by an explicit rule instead, and gives the same consensus as `36684bb` on the alignments that the biofrills section of `AGENTS.md` lists as measured.
- **biofrills version.**
  The ported CladeCompare runs against the ported biofrills rather than `36684bb`.
  On the alignments that the biofrills section of `AGENTS.md` lists as measured, the values biofrills returns to CladeCompare are identical or differ by at most a few units in the last place.
  End to end, with `ff11246` and the code of that time, every output file of every command-line case and the response of every successful web submission reproduced this reference, including `multi.gtest` and the cases that align with hmmalign and MUSCLE, under HMMER 3.0 and 3.4 alike, apart from the SciPy and hypg differences above, the deliberate changes to `-s phospho`, `unequal.gtest`, `web.pair.hypg` and `web.pdb.1HCK.gtest` under "Recorded defects", and the differences listed under "Python 3, bottle and Biopython" below.
  With `4aa8b57`, the commit now pinned, the same holds for every unweighted case, measured with HMMER 3.4, and the seven weighted cases listed under "Regenerated outputs in `current/`" above differ as that section describes.
- **Raw p-values.**
  The `.out` files hold unadjusted p-values, and their stars compare those raw values with alpha.
  The Benjamini-Hochberg procedure, limited to the 50 most significant sites, only selects the sites in the `.pttrn` file, which also name the pair logos and the PyMOL script.
- **Log format.**
  Standard error and the web logs come from logging with the format `%(module)s [@%(lineno)s]`, so editing a module above a logging call changes them without changing behavior.
  The legacy cases log in the default format of the root logger instead (`WARNING:root:...`) and record no `INFO` lines.
  Three messages of the current code also differ in wording: the duplicate-sequence message names the set (`Purging N identical sequences from the foreground` or `from the background`, where this reference has `from the alignment`), and the strategy lines read `Using ball-in-urn model of consensus residue conservation` and `Using urn model for phosphorylatable residues (S/T/Y)`, where this reference has `Using ball-in-urn statistical model` and `Using urn model for phosphorylatable residues`.
  One message is new: with `-w`, `single.urn.w` logs `This strategy does not weight the sequences of one alignment`, which the reference does not record.
  At both captured commits, `pairlogo.py` looked for the DejaVu font first under `/usr/share/fonts/truetype/ttf-dejavu/` and then for a bundled copy, and used Courier if it found neither.
  `e265b1f` added that copy as `data/DejaVuSansMono-Bold.ttf` and pointed the fallback at it; at the legacy commit the fallback pointed inside `cladecomparelib/`, where no font was ever committed.
  On the capture host the system path did not exist either, so at the legacy commit `pairlogo.py` warned while being imported, and that first call configured logging before `cladecompare.py` could.
  On a host that has the font at that system path, the legacy standard error would differ from this reference; at `e265b1f` a font is always found, at least the bundled copy, so the current cases do not depend on the host in this way.
  The current code ships the same font file as package data and always loads that copy, with no system path and no fallback.
- **Python 3, bottle and Biopython.**
  Under Python 3, tracebacks carry Python 3's own formatting and exception text; `web.pair.hypg`, for instance, logged `KeyError: 'hypg'` under Python 3, before the form offered hypg, where this reference has `KeyError: u'hypg'`.
  bottle 0.13 words its error page slightly differently from 0.12 (`font-family: sans-serif` for `sans`), which changes the response body of the failing web cases.
  bottle 0.13 also spools a request body or upload larger than 100 KiB to a named temporary file, which remains in the temporary directory while the request is still referenced, so `leftover_temp_files` is larger by one for the request and one for each such upload (`CDK.hmm`, `1HCK.pdb`) when counted inside the capture process; bottle 0.12.25 gives the recorded counts.
  Biopython 1.83 and 1.88 warn about the `END` record of `1FIN-AB-CA.pdb` (`PDBConstructionWarning`), which adds two lines to the standard error of `pdb.1FIN-AB-CA.gtest`.
- **Not captured.**
  MAPGAPS (`--mapgaps`): `run_gaps` is not available from conda-forge or bioconda.
  The pair-logo PDFs, beyond their names.
  The `--pmlout` option, which only renames the PyMOL script.
  The form page itself (`GET /cladecompare`), which the code generates from the strategy registry; `tests/test_web.py` checks its strategy buttons.
