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

Three entry points sit at the top level and share one library of about 2,100
lines under `cladecomparelib/`:

| script | role |
|---|---|
| `cladecompare.py` | command line; the full option surface |
| `cladeweb.py` | local browser form; one server instance per machine |
| `cladereport.py` | renders an existing `.out` file to HTML, no statistics |

## Five contrast statistics behind one flag

`-s/--strategy` selects one statistics module per run:

- `gtest` — G-test of column amino-acid frequencies
- `urn` — CHAIN's ball-in-urn statistic
- `jsd` — Jensen–Shannon divergence
- `hypg` — hypergeometric
- `phospho` — Ser/Thr/Tyr phosphorylation sites specifically

A strategy is not registered in one place.
Adding a sixth means writing the module under `cladecomparelib/` and then listing it in four more places: the import and the `if`/`elif` chain in `cladecomparelib/core.py` (lines 24 and 76-92), the `choices` tuple and help text of `-s` in `cladecompare.py` (lines 56-65), the strategy table in `cladeweb.py` (line 209, with its import at line 25), and a radio button in that file's form (lines 55-79).
`hypg` shows the cost of missing one: it reached `core.py` and `cladecompare.py` but never the web application, whose form does not offer it and whose strategy table raises `KeyError` on a request that names it.

## `.out` p-values are raw; Benjamini–Hochberg only selects the pattern

The `.out` file holds each column's unadjusted p-value.
`write_pvalues` (`cladecomparelib/core.py:230`) writes the value as the statistic returned it, and its asterisks compare that raw value with `-a/--alpha`, one asterisk per order of magnitude.
`cladereport.py` and the cladeweb result page color each column from the same raw value (`p2class` in `cladecomparelib/report.py`), at fixed cut points that do not depend on alpha.

The multiple-testing correction exists only in `top_hits` (`core.py:431`): a Benjamini–Hochberg step-up procedure at false discovery rate alpha, with m equal to the number of columns in `hits`, including the columns that `process_pair` and `process_one` skip as gaps and assign p = 1.
In paired mode that count is taken after `clean_alignments` has dropped the columns that are gaps in every sequence of both sets; the column numbers in `.out` and `.pttrn` refer to that cleaned alignment.
Its result is the set of significant sites, which becomes the `.pttrn` file when `-p` is given, names the pair logos, and lists the residues highlighted in the PyMOL script.
The step-up search examines only the 50 most significant columns, so at most 50 sites are selected, and where the full procedure would accept more than 50, the result can be smaller than 50.

The consequence for a local change: altering how many columns enter `hits`, for instance with a filter that drops columns before scoring, does not rescale any value in `.out`, but it changes m and therefore every threshold `alpha * k / m` (k being a column's rank by p-value).
A gap-handling rule that only changes which columns receive p = 1 leaves m unchanged but changes the ranking.
Either can move which columns reach the pattern, the logos and the PyMOL script.
The raw values are not purely local either: pseudocounts come from amino-acid frequencies over the whole alignment (`combined_frequencies` in `cladecomparelib/shared.py` in paired mode, `alnutils.aa_frequencies` in single mode), so a change to which sequences or residues enter those frequencies moves the p-value of every column.
Judge such a change on the whole output, pattern included, never on the column you were looking at.

## The library does not currently import

`cladecomparelib/core.py` uses two Biopython APIs that no longer exist:

```
18: from Bio.Alphabet import generic_protein   # Bio.Alphabet removed in 1.78
20: from Bio._py3k import StringIO             # removed in 1.76
```

`generic_protein` is then passed at four call sites (`core.py:276`, `:302`, `:304`, `:345`).
Because the failure is at import time, *nothing* in the package runs on a current Biopython.

A port nevertheless has a reference to be checked against.
`tests/data/py27/` holds the outputs of the unmodified code, captured under Python 2.7 and Biopython 1.76 from commit `e265b1f`, for the command line, `cladereport.py` and the web application.
`tests/data/README.md` records the environment and how to regenerate it, the inputs, and the defects the reference preserves on purpose.
A ported program should reproduce those outputs, or differ from them only by a documented, deliberate change; the statistics are exactly the kind of output where a wrong answer looks entirely plausible.

## biofrills is not the blocker it appears to be

Every statistics module imports [biofrills](https://github.com/etal/biofrills),
last touched in 2013. It ships a Cython extension, `biofrills/cpairutils.pyx`,
which fails to build on Python 3.12 and later — `PyCode_New`'s signature changed
and thread-state fields it uses were removed — and because that repository's
`setup.py` adds the extension to `ext_modules` unconditionally, `pip install
biofrills` fails outright. The surface reading is that a Cython port blocks this
project.

It does not. Measured: this project's entire biofrills surface is four
pure-Python modules — `alnutils`, `consensus`, `logoutils`, and
`stats.chisq.chisqprob` — and the extension is imported in exactly one place,
`biofrills/pairutils.py`, inside a `try`/`except ImportError` guarding an
optional speedup, in a module nothing here reaches. The requirement is "import
four pure-Python modules", not "port Cython".

## `test/` is a smoke harness; `tests/data/` holds the reference

`cd test && make` runs the program in eight configurations — four statistics,
paired and single-alignment — and writes `.out`, `.pttrn` and `.html` files for a
human to look at. **There are no assertions.** A clean `make` means the program
ran to completion, and says nothing about whether the numbers are right. Note
also that the targets pass different `-a` thresholds (`.001`, `.2`, default), so
they are not comparable to each other.

The fixtures are real alignments from Talevich, Mirza & Kannan (2011), *BMC
Evolutionary Biology* 11:321.

`tests/data/` is the committed counterpart: the reference outputs described above, covering the eight Makefile targets and further command-line and web cases.
No test compares the program against it yet.

## Conventions

- The code predates this project's current Python baseline; it is 2013-era
  Python 2/3-straddling style. Modernize deliberately and in separate commits
  from behaviour changes, so a numerical difference is never buried in a
  reformat.
- Output files under `test/` are generated, gitignored, and should stay so.
  The files under `tests/data/py27/` are generated too but committed deliberately, as the reference; regenerate them with `tests/data/capture_py27.py` rather than editing them.
- Maintainer planning for this repository lives in a local beads store that is excluded from git and is not present in a clone.
  If you have `bd`, run `bd prime --mcp`, which renders each memory as a one-line summary because the full text printed by plain `bd prime` can exceed an agent's tool-output limits; `bd memories <keyword>` searches the full text, and `bd recall <key>` prints one memory in full.
  If you do not have `bd`, nothing here depends on it.
