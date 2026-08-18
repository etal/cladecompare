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

`-s/--strategy` selects one module, and adding a sixth means adding a module
rather than a branch (`cladecompare.py:56-58`):

- `gtest` — G-test of column amino-acid frequencies
- `urn` — CHAIN's ball-in-urn statistic
- `jsd` — Jensen–Shannon divergence
- `hypg` — hypergeometric
- `phospho` — Ser/Thr/Tyr phosphorylation sites specifically

**p-values are FDR-corrected across the whole alignment**, Benjamini–Hochberg
step-up, in `cladecomparelib/core.py:440`. The consequence is easy to miss when
changing something local: anything that alters how many columns reach the
correction — a filter, a gap-handling rule, a change to which columns are scored
— moves the reported p-value of *every* column, not only the ones it touched.
Judge such a change on the whole output, never on the column you were looking at.

## The library does not currently import

`cladecomparelib/core.py` uses two Biopython APIs that no longer exist:

```
18: from Bio.Alphabet import generic_protein   # Bio.Alphabet removed in 1.78
20: from Bio._py3k import StringIO             # removed in 1.76
```

`generic_protein` is then passed at four call sites (`core.py:276`, `:302`,
`:304`, `:345`). Because the failure is at import time, *nothing* in the package
runs on a current Biopython, and there is therefore no working baseline on this
machine to compare a port against.

That ordering matters more than the port, which is small. **Capture reference
outputs from an environment where the code does run, before changing any
behaviour.** Without them a port has no oracle, and the statistics are exactly
the kind of output where a wrong answer looks entirely plausible.

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

## `test/` is a smoke harness, not a regression suite

`cd test && make` runs the program in eight configurations — four statistics,
paired and single-alignment — and writes `.out`, `.pttrn` and `.html` files for a
human to look at. **There are no assertions.** A clean `make` means the program
ran to completion, and says nothing about whether the numbers are right. Note
also that the targets pass different `-a` thresholds (`.001`, `.2`, default), so
they are not comparable to each other.

The fixtures are real alignments from Talevich, Mirza & Kannan (2011), *BMC
Evolutionary Biology* 11:321.

## Conventions

- The code predates this project's current Python baseline; it is 2013-era
  Python 2/3-straddling style. Modernize deliberately and in separate commits
  from behaviour changes, so a numerical difference is never buried in a
  reformat.
- Output files under `test/` are generated. They are gitignored and should stay
  that way.
- Maintainer planning for this repository lives in a local beads store that is
  excluded from git and is not present in a clone. If you have `bd`, run
  `bd prime`; if you do not, nothing here depends on it.
