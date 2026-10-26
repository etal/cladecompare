"""Regression tests that replay the recorded command-line cases.

Each case replays a command-line case recorded in ``tests/data/py27/cli``, in
a scratch directory holding copies of ``tests/data/inputs``, with the installed
``cladecompare`` and ``cladereport`` commands in place of the scripts the
record names. Each command must return the recorded code, and a failing one
must fail with the recorded exception. The text files the commands write must
equal the recorded ones; the reference records the pair-logo PDFs by name
only, so those must have the recorded names and embed the font that
CladeCompare bundles. Every p-value written must be finite and within [0, 1].

Where the code departs from the reference by design, as
``tests/data/README.md`` lists, the difference is declared here: a value
substitution in ``SUBSTITUTIONS``, a file in ``NOT_WRITTEN``, or a test of
its own.

The cases that run hmmalign or MUSCLE carry the ``integration`` mark and are
skipped unless the programs are on ``PATH``. ``pixi run test-integration``
provides them, runs these tests alone and sets ``CLADECOMPARE_REQUIRE_TOOLS``,
which turns each of those skips into a failure.
"""

import json
import math
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from cladecomparelib import report

DATA = Path(__file__).resolve().parent / "data"
REFERENCE = DATA / "py27" / "cli"

# Recorded cases replayed as recorded, and the programs each needs.
CASES = {
    "multi.gtest": [],
    "pair.gtest": [],
    "pair.gtest.w": [],
    "pair.hypg": [],
    "pair.jsd.a.2": [],
    "pair.jsd.w.a.2": [],
    "pair.urn": [],
    "pair.urn.w": [],
    "single.gtest.a.001": [],
    "single.gtest.w.a.001": [],
    "single.hypg.a.2": [],
    "single.hypg.w.a.2": [],
    "single.jsd.a.2": [],
    "single.jsd.w.a.2": [],
    "single.urn": [],
    "single.urn.w": [],
    "hmm.pair.gtest": ["hmmalign"],
    "hmm.single.gtest": ["hmmalign"],
    "pdb.1FIN-AB-CA.gtest": ["hmmalign"],
    "pdb.1HCK.gtest": ["hmmalign"],
    "unequal.gtest": ["muscle"],
}

# SciPy 1.2.1, which computed the reference, gave three hypergeometric tails
# of pair.hypg less accurately than later releases do; exact rational
# arithmetic agrees with the later values. Each (reference, current, count)
# replaces a p-value printed at `count` sites, in the .out and .html files.
SUBSTITUTIONS = {
    "pair.hypg": [
        ("6.8898e-10", "6.88977e-10", 5),
        ("5.74155e-11", "5.74147e-11", 5),
        ("2.12435e-09", "2.12434e-09", 3),
    ],
}

# The reference case unequal.gtest also wrote this file, which the code no
# longer writes.
NOT_WRITTEN = {"_cc_combined.seq"}

# Recorded cases with a test of their own below, and the legacy cases, which
# those tests use as their expected output.
OWN_TESTS = {
    "pair.hypg.w",
    "pair.phospho",
    "pair.phospho.w",
    "single.phospho",
    "single.phospho.w",
}
LEGACY = {"legacy.pair.phospho.nw", "legacy.single.phospho.nw"}


def load(case):
    text = (REFERENCE / case / "run.json").read_text(encoding="utf-8")
    return json.loads(text)


def command(script):
    """The installed console script for a recorded script name.

    A record names the script as, for example, ``cladecompare.py``.
    """
    path = Path(sys.executable).with_name(script.removesuffix(".py"))
    assert path.exists(), (
        f"{path.name} is not installed; run 'uv sync' or 'pixi install'"
    )
    return str(path)


def run_commands(commands, tmp_path, path_prefix=None):
    """Run recorded commands in a copy of the inputs, stopping at a failure.

    Return the working directory, the names of the files the commands wrote,
    the completed processes and the temporary directory they used.
    """
    workdir = tmp_path / "work"
    shutil.copytree(DATA / "inputs", workdir)
    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()
    inputs = {path.name for path in workdir.iterdir()}
    env = os.environ | {"TMPDIR": str(tmpdir)}
    if path_prefix:
        env["PATH"] = str(path_prefix) + os.pathsep + env["PATH"]
    results = []
    for argv in commands:
        result = subprocess.run(
            [command(argv[0]), *argv[1:]],
            cwd=workdir,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        results.append(result)
        if result.returncode:
            break
    written = {path.name for path in workdir.iterdir()} - inputs
    return workdir, written, results, tmpdir


def assert_valid_pvalues(path):
    """Every p-value in a .out file is finite and within [0, 1]."""
    pvalues = [prob for _fg, _bg, _posn, prob, _stars in report.parse_out(path)]
    invalid = [p for p in pvalues if not (math.isfinite(p) and 0.0 <= p <= 1.0)]
    assert pvalues and not invalid, f"{path.name}: {invalid[:5]}"


def assert_pdfs(workdir, pdfs, expected):
    assert pdfs == set(expected)
    for name in sorted(pdfs):
        assert b"DejaVuSansMono-Bold" in (workdir / name).read_bytes(), name


def expected_text(case, name):
    expected = (REFERENCE / case / name).read_bytes()
    if name.endswith((".out", ".html")):
        for old, new, count in SUBSTITUTIONS.get(case, []):
            old, new = b"=" + old.encode(), b"=" + new.encode()
            assert expected.count(old) == count, (name, old)
            expected = expected.replace(old, new)
    return expected


def check_case(case, tmp_path, path_prefix=None, commands=None):
    """Replay a recorded case and compare it with the record."""
    record = load(case)
    steps = record["commands"]
    workdir, written, results, tmpdir = run_commands(
        commands or [step["argv"] for step in steps], tmp_path, path_prefix
    )

    assert [r.returncode for r in results] == [s["returncode"] for s in steps], results[
        -1
    ].stderr
    for result, step in zip(results, steps):
        if step["returncode"]:
            error = step["stderr"].strip().splitlines()[-1].split(":")[0]
            assert error in result.stderr, result.stderr

    pdfs = {name for name in written if name.endswith(".pdf")}
    assert_pdfs(workdir, pdfs, record["pdfs"])
    expected = set(record["written_files"]) - NOT_WRITTEN
    assert written - pdfs == expected
    for name in sorted(expected):
        assert (workdir / name).read_bytes() == expected_text(case, name), name
        if name.endswith(".out"):
            assert_valid_pvalues(workdir / name)
    leftovers = list(tmpdir.iterdir())
    assert len(leftovers) == record["leftover_temp_files"], leftovers


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(case, marks=[pytest.mark.integration] if programs else [])
        for case, programs in sorted(CASES.items())
    ],
)
def test_reference_case(case, tmp_path, require):
    require(*CASES[case])
    check_case(case, tmp_path)


def test_every_recorded_case_is_tested():
    recorded = {path.name for path in REFERENCE.iterdir()}
    assert recorded == set(CASES) | OWN_TESTS | LEGACY


@pytest.mark.parametrize("mode", ["pair", "single"])
def test_phospho_matches_legacy(mode, tmp_path):
    """Unweighted phospho reproduces the last commit at which it ran.

    The reference records the TypeError that the alignment size caused
    before it was rounded to an integer; the legacy case records the same
    comparison from before the size became a floating-point sum.
    """
    case = mode + ".phospho"
    record = load(case)
    workdir, written, results, _tmpdir = run_commands(
        [step["argv"] for step in record["commands"]], tmp_path
    )
    assert [r.returncode for r in results] == [0]
    legacy = f"legacy.{mode}.phospho.nw"
    expected = load(legacy)["pdfs"]
    assert_pdfs(
        workdir,
        {name for name in written if name.endswith(".pdf")},
        [name.replace(legacy, case) for name in expected],
    )
    for suffix in [".out", ".pttrn"]:
        assert (workdir / (case + suffix)).read_bytes() == (
            REFERENCE / legacy / (legacy + suffix)
        ).read_bytes()
    assert_valid_pvalues(workdir / (case + ".out"))


@pytest.mark.parametrize("mode", ["pair", "single"])
def test_phospho_weighted_completes(mode, tmp_path):
    """Weighted phospho, which has never run before, completes.

    Its values have no reference, and the rule for its sample size and
    trial count is not settled, so only their range is checked.
    """
    case = mode + ".phospho.w"
    workdir, _written, results, _tmpdir = run_commands(
        [step["argv"] for step in load(case)["commands"]], tmp_path
    )
    assert [r.returncode for r in results] == [0]
    assert_valid_pvalues(workdir / (case + ".out"))


def test_hypg_weighted(tmp_path):
    """Weighted hypg gives finite p-values, and on these inputs no pattern.

    The reference values came from SciPy 1.2.1, which accepted a non-integral
    population count, and have no counterpart in an integer model. With the
    counts rounded up, 6 sequences are drawn from 13, and the smallest
    attainable p-value, 1 / C(13, 6) or about 5.8e-4, is attained by only 16
    columns, too few to pass the Benjamini-Hochberg threshold at the default
    alpha, which first reaches that value at rank 29; every larger p-value
    exceeds the threshold at its rank.
    """
    case = "pair.hypg.w"
    workdir, written, results, _tmpdir = run_commands(
        [step["argv"] for step in load(case)["commands"]], tmp_path
    )
    assert [r.returncode for r in results] == [0, 0], results[-1].stderr
    assert_valid_pvalues(workdir / (case + ".out"))
    assert (workdir / (case + ".pttrn")).read_text() == "1:"
    assert written == {case + suffix for suffix in [".out", ".pttrn", ".html"]}


RUN_GAPS = '''\
#!{python}
"""Stand-in for MAPGAPS run_gaps: run_gaps PROFILE FASTA.

Aligns FASTA to PROFILE.hmm with hmmalign and writes FASTA_aln.cma, one
{{()ROW()}}* record per sequence, with the inserts in lowercase.
"""
import sys
from io import StringIO

from Bio import SeqIO

from cladecomparelib.core import hmmalign_a2m

profile, fasta = sys.argv[1:]
with open(fasta + "_aln.cma", "w") as out:
    a2m = hmmalign_a2m(profile + ".hmm", fasta)
    for rec in SeqIO.parse(StringIO(a2m), "fasta"):
        row = str(rec.seq).replace(".", "")
        out.write(">%s\\n{{()%s()}}*\\n" % (rec.id, row))
'''


@pytest.mark.integration
def test_mapgaps_matches_hmm_alignment(tmp_path, require):
    """--mapgaps reads run_gaps output as --hmm reads hmmalign output.

    MAPGAPS is not packaged, so a stand-in run_gaps aligns with hmmalign and
    writes CMA records; the outputs must then equal those of the recorded
    --hmm case, PyMOL script included.
    """
    require("hmmalign")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    run_gaps = bindir / "run_gaps"
    run_gaps.write_text(RUN_GAPS.format(python=sys.executable))
    run_gaps.chmod(run_gaps.stat().st_mode | stat.S_IXUSR)

    case = "pdb.1HCK.gtest"
    commands = [step["argv"] for step in load(case)["commands"]]
    hmm = commands[0].index("--hmm")
    commands[0][hmm : hmm + 2] = ["--mapgaps", "CDK"]
    check_case(case, tmp_path, path_prefix=bindir, commands=commands)
