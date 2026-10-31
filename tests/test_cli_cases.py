"""Regression tests that replay the recorded command-line cases.

Each case replays a command-line case recorded in ``tests/data/py27/cli``, in
a scratch directory holding copies of ``tests/data/inputs``, with the installed
``cladecompare`` and ``cladereport`` commands in place of the scripts the
record names. Each command must return the recorded code, and a failing one
must fail with the recorded exception. The text files the commands write must
equal the recorded ones, or the regenerated ones for the cases in
``REGENERATED``; the reference records the pair-logo PDFs by name
only, so those must have the recorded names and embed the font that
CladeCompare bundles. Every p-value written must be finite and within [0, 1].

Where the code departs from the reference by design, as
``tests/data/README.md`` lists, the difference is declared here: a value
substitution in ``SUBSTITUTIONS``, a file in ``NOT_WRITTEN``, a return code
in ``RETURNCODES``, or a case in ``REGENERATED``, whose expected output is the
current output that ``tests/data/replay.py`` writes to ``tests/data/current``.
The cases in ``LEGACY`` are kept as a record and not replayed.

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
CURRENT = DATA / "current"

# Recorded cases replayed as recorded, and the programs each needs.
CASES = {
    "multi.gtest": [],
    "pair.gtest": [],
    "pair.gtest.w": [],
    "pair.hypg": [],
    "pair.hypg.w": [],
    "pair.jsd.a.2": [],
    "pair.jsd.w.a.2": [],
    "pair.phospho": [],
    "pair.phospho.w": [],
    "pair.urn": [],
    "pair.urn.w": [],
    "single.gtest.a.001": [],
    "single.gtest.w.a.001": [],
    "single.hypg.a.2": [],
    "single.hypg.w.a.2": [],
    "single.jsd.a.2": [],
    "single.jsd.w.a.2": [],
    "single.phospho": [],
    "single.phospho.w": [],
    "single.urn": [],
    "single.urn.w": [],
    "hmm.pair.gtest": ["hmmalign"],
    "hmm.single.gtest": ["hmmalign"],
    "pdb.1FIN-AB-CA.gtest": ["hmmalign"],
    "pdb.1HCK.gtest": ["hmmalign"],
    "unequal.gtest": ["muscle"],
}

# Printed p-values that differ from the reference because the reference
# computed them less accurately. SciPy 1.2.1 gave three hypergeometric tails
# of pair.hypg less accurately than later releases do, and exact rational
# arithmetic agrees with the later values; biofrills' chisqprob, which gtest
# called until October 2026, gave one chi-squared tail of multi.gtest less
# accurately than SciPy does. Each (reference, current, counts) replaces a
# p-value printed at the number of sites that `counts` gives for each file.
SUBSTITUTIONS = {
    "multi.gtest": [
        ("0.528722", "0.528721", {"cdc2-b.aln.fasta.out": 1}),
    ],
    "pair.hypg": [
        ("6.8898e-10", "6.88977e-10", {"pair.hypg.out": 5, "pair.hypg.html": 5}),
        ("5.74155e-11", "5.74147e-11", {"pair.hypg.out": 5, "pair.hypg.html": 5}),
        ("2.12435e-09", "2.12434e-09", {"pair.hypg.out": 3, "pair.hypg.html": 3}),
    ],
}

# The reference case unequal.gtest also wrote this file, which the code no
# longer writes.
NOT_WRITTEN = {"_cc_combined.seq"}

# The reference records -s phospho failing with TypeError in every case;
# each now completes ("Recorded defects" in tests/data/README.md).
RETURNCODES = {
    "pair.phospho": [0],
    "pair.phospho.w": [0],
    "single.phospho": [0],
    "single.phospho.w": [0],
}

# Recorded cases whose expected output is the current output in CURRENT,
# written by tests/data/replay.py. Five weighted cases moved with biofrills
# 72282be, whose sequence weights keep their ratios but sum to an estimate of
# the effective number of independent sequences; CladeCompare takes that sum
# as the sample size, while the pseudocounts keep a fixed total, so every
# weighted value that depends on the scale moved. pair.hypg.w also departs
# from its record for the reason that "SciPy version" in tests/data/README.md
# gives. The two urn cases moved because urn now weights the sequences as
# CHAIN does ("CHAIN's weighting in urn" in tests/data/README.md). The
# phospho cases fail in the reference. The three single-alignment gtest cases
# moved because single-mode gtest now scales its expected counts to the
# column total ("The G-test" in tests/data/README.md); single.gtest.w.a.001
# moved for that reason and with the weights.
REGENERATED = {
    "hmm.single.gtest",
    "pair.gtest.w",
    "pair.hypg.w",
    "pair.jsd.w.a.2",
    "pair.phospho",
    "pair.phospho.w",
    "pair.urn.w",
    "single.gtest.a.001",
    "single.gtest.w.a.001",
    "single.hypg.w.a.2",
    "single.phospho",
    "single.phospho.w",
    "single.urn.w",
}

# Recorded at the last commit at which -s phospho ran, with an option that no
# longer exists; they document that behavior and are not replayed.
LEGACY = {"legacy.pair.phospho.nw", "legacy.single.phospho.nw"}


def load(case, directory=REFERENCE):
    text = (directory / case / "run.json").read_text(encoding="utf-8")
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
    """Every p-value in a .out file is finite and within [0, 1]; return them."""
    pvalues = [prob for _fg, _bg, _posn, prob, _stars in report.parse_out(path)]
    invalid = [p for p in pvalues if not (math.isfinite(p) and 0.0 <= p <= 1.0)]
    assert pvalues and not invalid, f"{path.name}: {invalid[:5]}"
    return pvalues


def assert_pdfs(workdir, pdfs, expected):
    assert pdfs == set(expected)
    for name in sorted(pdfs):
        assert b"DejaVuSansMono-Bold" in (workdir / name).read_bytes(), name


def reference_text(case, name):
    """The recorded bytes of a file, with the SUBSTITUTIONS applied."""
    expected = (REFERENCE / case / name).read_bytes()
    for old, new, counts in SUBSTITUTIONS.get(case, []):
        old, new = b"=" + old.encode(), b"=" + new.encode()
        assert expected.count(old) == counts.get(name, 0), (name, old)
        expected = expected.replace(old, new)
    return expected


def expectation(case, directory=None):
    """The names of the text files and the PDFs that a case must write.

    By default the regenerated record in CURRENT for the cases in
    REGENERATED, and the reference otherwise.
    """
    if directory is None:
        directory = CURRENT if case in REGENERATED else REFERENCE
    record = load(case, directory)
    written = set(record["written_files"])
    if directory == REFERENCE:
        written -= NOT_WRITTEN
    return written, set(record["pdfs"])


def expected_text(case, name):
    """The bytes that a case must write to a file."""
    if case in REGENERATED:
        return (CURRENT / case / name).read_bytes()
    return reference_text(case, name)


def expected_codes(case, steps):
    """The return codes that the recorded commands of a case must give."""
    return RETURNCODES.get(case, [step["returncode"] for step in steps])


def check_case(case, tmp_path, path_prefix=None, commands=None):
    """Replay a recorded case and compare it with the record."""
    record = load(case)
    steps = record["commands"]
    workdir, written, results, tmpdir = run_commands(
        commands or [step["argv"] for step in steps], tmp_path, path_prefix
    )

    codes = expected_codes(case, steps)
    assert [r.returncode for r in results] == codes, results[-1].stderr
    for result, step, code in zip(results, steps, codes, strict=True):
        if code:
            error = step["stderr"].strip().splitlines()[-1].split(":")[0]
            assert error in result.stderr, result.stderr

    pdfs = {name for name in written if name.endswith(".pdf")}
    expected, expected_pdfs = expectation(case)
    assert_pdfs(workdir, pdfs, expected_pdfs)
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
    assert recorded == set(CASES) | LEGACY


def test_regenerated_cases_are_declared():
    """CURRENT holds exactly the REGENERATED cases, each unlike its record."""
    assert {path.name for path in CURRENT.iterdir()} == REGENERATED
    assert REGENERATED <= set(CASES)
    assert REGENERATED.isdisjoint(SUBSTITUTIONS)
    for case in sorted(REGENERATED):
        written, pdfs = expectation(case)
        assert (written, pdfs) != expectation(case, REFERENCE) or any(
            expected_text(case, name) != reference_text(case, name) for name in written
        ), case


@pytest.mark.parametrize("strategy", ["urn", "phospho"])
def test_weighting_leaves_single_alignment_unchanged(strategy):
    """urn and phospho weight no sequence of a single alignment.

    A single alignment plays the part of CHAIN's query family, so the
    expected .out and .pttrn files of the weighted case equal those of the
    unweighted one. Together with the replay of the weighted case against
    those expected files, this holds the command to the unweighted output.
    """
    case = f"single.{strategy}"
    for suffix in [".out", ".pttrn"]:
        assert expected_text(case + ".w", case + ".w" + suffix) == (
            expected_text(case, case + suffix)
        )


def test_output_defaults_to_standard_output(tmp_path):
    """Without -o, the table that pair.gtest wrote to a file goes to stdout."""
    argv = load("pair.gtest")["commands"][0]["argv"]
    assert argv[-2] == "-o"
    _workdir, written, [result], _tmpdir = run_commands([argv[:-2]], tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected_text("pair.gtest", "pair.gtest.out").decode()
    assert not any(name.endswith(".out") for name in written)


def test_more_than_two_alignments_ignore_output_options(tmp_path):
    """multi.gtest given -o and -p writes what it recorded, and says so."""
    argv = load("multi.gtest")["commands"][0]["argv"]
    _workdir, written, [result], _tmpdir = run_commands(
        [argv + ["-o", "x.out", "-p", "x.pttrn"]], tmp_path
    )
    assert result.returncode == 0, result.stderr
    assert "ignoring -o and -p" in result.stderr
    assert written - {name for name in written if name.endswith(".pdf")} == set(
        load("multi.gtest")["written_files"]
    )


def test_input_error_prints_a_message_without_traceback(tmp_path):
    _workdir, written, [result], _tmpdir = run_commands(
        [["cladecompare.py", "cdc2-domain.aln.fasta", "cdc2-domain.aln.fasta"]],
        tmp_path,
    )
    assert result.returncode == 1
    assert result.stderr.endswith(
        "cladecompare: error: Background set is a subset of the foreground. "
        "To fix this, use a background containing sequences not in the "
        "foreground.\n"
    )
    assert "Traceback" not in result.stderr
    assert written == set()


@pytest.mark.parametrize(
    "alpha, message",
    [
        ("1.5", "must be between 0 and 1"),
        ("-0.1", "must be between 0 and 1"),
        ("nan", "must be between 0 and 1"),
        ("x", "invalid alpha value"),
    ],
)
def test_rejected_alpha_values_are_usage_errors(alpha, message, tmp_path):
    _workdir, _written, [result], _tmpdir = run_commands(
        [["cladecompare.py", "cdc2-a.aln.fasta", "cdc2-b.aln.fasta", "-a", alpha]],
        tmp_path,
    )
    assert result.returncode == 2
    assert f"argument -a/--alpha: {message}" in result.stderr
    assert "Traceback" not in result.stderr


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
