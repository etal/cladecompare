"""End-to-end tests of the command-line cases that run external programs.

Each case replays a command-line case recorded in ``tests/data/py27/cli``, in
a scratch directory holding copies of ``tests/data/inputs``, with the installed
``cladecompare`` and ``cladereport`` commands in place of the scripts the
record names. The text files it writes must equal the recorded ones; the
reference records the pair-logo PDFs by name only, so those must have the
recorded names and embed the font that CladeCompare bundles.

Every test here carries the ``integration`` mark and is skipped unless the
programs it needs are on ``PATH``. ``pixi run test-integration`` provides them,
runs these tests alone and sets ``CLADECOMPARE_REQUIRE_TOOLS``, which turns
each of those skips into a failure.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "data"

pytestmark = pytest.mark.integration

# Recorded cases that complete with hmmalign or MUSCLE, and the programs each
# needs. pdb.1FIN-AB-CA.gtest is left out: it records a failure that the
# current code still reproduces.
CASES = {
    "hmm.pair.gtest": ["hmmalign"],
    "hmm.single.gtest": ["hmmalign"],
    "pdb.1HCK.gtest": ["hmmalign"],
    "unequal.gtest": ["muscle"],
}

# The reference case unequal.gtest also wrote this file, which the code no
# longer writes; see "Recorded defects" in tests/data/README.md.
NOT_WRITTEN = {"_cc_combined.seq"}


def skip_unless_required(reason):
    """Skip the test, or fail it if CLADECOMPARE_REQUIRE_TOOLS is non-empty."""
    if os.environ.get("CLADECOMPARE_REQUIRE_TOOLS"):
        pytest.fail(reason)
    pytest.skip(reason)


def require(*programs):
    """Skip the test unless every program is on PATH; see skip_unless_required."""
    missing = [name for name in programs if shutil.which(name) is None]
    if missing:
        skip_unless_required("not on PATH: " + ", ".join(missing))


def command(script):
    """The installed console script, beside the interpreter, for a script name."""
    path = Path(sys.executable).with_name(script.removesuffix(".py"))
    assert path.exists(), f"{path.name} is not installed; run 'uv sync' or 'pixi install'"
    return str(path)


@pytest.mark.parametrize("case", sorted(CASES))
def test_reference_case(case, tmp_path):
    require(*CASES[case])
    recorded_dir = DATA / "py27" / "cli" / case
    record = json.loads((recorded_dir / "run.json").read_text(encoding="utf-8"))

    workdir = tmp_path / "work"
    shutil.copytree(DATA / "inputs", workdir)
    inputs = {path.name for path in workdir.iterdir()}
    env = os.environ | {"TMPDIR": str(tmp_path)}
    for step in record["commands"]:
        argv = [command(step["argv"][0]), *step["argv"][1:]]
        result = subprocess.run(argv, cwd=workdir, env=env, capture_output=True,
                                text=True, check=False)
        assert result.returncode == step["returncode"], result.stderr

    written = {path.name for path in workdir.iterdir()} - inputs
    pdfs = {name for name in written if name.endswith(".pdf")}
    assert pdfs == set(record["pdfs"])
    for name in sorted(pdfs):
        assert b"DejaVuSansMono-Bold" in (workdir / name).read_bytes(), name
    expected = set(record["written_files"]) - NOT_WRITTEN
    assert written - pdfs == expected
    for name in sorted(expected):
        assert (workdir / name).read_bytes() == (recorded_dir / name).read_bytes(), name
    leftovers = [path for path in tmp_path.iterdir() if path != workdir]
    assert len(leftovers) == record["leftover_temp_files"], leftovers
