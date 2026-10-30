#!/usr/bin/env python3
"""Write the current output of every command-line case that departs from py27.

Each case of tests/test_cli_cases.py that needs no external program is
replayed with the installed cladecompare and cladereport commands, as that
module replays it.  Where the text files written, their contents or the names
of the pair-logo PDFs differ from the reference in py27/cli, with the
substitutions that module declares applied, the case's directory under
current/ is rewritten: the text files the commands wrote, and a run.json
holding the names of those files ("written_files") and of the PDFs ("pdfs"),
with the meaning those keys have in the reference.  Where they no longer
differ, the directory is removed.  The cases that need hmmalign or MUSCLE are
not replayed, and their directories, of which there are none, are left alone.

The files under current/ are generated; regenerate them with this script and
review the change rather than editing them.  tests/test_cli_cases.py expects
the cases it lists in REGENERATED to match these files and fails if that list
and the directories disagree.  README.md beside this script explains each
regenerated case.

Usage, from the repository root:

    uv run python tests/data/replay.py
"""

import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from test_cli_cases import (  # noqa: E402
    CASES,
    CURRENT,
    REFERENCE,
    expectation,
    load,
    reference_text,
    run_commands,
)


def replay(case):
    """Replay a case; return its current files and record, or None if they
    equal the reference."""
    record = load(case)
    with tempfile.TemporaryDirectory() as tmp:
        workdir, written, results, _tmpdir = run_commands(
            [step["argv"] for step in record["commands"]], Path(tmp)
        )
        codes = [result.returncode for result in results]
        expected_codes = [step["returncode"] for step in record["commands"]]
        if codes != expected_codes:
            sys.exit(f"{case}: return codes {codes}, recorded {expected_codes}")
        pdfs = {name for name in written if name.endswith(".pdf")}
        names = written - pdfs
        files = {name: (workdir / name).read_bytes() for name in sorted(names)}
    if (names, pdfs) == expectation(case, REFERENCE) and all(
        data == reference_text(case, name) for name, data in files.items()
    ):
        return None
    return files, {"written_files": sorted(names), "pdfs": sorted(pdfs)}


def main():
    for case in sorted(CASES):
        if CASES[case]:
            continue
        current = replay(case)
        outdir = CURRENT / case
        if outdir.exists():
            shutil.rmtree(outdir)
        if current is None:
            continue
        files, record = current
        outdir.mkdir(parents=True)
        for name, data in files.items():
            (outdir / name).write_bytes(data)
        (outdir / "run.json").write_text(
            json.dumps(record, indent=2) + "\n", encoding="utf-8"
        )
        print("wrote", outdir.relative_to(HERE))


if __name__ == "__main__":
    main()
