#!/usr/bin/env python3
"""Write the current output of every command-line case that departs from py27.

Each case of tests/test_cli_cases.py is replayed with the installed
cladecompare and cladereport commands, as that module replays it, and must
give the return codes that the module expects.
Where the text files written, their contents or the names
of the pair-logo PDFs differ from the reference in py27/cli, with the
substitutions that module declares applied, the case's directory under
current/ is rewritten: the text files the commands wrote, and a run.json
holding the names of those files ("written_files") and of the PDFs ("pdfs"),
with the meaning those keys have in the reference.  Where they no longer
differ, the directory is removed.  A case that needs hmmalign or MUSCLE is
replayed only when the program is on PATH, as it is under ``pixi run``;
otherwise the script names the case and leaves its directory alone.

The files under current/ are generated; regenerate them with this script and
review the change rather than editing them.  tests/test_cli_cases.py expects
the cases it lists in REGENERATED to match these files and fails if that list
and the directories disagree.  README.md beside this script explains each
regenerated case.  The script refuses to run while PYTHONPATH is set, which
the commands it runs would inherit, so that the files record the packages
that the lock files install, biofrills at its pinned commit among them, and
not a checkout on that path.

Usage, from the repository root:

    pixi run python tests/data/replay.py

or ``uv run python tests/data/replay.py``, which skips the cases that need
hmmalign or MUSCLE unless those are installed.
"""

import json
import os
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
    expected_codes,
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
        expected = expected_codes(case, record["commands"])
        if codes != expected:
            sys.exit(f"{case}: return codes {codes}, expected {expected}")
        pdfs = {name for name in written if name.endswith(".pdf")}
        names = written - pdfs
        files = {name: (workdir / name).read_bytes() for name in sorted(names)}
    if (names, pdfs) == expectation(case, REFERENCE) and all(
        data == reference_text(case, name) for name, data in files.items()
    ):
        return None
    return files, {"written_files": sorted(names), "pdfs": sorted(pdfs)}


def main():
    if os.environ.get("PYTHONPATH"):
        sys.exit("PYTHONPATH is set; unset it so that the installed packages run")
    for case in sorted(CASES):
        missing = [program for program in CASES[case] if not shutil.which(program)]
        if missing:
            print("skipped", case, "for want of", ", ".join(missing), file=sys.stderr)
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
