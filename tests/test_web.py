"""Regression tests that replay the recorded cladeweb form submissions.

Each case submits the form fields and uploads recorded in
``tests/data/py27/web/<case>/run.json`` to the cladeweb WSGI application,
with the working and temporary directories redirected to scratch
directories. The response status must be the recorded one; a successful
response must equal the recorded body, and a failed one must have logged the
recorded exception, or the one that ``ERRORS`` names. No case may write into
the working directory.

The cases that run hmmalign carry the ``integration`` mark; see
``test_cli_cases.py``.
"""

import io
import json
import os
import re
import tempfile
import wsgiref.util
from pathlib import Path

import bottle
import pytest
from test_cli_cases import expected_text

import cladecomparelib.cladeweb  # registers the routes
from cladecomparelib import core

DATA = Path(__file__).resolve().parent / "data"
REFERENCE = DATA / "py27" / "web"

# Recorded submissions replayed as recorded, and the programs each needs.
CASES = {
    "web.hmm.pair.gtest": ["hmmalign"],
    "web.pair-cma.gtest": [],
    "web.pair.gtest": [],
    "web.pair.jsd": [],
    "web.pair.urn": [],
    "web.pdb.1HCK.gtest": ["hmmalign"],
    "web.single.gtest": [],
    "web.single.jsd": [],
    "web.single.urn": [],
}

# bottle 0.13 spools a request body or upload larger than 100 KiB to a named
# temporary file, which is still present when the files are counted, so the
# count of files left in the temporary directory is not comparable for the
# submissions that upload CDK.hmm.
SPOOLED = {"web.hmm.pair.gtest", "web.pdb.1HCK.gtest"}

# Recorded failures that now raise a different exception by design, and the
# exception each raises: pmlscript reports a file name it cannot read with
# InputError instead of an assertion.
ERRORS = {"web.pdb.1HCK.gtest": "InputError"}

# Recorded submissions with a test of their own below.
OWN_TESTS = {"web.pair.hypg", "web.pair.phospho", "web.single.phospho"}

BOUNDARY = "cladecompare-test-boundary"


def multipart(fields, files):
    parts = []
    for name, value in fields:
        parts.append(
            f"--{BOUNDARY}\r\nContent-Disposition: form-data; "
            f'name="{name}"\r\n\r\n{value}\r\n'.encode()
        )
    for name, fname in files:
        parts.append(
            f"--{BOUNDARY}\r\nContent-Disposition: form-data; "
            f'name="{name}"; filename="{fname}"\r\n'
            "Content-Type: application/octet-stream\r\n\r\n".encode()
            + (DATA / "inputs" / fname).read_bytes()
            + b"\r\n"
        )
    return b"".join(parts) + f"--{BOUNDARY}--\r\n".encode()


def submit(record, tmp_path, monkeypatch):
    """POST a recorded submission; return the working and temporary
    directories, the status, the body and what the application logged."""
    workdir = tmp_path / "work"
    tmpdir = tmp_path / "tmp"
    workdir.mkdir()
    tmpdir.mkdir()
    monkeypatch.chdir(workdir)
    monkeypatch.setattr(tempfile, "tempdir", str(tmpdir))

    body = multipart(record["fields"], record["uploads"])
    errors = io.StringIO()
    environ = {}
    wsgiref.util.setup_testing_defaults(environ)
    environ.update(
        {
            "REQUEST_METHOD": "POST",
            "PATH_INFO": "/cladecompare",
            "CONTENT_TYPE": "multipart/form-data; boundary=" + BOUNDARY,
            "CONTENT_LENGTH": str(len(body)),
            "wsgi.input": io.BytesIO(body),
            "wsgi.errors": errors,
        }
    )
    status = []
    response = b"".join(
        bottle.default_app()(
            environ, lambda stat, headers, exc_info=None: status.append(stat)
        )
    )
    return workdir, tmpdir, status[0], response, errors.getvalue()


def load(case):
    text = (REFERENCE / case / "run.json").read_text(encoding="utf-8")
    return json.loads(text)


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(case, marks=[pytest.mark.integration] if programs else [])
        for case, programs in sorted(CASES.items())
    ],
)
def test_reference_submission(case, tmp_path, monkeypatch, require):
    require(*CASES[case])
    record = load(case)
    workdir, tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == record["status"], errors
    if record["error"]:
        error = record["error"].strip().splitlines()[-1].split(":")[0]
        assert ERRORS.get(case, error) in errors
    else:
        assert response == (REFERENCE / case / "response.html").read_bytes()
    assert os.listdir(workdir) == record["written_files"] == []
    if case not in SPOOLED:
        leftovers = os.listdir(tmpdir)
        assert len(leftovers) == record["leftover_temp_files"], leftovers


def test_every_recorded_submission_is_tested():
    recorded = {path.name for path in REFERENCE.iterdir()}
    assert recorded == set(CASES) | OWN_TESTS


def test_form_offers_every_strategy():
    """The form offers each strategy that -s accepts, with the default checked."""
    page = cladecomparelib.cladeweb.form()
    assert re.findall(r'name="strategy" value="(\w+)"', page) == list(core.STRATEGIES)
    assert re.findall(r'name="strategy" value="(\w+)" checked="checked"', page) == [
        core.DEFAULT_STRATEGY
    ]


@pytest.mark.parametrize("mode", ["pair", "single"])
def test_phospho_submission_completes(mode, tmp_path, monkeypatch):
    """phospho, which the reference records failing with TypeError, runs."""
    workdir, _tmpdir, status, response, errors = submit(
        load(f"web.{mode}.phospho"), tmp_path, monkeypatch
    )
    assert status == "200 OK", errors
    assert b"<table>" in response
    assert os.listdir(workdir) == []


def test_hypg_submission_matches_the_command_line(tmp_path, monkeypatch):
    """hypg, which the reference records failing with KeyError, runs.

    The form now offers every strategy that the command line accepts. The
    report equals the one that the recorded command-line case pair.hypg
    wrote from the same alignments in CMA form, apart from the title.
    """
    workdir, _tmpdir, status, response, errors = submit(
        load("web.pair.hypg"), tmp_path, monkeypatch
    )
    assert status == "200 OK", errors
    expected = expected_text("pair.hypg", "pair.hypg.html")
    assert expected.count(b"pair.hypg.out") == 2
    assert response == expected.replace(b"pair.hypg.out", b"web-pair-hypg")
    assert os.listdir(workdir) == []
