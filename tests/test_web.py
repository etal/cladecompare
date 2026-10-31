"""Regression tests that replay the recorded cladeweb form submissions.

Each case submits the form fields and uploads recorded in
``tests/data/py27/web/<case>/run.json`` to the cladeweb WSGI application,
with the working and temporary directories redirected to scratch
directories. The response status must be the recorded one; a successful
response must equal the recorded body, and a failed one must have logged the
recorded exception, or the one that ``ERRORS`` names. No case may write into
the working directory or leave a file of its own in the temporary directory,
which the reference records it doing.

The cases that run hmmalign carry the ``integration`` mark; see
``test_cli_cases.py``.
"""

import io
import json
import logging
import os
import re
import socket
import tempfile
import webbrowser
import wsgiref.util
from pathlib import Path
from wsgiref.simple_server import WSGIRequestHandler

import bottle
import pytest
from test_cli_cases import expected_text, reference_text

import cladecomparelib.cladeweb  # registers the routes
from cladecomparelib import core, report

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
    "web.single.jsd": [],
    "web.single.urn": [],
}

# Under bottle 0.13, the submissions that upload CDK.hmm leave the named
# temporary files to which bottle spools a request body or upload larger than
# 100 KiB; they are bottle's own, outside the request's directory.
SPOOLED = {"web.hmm.pair.gtest", "web.pdb.1HCK.gtest"}

# Recorded failures that now raise a different exception by design, and the
# exception each raises: pmlscript reports a file name it cannot read with
# InputError instead of an assertion.
ERRORS = {"web.pdb.1HCK.gtest": "InputError"}

# Recorded submissions with a test of their own below.
OWN_TESTS = {
    "web.pair.hypg",
    "web.pair.phospho",
    "web.single.gtest",
    "web.single.phospho",
}

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
    leftovers = os.listdir(tmpdir)
    if case in SPOOLED:
        leftovers = [name for name in leftovers if name.startswith("cladeweb-")]
    assert leftovers == []


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


def assert_matches_cli_report(response, case, title):
    """Assert that a response equals the report of a recorded command-line
    case, whose title names its .out file, apart from the title."""
    expected = expected_text(case, case + ".html")
    outfname = (case + ".out").encode()
    assert expected.count(outfname) == 2
    assert response == expected.replace(outfname, title.encode())


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
    assert_matches_cli_report(response, "pair.hypg", "web-pair-hypg")
    assert os.listdir(workdir) == []


def test_weighted_submission_matches_the_command_line(tmp_path, monkeypatch):
    """The form's weight checkbox gives the report of cladecompare -w."""
    assert '<input type="checkbox" name="weight" />' in cladecomparelib.cladeweb.form()
    record = load("web.pair.gtest")
    record["fields"].append(["weight", "on"])
    _workdir, _tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert_matches_cli_report(response, "pair.gtest.w", "web-pair-gtest")


def test_single_gtest_submission_matches_the_command_line(tmp_path, monkeypatch):
    """Single-alignment gtest gives the report of the command line.

    The form submits the alignment of single.gtest.a.001 in FASTA form with
    alpha 0.005, and the reference records the report of that case, apart
    from the title, because the Benjamini-Hochberg search stops at the same
    50 sites at both alphas. The p-values have since moved with the expected
    counts of single-mode gtest, and the response follows the case.
    """
    recorded = (REFERENCE / "web.single.gtest" / "response.html").read_bytes()
    report_html = reference_text("single.gtest.a.001", "single.gtest.a.001.html")
    assert recorded == report_html.replace(
        b"single.gtest.a.001.out", b"web-single-gtest"
    )
    workdir, tmpdir, status, response, errors = submit(
        load("web.single.gtest"), tmp_path, monkeypatch
    )
    assert status == "200 OK", errors
    assert_matches_cli_report(response, "single.gtest.a.001", "web-single-gtest")
    assert os.listdir(workdir) == []
    assert os.listdir(tmpdir) == []


def test_submission_name_is_escaped(tmp_path, monkeypatch):
    """The submission name, which titles the report, is text, not HTML."""
    record = load("web.pair.gtest")
    record["fields"][0] = ["name", "<b>R&D</b>"]
    _workdir, _tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert b"<b>" not in response
    assert response.count(b"&lt;b&gt;R&amp;D&lt;/b&gt;") == 2


def test_a_failing_submission_leaves_no_temporary_files(tmp_path, monkeypatch, caplog):
    """An exception after every output file is written removes them all."""

    def fail(*_args):
        raise RuntimeError("report failed")

    monkeypatch.setattr(report, "do_single", fail)
    caplog.set_level(logging.INFO, logger="cladecomparelib.pairlogo")
    _workdir, tmpdir, status, _response, errors = submit(
        load("web.pair.gtest"), tmp_path, monkeypatch
    )
    assert status.startswith("500"), errors
    assert "report failed" in errors
    # The pair-logo PDFs were written before the failure.
    assert any(message.endswith(".pdf") for message in caplog.messages)
    assert os.listdir(tmpdir) == []


def test_a_rejected_submission_leaves_no_temporary_files(tmp_path, monkeypatch):
    """A submission rejected with an error message removes its uploads."""
    record = load("web.pair.gtest")
    record["fields"][2] = ["alpha", "1.5"]
    _workdir, tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert response == b"Error: alpha must be a number between 0 and 1"
    assert os.listdir(tmpdir) == []


def test_browser_opens_once_the_server_listens(monkeypatch):
    """The browser opens the form's address only when it can connect."""
    opened = []

    def open_url(url):
        port = int(url.split(":")[2].split("/")[0])
        with socket.create_connection(("localhost", port), timeout=5):
            opened.append(url)

    monkeypatch.setattr(webbrowser, "open", open_url)
    server = cladecomparelib.cladeweb.BrowserServer(
        ("localhost", 0), WSGIRequestHandler
    )
    server.server_close()
    assert opened == ["http://localhost:%d/cladecompare" % server.server_port]


def test_browser_does_not_open_when_the_port_is_taken(monkeypatch):
    """A server already listening on the same address and port stops cladeweb
    before the browser opens."""
    opened = []
    monkeypatch.setattr(webbrowser, "open", opened.append)
    with socket.create_server(("localhost", 0)) as taken:
        with pytest.raises(OSError):
            cladecomparelib.cladeweb.BrowserServer(
                taken.getsockname()[:2], WSGIRequestHandler
            )
    assert opened == []
