"""Regression tests that replay the recorded cladeweb form submissions.

Each case submits the form fields and uploads recorded in
``tests/data/py27/web/<case>/run.json`` to the cladeweb WSGI application,
with the working and temporary directories redirected to scratch
directories. The response status must be the recorded one; a successful
response must equal the recorded body, or, for the submissions in
``REPORTS``, the current report of the command-line case whose report the
reference records it equal to; and a failed one must have logged the
recorded exception. A recorded failure in ``REPAIRED`` now succeeds, and its
response must equal the current report of the command-line case given
there. No case may write into the working directory or leave a file of its
own in the temporary directory, which the reference records it doing.

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
from test_core import aln

import cladecompare.cladeweb  # registers the routes
from cladecompare import core, report

DATA = Path(__file__).resolve().parent / "data"
REFERENCE = DATA / "py27" / "web"

# Recorded submissions, each compared with its recorded response or, for those
# in REPORTS, with the current report of a command-line case, and the programs
# each needs.
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

# Under bottle 0.13, the submissions that upload CDK.hmm leave the named
# temporary files to which bottle spools a request body or upload larger than
# 100 KiB; they are bottle's own, outside the request's directory.
SPOOLED = {"web.hmm.pair.gtest", "web.pdb.1HCK.gtest"}

# Recorded failures that now succeed, and the command-line case whose current
# report, apart from the title, each response must equal. The application saved
# an uploaded structure under the fixed name pdbfile, from which pmlscript
# could not read a PDB ID; the PDB ID now comes from the structure's HEADER
# record, and the upload keeps the name it was uploaded with.
REPAIRED = {"web.pdb.1HCK.gtest": "pdb.1HCK.gtest"}

# Recorded submissions whose response the reference records equal, apart from
# the title, to the report of a recorded command-line case on the same
# alignments. The G-test and urn have since changed, so each response must
# equal the current report of its case. web.single.gtest submits alpha 0.005
# where single.gtest.a.001 used 0.001, and its report is the same because the
# Benjamini-Hochberg procedure accepts more than 50 columns at both (239 at
# 0.005, 237 at 0.001), and the selection keeps the same 50 most significant.
REPORTS = {
    "web.hmm.pair.gtest": "hmm.pair.gtest",
    "web.pair-cma.gtest": "pair.gtest",
    "web.pair.gtest": "pair.gtest",
    "web.pair.urn": "pair.urn",
    "web.single.gtest": "single.gtest.a.001",
    "web.single.urn": "single.urn",
}

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
    assert status == ("200 OK" if case in REPAIRED else record["status"]), errors
    title = dict(record["fields"])["name"]
    if case in REPAIRED:
        assert response == cli_report(REPAIRED[case], title)
    elif record["error"]:
        error = record["error"].strip().splitlines()[-1].split(":")[0]
        assert error in errors
    elif case in REPORTS:
        assert response == cli_report(REPORTS[case], title)
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


def test_recorded_responses_are_the_recorded_reports():
    """Each submission in REPORTS recorded the report of its command-line case,
    apart from the title."""
    for case, cli_case in sorted(REPORTS.items()):
        title = dict(load(case)["fields"])["name"]
        recorded = (REFERENCE / case / "response.html").read_bytes()
        assert recorded == cli_report(cli_case, title, reference_text), case


def test_form_offers_every_strategy():
    """The form offers each strategy that -s accepts, with the default checked."""
    page = cladecompare.cladeweb.form()
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


def test_single_alignment_hypg_submission_is_rejected(tmp_path, monkeypatch):
    """hypg on one alignment returns an error message and leaves no files."""
    record = load("web.single.phospho")
    record["fields"][1] = ["strategy", "hypg"]
    workdir, tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert response == (
        b"Error: hypg compares two alignments; give a foreground and a background"
    )
    assert os.listdir(workdir) == []
    assert os.listdir(tmpdir) == []


def cli_report(case, title, source=expected_text):
    """The report of a recorded command-line case, whose title names its .out
    file, retitled as the report of a submission."""
    expected = source(case, case + ".html")
    outfname = (case + ".out").encode()
    assert expected.count(outfname) == 2
    return expected.replace(outfname, title.encode())


def test_hypg_submission_matches_the_command_line(tmp_path, monkeypatch):
    """hypg, which the reference records failing with KeyError, runs.

    The form now offers every strategy that the command line accepts. The
    report equals the current report of the command-line case pair.hypg, on
    the same alignments in CMA form, apart from the title.
    """
    workdir, _tmpdir, status, response, errors = submit(
        load("web.pair.hypg"), tmp_path, monkeypatch
    )
    assert status == "200 OK", errors
    assert response == cli_report("pair.hypg", "web-pair-hypg")
    assert os.listdir(workdir) == []


def test_weighted_submission_matches_the_command_line(tmp_path, monkeypatch):
    """The form's weight checkbox gives the report of cladecompare -w."""
    assert '<input type="checkbox" name="weight" />' in cladecompare.cladeweb.form()
    record = load("web.pair.gtest")
    record["fields"].append(["weight", "on"])
    _workdir, _tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert response == cli_report("pair.gtest.w", "web-pair-gtest")


def test_profile_submission_keeps_unoccupied_match_states(tmp_path, monkeypatch):
    """With a profile, a match state that no sequence occupies keeps its
    column, as on the command line."""
    aligned = [
        aln(("a", "A-WWK"), ("b", "C-WWK")),
        aln(("x", "E-GGK"), ("y", "F-GGK")),
    ]
    monkeypatch.setattr(core, "hmm_align_and_read", lambda *_args: aligned.pop(0))
    monkeypatch.setattr(
        report, "do_single", lambda out, _pattern: (None, Path(out).read_text())
    )
    _workdir, _tmpdir, status, response, errors = submit(
        load("web.hmm.pair.gtest"), tmp_path, monkeypatch
    )
    assert status == "200 OK", errors
    assert b"- (-) 2 : prob=1\t\n" in response
    assert b" 5 : prob=" in response


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
    caplog.set_level(logging.INFO, logger="cladecompare.pairlogo")
    _workdir, tmpdir, status, _response, errors = submit(
        load("web.pair.gtest"), tmp_path, monkeypatch
    )
    assert status.startswith("500"), errors
    assert "report failed" in errors
    # The pair-logo PDFs were written before the failure.
    assert any(message.endswith(".pdf") for message in caplog.messages)
    assert os.listdir(tmpdir) == []


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("alpha", "1.5", b"Error: alpha must be a number between 0 and 1, not '1.5'"),
        (
            "strategy",
            "nosuch",
            b"Error: Unknown strategy 'nosuch'; choose one of "
            + ", ".join(core.STRATEGIES).encode(),
        ),
        (
            "seqfile1",
            "scttlre-domain.fasta",
            b"Error: Cannot read scttlre-domain.fasta as a fasta alignment: "
            b"Sequences must all be the same length",
        ),
    ],
)
def test_a_rejected_submission_leaves_no_temporary_files(
    tmp_path, monkeypatch, field, value, message
):
    """A submission rejected with an error message removes its uploads, and
    the message names an upload by its file name."""
    record = load("web.pair.gtest")
    for entries in record["fields"], record["uploads"]:
        for entry in entries:
            if entry[0] == field:
                entry[1] = value
    _workdir, tmpdir, status, response, errors = submit(record, tmp_path, monkeypatch)
    assert status == "200 OK", errors
    assert response == message
    assert os.listdir(tmpdir) == []


def test_browser_opens_once_the_server_listens(monkeypatch):
    """The browser opens the form's address only when it can connect."""
    opened = []

    def open_url(url):
        port = int(url.split(":")[2].split("/")[0])
        with socket.create_connection(("localhost", port), timeout=5):
            opened.append(url)

    monkeypatch.setattr(webbrowser, "open", open_url)
    server = cladecompare.cladeweb.BrowserServer(("localhost", 0), WSGIRequestHandler)
    server.server_close()
    assert opened == ["http://localhost:%d/cladecompare" % server.server_port]


def test_browser_does_not_open_when_the_port_is_taken(monkeypatch):
    """A server already listening on the same address and port stops cladeweb
    before the browser opens."""
    opened = []
    monkeypatch.setattr(webbrowser, "open", opened.append)
    with socket.create_server(("localhost", 0)) as taken:
        with pytest.raises(OSError):
            cladecompare.cladeweb.BrowserServer(
                taken.getsockname()[:2], WSGIRequestHandler
            )
    assert opened == []
