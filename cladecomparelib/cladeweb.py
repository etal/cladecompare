"""Web browser interface for CladeCompare.

Input (GET):
    Form for submitting FG, [BG, HMM].
Output (POST):
    CladeReport (heat map) of submission.

"""

# ENH (in the report):
#   - asterisks up top link to PDF pairlogos
#   - below title, link to PyMOL script of PDB(s)

import argparse
import html
import logging
import os
import tempfile
import webbrowser
from wsgiref.simple_server import WSGIServer

import bottle

from cladecomparelib import InputError, configure_logging, core, report

logger = logging.getLogger(__name__)

STRATEGY_RADIO = """\
  <p>
    <label>
      <input type="radio" name="strategy" value="%(name)s"%(checked)s />
      %(label)s
    </label>
  </p>
"""

# ENH: include textarea as alternative to each file upload input
# ENH: dropdown alpha (or textinput w/ JS validation)
FORM_HTML = """\
<html>
<body>
<h1>CladeCompare</h1>

<form action="cladecompare" method="post" enctype="multipart/form-data">
  <p>
    Submission name:
    <input type="text" name="name" />
  </p>

  <h2>Sequences</h2>
  <p>
    Sequence file 1 (required):
    <br />
    <input type="file" name="seqfile1" size=50 />
  </p>

  <p>
    Sequence file 2:
    <br />
    <input type="file" name="seqfile2" size=50 />
  </p>

  <h2>Statistical strategy</h2>
%(strategies)s
  <p>
    Significance cutoff (alpha):
    <input type="text" name="alpha" value="0.005" />
  </p>

  <p>
    <label>
      <input type="checkbox" name="weight" />
      Weight the sequences for the statistical test (cladecompare -w)
    </label>
  </p>


  <h2>Alignment profile</h2>
  <p>
    HMM (.hmm) profile:
    <br />
    <input type="file" name="profile" size=50 />
  </p>

  <p />

  <p><input type="submit" /></p>
</form>

<hr />

<p>Project page: <a
href="https://github.com/etal/cladecompare">https://github.com/etal/cladecompare</a></p>

<p>If you use this software in a publication, please cite our paper that
describes it:</p>

<blockquote>Talevich, E. & Kannan, N. (2013)
<a href="https://doi.org/10.1186/1471-2148-13-117">Structural and
evolutionary adaptation of rhoptry kinases and pseudokinases, a family of
coccidian virulence factors</a>.
<i>BMC Evolutionary Biology</i> 13:117 doi:10.1186/1471-2148-13-117
</blockquote>

</body>
</html>
""" % {
    "strategies": "".join(
        STRATEGY_RADIO
        % {
            "name": name,
            "checked": ' checked="checked"' if name == core.DEFAULT_STRATEGY else "",
            "label": description[0].upper() + description[1:],
        }
        for name, (description, _module) in core.STRATEGIES.items()
    )
}


# --- Routes ---


@bottle.get("/cladecompare")
def form():
    return FORM_HTML


# TODO - routes for downloading .pml, .pdf -- use static_file. The files of a
# request are removed when it ends, so they need a place that outlives it.


@bottle.post("/cladecompare")
def form_submit():
    # The uploads and every file the run writes, including the pair-logo PDFs
    # that process_output writes beside the pattern file, go into one
    # directory, which is removed when the request ends, however it ends.
    with tempfile.TemporaryDirectory(prefix="cladeweb-") as tmpdir:
        uploads = {}
        try:
            return compare_submission(tmpdir, uploads)
        except InputError as exc:
            # Name each upload by its file name and every other file in
            # tmpdir by its own name, not by a temporary path.
            message = str(exc)
            for fname, upload_name in uploads.items():
                message = message.replace(fname, upload_name)
            message = message.replace(tmpdir + os.sep, "")
            logger.warning("Rejected the submission: %s", message)
            return "Error: " + html.escape(message, quote=False)


def compare_submission(tmpdir, uploads):
    """Run the submitted form, writing files in tmpdir; return the response.

    Each upload is saved in tmpdir, and uploads maps the name it is saved
    under to its file name, as bottle normalizes it. Unusable input raises
    InputError.
    """
    # ENH: pick a unique, informative name -- e.g. date or hostname
    name = bottle.request.forms.name

    strategy = bottle.request.forms.strategy
    if strategy not in core.STRATEGIES:
        raise InputError(
            "Unknown strategy %r; choose one of %s"
            % (strategy, ", ".join(core.STRATEGIES))
        )
    _description, stat_module = core.STRATEGIES[strategy]
    try:
        alpha = core.parse_alpha(bottle.request.forms.alpha)
    except ValueError:
        raise InputError(
            "alpha must be a number between 0 and 1, not %r"
            % bottle.request.forms.alpha
        ) from None
    do_weight = bool(bottle.request.forms.weight)

    seqfile1 = bottle.request.files.seqfile1
    if not hasattr(seqfile1, "file"):
        raise InputError("You need to specify at least one sequence file.")
    suffix1 = ".cma" if seqfile1.filename.endswith(".cma") else ".seq"
    seq1fname = save_upload(seqfile1, tmpdir, "seqfile1" + suffix1, uploads)

    # Optional second sequence set -- if missing, do single mode
    seqfile2 = bottle.request.files.seqfile2
    if hasattr(seqfile2, "file"):
        suffix2 = ".cma" if seqfile2.filename.endswith(".cma") else ".seq"
        seq2fname = save_upload(seqfile2, tmpdir, "seqfile2" + suffix2, uploads)
        if not name:
            name = "%s-vs-%s" % (
                seqfile1.filename.rsplit(".", 1)[0],
                seqfile2.filename.rsplit(".", 1)[0],
            )
    else:
        seq2fname = ""
        if not name:
            name = seqfile1.filename

    # Optional HMM profile for alignment
    profile = bottle.request.files.profile
    if hasattr(profile, "file"):
        if not profile.filename.endswith(".hmm"):
            raise InputError("HMM profile file name must end in .hmm")
        profname = save_upload(profile, tmpdir, "profile.hmm", uploads)
        logger.info("Aligning %s with profile %s", seq1fname, profname)
        fg_aln = core.hmm_align_and_read(profname, seq1fname)
        if seq2fname:
            logger.info("Aligning %s with profile %s", seq2fname, profname)
            bg_aln = core.hmm_align_and_read(profname, seq2fname)
    else:
        profname = ""
        fg_aln = core.read_aln(seq1fname, "fasta")
        if seq2fname:
            bg_aln = core.read_aln(seq2fname, "fasta")

    # The form does not offer this field, and the PyMOL script that a
    # structure gives is written in tmpdir but not served.
    pdbfile = bottle.request.files.pdbfile
    if hasattr(pdbfile, "file"):
        if not profname:
            raise InputError(
                "to generate a PyMOL script for a PDB file you must "
                "also specify an HMM profile"
            )
        # The structure keeps its uploaded file name, in a directory of its
        # own, so that the script loads a file of that name and names a
        # structure without a PDB ID after it.
        os.mkdir(os.path.join(tmpdir, "structure"))
        pdbfname = save_upload(
            pdbfile, tmpdir, os.path.join("structure", pdbfile.filename), uploads
        )
        logger.info("Aligning %s with profile %s", pdbfile.filename, profname)
        pdb_rec, pdb_resnums, pdb_inserts = core.pdb_hmm(profname, pdbfname)
        pdb_data = [(pdbfname, pdb_rec, pdb_resnums, pdb_inserts)]
    else:
        pdb_data = []

    tmp_output = os.path.join(tmpdir, "cladecompare.out")
    tmp_pattern = os.path.join(tmpdir, "cladecompare.pttrn")

    # Run the algorithms...

    if seq2fname:
        # Pair mode
        fg_clean, bg_clean, hits = core.process_pair(
            fg_aln, bg_aln, stat_module, do_weight, bool(profname)
        )
        core.process_output(
            fg_clean, bg_clean, hits, alpha, tmp_output, tmp_pattern, pdb_data
        )
    else:
        # Single mode; a strategy that compares two alignments refuses it
        aln, hits = core.process_one(fg_aln, stat_module, do_weight)
        core.process_output(aln, None, hits, alpha, tmp_output, tmp_pattern, pdb_data)

    # Get the HTML report data
    contents = report.do_single(tmp_output, tmp_pattern)[1]
    return report.html_page(name, contents)


# --- Helpers ---


def save_upload(upload, dirname, basename, uploads):
    """Save an uploaded file as dirname/basename; return that file name.

    uploads maps the saved file's name to the upload's file name.
    """
    fname = os.path.join(dirname, basename)
    upload.save(fname)
    uploads[fname] = upload.filename
    return fname


class FormServer(WSGIServer):
    """A WSGI server that does not bind an address and port another one listens on.

    On Windows, SO_REUSEADDR lets a socket bind the address and port that
    another socket is listening on, so the option is set only on other
    systems, where it lets a new server reuse a port left in TIME_WAIT.
    """

    allow_reuse_address = os.name != "nt"


class BrowserServer(FormServer):
    """A form server that opens the form in a web browser once it listens.

    If the port cannot be bound, OSError is raised before this point, and
    no browser opens.
    """

    def server_activate(self):
        super().server_activate()
        webbrowser.open("http://localhost:%d/cladecompare" % self.server_port)


# --- Run ---


def main():
    AP = argparse.ArgumentParser(
        description="Serve the CladeCompare form on this computer."
    )
    AP.add_argument(
        "--port",
        type=int,
        default=8080,
        help="Port on localhost to serve the form from. (Default: 8080)",
    )
    AP.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open the form in a web browser.",
    )
    args = AP.parse_args()

    configure_logging()
    bottle.run(
        host="localhost",
        port=args.port,
        server_class=FormServer if args.no_browser else BrowserServer,
    )
