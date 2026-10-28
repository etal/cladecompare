"""Web browser interface for CladeCompare.

Input (GET):
    Form for submitting FG, [BG, HMM].
Output (POST):
    CladeReport (heat map) of submission.

"""

# TODO:
#   - take boolean do_weight as a form option (checkbox?)
# ENH (in the report):
#   - asterisks up top link to PDF pairlogos
#   - below title, link to PyMOL script of PDB(s)

import logging
import os
import tempfile
import webbrowser

import bottle

from cladecomparelib import configure_logging, core, report

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
href="http://github.com/etal/cladecompare">http://github.com/etal/cladecompare</a></p>

<p>If you use this software in a publication, please cite our paper that
describes it:</p>

<blockquote>Talevich, E. & Kannan, N. (2013)
<a href="http://www.biomedcentral.com/1471-2148/13/117">Structural and
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
        return compare_submission(tmpdir)


def compare_submission(tmpdir):
    """Run the submitted form, writing files in tmpdir; return the response."""
    # ENH: pick a unique, informative name -- e.g. date or hostname
    name = bottle.request.forms.name

    seqfile1 = bottle.request.files.seqfile1
    if not hasattr(seqfile1, "file"):
        return "Error: You need to specify at least one sequence file."
    suffix1 = ".cma" if seqfile1.filename.endswith(".cma") else ".seq"
    seq1fname = save_upload(seqfile1, tmpdir, "seqfile1" + suffix1)

    # Optional second sequence set -- if missing, do single mode
    seqfile2 = bottle.request.files.seqfile2
    if hasattr(seqfile2, "file"):
        suffix2 = ".cma" if seqfile2.filename.endswith(".cma") else ".seq"
        seq2fname = save_upload(seqfile2, tmpdir, "seqfile2" + suffix2)
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
            return "HMM profile file name must end in .hmm"
        profname = save_upload(profile, tmpdir, "profile.hmm")
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

    # The form does not offer this field. The structure is saved under a
    # fixed name, from which pmlscript cannot read a PDB ID, so a request that
    # sends one with a profile fails with InputError.
    pdbfile = bottle.request.files.pdbfile
    if hasattr(pdbfile, "file"):
        if not profname:
            return (
                "Error: to generate a PyMOL script for a PDB file you must "
                "also specify an HMM profile"
            )
        pdbfname = save_upload(pdbfile, tmpdir, "pdbfile")
        logger.info("Aligning %s with profile %s", pdbfile.filename, profname)
        pdb_rec, pdb_resnums, pdb_inserts = core.pdb_hmm(profname, pdbfname)
        pdb_data = [(pdbfname, pdb_rec, pdb_resnums, pdb_inserts)]
    else:
        pdb_data = None

    _description, stat_module = core.STRATEGIES[bottle.request.forms.strategy]
    try:
        alpha = float(bottle.request.forms.alpha)
        if not 0.0 <= alpha <= 1.0:
            raise ValueError
    except ValueError:
        return "Error: alpha must be a number between 0 and 1"

    tmp_output = os.path.join(tmpdir, "cladecompare.out")
    tmp_pattern = os.path.join(tmpdir, "cladecompare.pttrn")

    # Run the algorithms...

    if seq2fname:
        # Pair mode
        fg_clean, bg_clean, hits = core.process_pair(fg_aln, bg_aln, stat_module, False)
        core.process_output(
            fg_clean, bg_clean, hits, alpha, tmp_output, tmp_pattern, pdb_data
        )
    else:
        # Single mode
        aln, hits = core.process_one(fg_aln, stat_module, False)
        core.process_output(aln, None, hits, alpha, tmp_output, tmp_pattern, pdb_data)

    # Get the HTML report data
    contents = report.do_single(tmp_output, tmp_pattern)[1]
    return report.html_page_tpl % dict(title=name, contents=contents)


# --- Helpers ---


def save_upload(upload, dirname, basename):
    """Save an uploaded file as dirname/basename; return that file name."""
    fname = os.path.join(dirname, basename)
    upload.save(fname)
    return fname


# --- Run ---


def main():
    configure_logging()
    webbrowser.open("http://localhost:8080/cladecompare")
    bottle.run(host="localhost", port=8080)
