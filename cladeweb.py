#!/usr/bin/env python

"""Web browser interface for CladeCompare.

Input (GET):
    Form for submitting FG, [BG, HMM].
Output (POST):
    CladeReport (heat map) of submission.

"""
# ENH (report):
#   - asterisks up top link to PDF pairlogos
#   - below title, link to PyMOL script of PDB(s)

import logging
import os
import tempfile
import webbrowser

import bottle

from cladecomparelib import (core, pairlogo, pmlscript, urn, gtest, jsd,
                             phospho, report)


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
  <p>
    <label>
      <input type="radio" name="strategy" value="gtest" checked="checked" />
      G-test (goodness-of-fit)
    </label>
  </p>
  <p>
    <label>
      <input type="radio" name="strategy" value="urn" />
      Ball-in-urn model (binomial)
    </label>
  </p>
  <p>
    <label>
      <input type="radio" name="strategy" value="jsd" />
      Jensen-Shannon divergence
    </label>
  </p>
  <p>
    <label>
      <input type="radio" name="strategy" value="phospho" />
      Phosphorylation site conservation
    </label>
  </p>

  <p>
    Significance cutoff (alpha):
    <input type="text" name="alpha" value="0.05" />
  </p>


  <h2>Alignment profile</h2>
  <p>
    HMM (.hmm) profile: 
    <br />
    <input type="file" name="profile" size=50 />
  </p>

  <!--
  <h2>Structure</h2>
  <p>
    PDB ID:
    <input type="text" name="pdbid" />
    <br />
    or upload a
    PDB file:
    <br />
    <input type="file" name="pdbfile" size=50 />
  </p>
  -->

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
"""


# --- Routes ---

@bottle.get('/cladecompare')
def form():
    return FORM_HTML


# TODO - routes for downloading .pml, .pdf -- use static_file


@bottle.post('/cladecompare')
def form_submit():
    # ENH: pick a unique, informative name -- e.g. date or hostname
    name = bottle.request.forms.name

    seqfile1 = bottle.request.files.seqfile1
    if not hasattr(seqfile1, 'file'):
        return "Error: You need to specify at least one sequence file."
    seq1fname = handle2temp(seqfile1.file,
                            suffix=('.cma' if seqfile1.filename.endswith('.cma')
                                    else '.seq'))

    # Optional second sequence set -- if missing, do single mode
    seqfile2 = bottle.request.files.seqfile2
    if hasattr(seqfile2, 'file'):
        seq2fname = handle2temp(seqfile2.file,
                                suffix=('.cma' if
                                        seqfile2.filename.endswith('.cma') else
                                        '.seq'))
        if not name:
            name = "%s-vs-%s" % (seqfile1.filename.rsplit('.', 1)[0],
                                 seqfile2.filename.rsplit('.', 1)[0])
    else:
        seq2fname = ''
        if not name:
            name = seqfile1.filename

    # Optional HMM profile for alignment
    profile = bottle.request.files.profile
    # Optional HMM profile for alignment
    profile = bottle.request.files.profile
    if hasattr(profile, 'file'):
        if not profile.filename.endswith('.hmm'):
            return "HMM profile file name must end in .hmm"
        profname = handle2temp(profile.file, suffix='.hmm')
        logging.info("Aligning %s with profile %s", seq1fname, profname)
        fg_aln = core.hmm_align_and_read(profname, seq1fname)
        if seq2fname:
            logging.info("Aligning %s with profile %s", seq2fname, profname)
            bg_aln = core.hmm_align_and_read(profname, seq2fname)
    else:
        profname = ''
        fg_aln = core.read_aln(seq1fname, 'fasta')
        if seq2fname:
            bg_aln = core.read_aln(seq2fname, 'fasta')

    pdbfile = bottle.request.files.pdbfile
    if hasattr(pdbfile, 'file'):
        if not profname:
            return ("Error: to generate a PyMOL script for a PDB file you must"
                    "also specify an HMM profile")
        pdbfname = handle2temp(pdbfile.file)
        logging.info("Aligning %s with profile %s", pdbfile.filename, profname)
        pdb_rec, pdb_resnums, pdb_inserts = core.pdb_hmm(profname,
                                                         pdbfname)
        pdb_data = [(pdbfname, pdb_rec, pdb_resnums, pdb_inserts)]
    else:
        pdbfname = ''
        pdb_data = None

    # Mutually exclusive with pdbfname (above):
    pdbid = bottle.request.forms.pdbid
    if pdbid:
        # If PDB ID: .pml should use "fetch" instead of "load"?
        # Can get this info w/o dl'ing actual PDB file (e.g. via FASTA)?
        pass

    stat_module = dict(gtest=gtest, urn=urn, jsd=jsd, phospho=phospho,
                      )[bottle.request.forms.strategy]
    try:
        alpha = float(bottle.request.forms.alpha)
        if not 0.0 <= alpha <= 1.0:
            raise ValueError
    except ValueError:
        return "Error: alpha must be a number between 0 and 1"

    _fdo, tmp_output = tempfile.mkstemp(suffix='.out')
    os.close(_fdo)
    _fdp, tmp_pattern = tempfile.mkstemp(suffix='.pttrn')
    os.close(_fdp)

    # Run the algorithms...

    if seq2fname:
        # Pair mode
        fg_clean, bg_clean, hits = core.process_pair(fg_aln, bg_aln,
                                                     stat_module)
        core.process_output(fg_clean, bg_clean, hits, alpha,
                            tmp_output, tmp_pattern,
                            pdb_data)
    else:
        # Single mode
        aln, hits = core.process_one(fg_aln, stat_module)
        core.process_output(aln, None, hits, alpha,
                            tmp_output, tmp_pattern,
                            pdb_data)

    # Get the HTML report data
    contents = report.do_single(tmp_output, tmp_pattern)[1]

    cleanup(seq1fname)
    cleanup(seq2fname)
    cleanup(profname)
    cleanup(tmp_output)
    cleanup(tmp_pattern)

    return report.html_page_tpl % dict(title=name, contents=contents)


# --- Helpers ---

def handle2temp(handle, suffix=''):
    """Write file handle contents to a temporary file, return tempfile name."""
    _fd, fname = tempfile.mkstemp(suffix=suffix)
    os.write(_fd, handle.read())
    os.close(_fd)
    return fname


def cleanup(fname):
    """Remove a temporary file that may or may not exist."""
    if os.path.isfile(fname):
        try:
            os.remove(fname)
            print "Cleaned up", fname
        except OSError:
            print "Failed to clean up", fname


# --- Run ---

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format="%(module)s [@%(lineno)s]: %(message)s")
    webbrowser.open("http://localhost:8080/cladecompare")
    bottle.run(host='localhost', port=8080)
