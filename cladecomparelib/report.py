"""Generate an HTML report summarizing CladeCompare output."""

import html
import logging
import re
from math import isfinite, log10, nan

from . import InputError, reading

logger = logging.getLogger(__name__)

html_page_tpl = """\
<html>
  <head>
    <title>%(title)s</title>
    <style>
    * { font-family: Verdana, sans; }
    body { margin: 0; }
    h1 { font: 1.6em bold Verdana, sans;
         background-color: #eee;
         border-bottom: 1px solid black;
         padding: 0 0 .5em 0; margin: 0 0 1em 0;
    }
    tr { border: 0; padding: 0; margin: 0; }
    td { font-family: Courier, monospace; border: 0; padding: 1px; margin: 0; }
    td.label { font-family: Verdana, sans; padding: 1px 3px; }
    td.del { background-color: #bbb }
    td.ins { background-color: #fe8 }
    td.b5 { background-color: #99f; }
    td.b4 { background-color: #bbf }
    td.b3 { background-color: #ccf }
    td.b2 { background-color: #ddf }
    td.b1 { background-color: #eef }
    td.wt { background-color: white }
    td.r1 { background-color: #fee }
    td.r2 { background-color: #fcc }
    td.r3 { background-color: #faa }
    td.r4 { background-color: #f88; color: white; }
    td.r5 { background-color: #f66; color: white; }
    td.r6 { background-color: #f44; color: white; }
    td.r7 { background-color: #e22; color: white; }
    td.r8 { background-color: #d11; color: white; }
    td.r9 { background-color: #c00; color: white; }
    tr.pos td { font: .9em bold Verdana, sans; padding: 0; }
    </style>
  </head>
  <body>
    <h1>%(title)s</h1>
    %(contents)s
  </body>
</html>
"""


def html_page(title, contents):
    """Render a report page; the title is plain text, the contents are HTML."""
    return html_page_tpl % dict(title=html.escape(title), contents=contents)


# A line of a '.out' file as core.write_pvalues writes it, without the
# surrounding white space, e.g. "Y (Y) 1 : prob=0.00363161\t**".
OUT_LINE = re.compile(
    r"(?P<fg>\S) \((?P<bg>\S)\) (?P<posn>\d+) : prob=(?P<prob>[^\s*]+)"
    r"(?:\s+(?P<stars>\*+))?"
)

# A site of a pattern file: a residue followed by its column number, e.g. D34.
PATTERN_SITE = re.compile(r"\D*(?P<posn>\d+)")


def parse_out(fname):
    """Read the sites of a '.out' file.

    Return a list of (foreground residue, background residue, column number,
    value, number of asterisks), one for each line. A file without sites, a
    line of any other shape, or a value that is not a finite number raises
    InputError.
    """
    sites = []
    with reading(fname), open(fname) as infile:
        for lineno, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue
            match = OUT_LINE.fullmatch(line)
            if match is None:
                raise InputError(
                    f"{fname}, line {lineno}: not a line of cladecompare output: "
                    f"{line!r}"
                )
            try:
                prob = float(match["prob"])
            except ValueError:
                prob = nan
            if not isfinite(prob):
                raise InputError(
                    f"{fname}, line {lineno}: the value {match['prob']!r} "
                    "is not a finite number"
                )
            nstars = len(match["stars"] or "")
            sites.append((match["fg"], match["bg"], int(match["posn"]), prob, nstars))
    if not sites:
        raise InputError(f"{fname} holds no lines of cladecompare output")
    return sites


def parse_pattern(fname):
    """Read the column numbers of the first pattern in a '.pttrn' file, sorted.

    A first line without the ':' that follows the pattern's number, or a site
    without a column number, raises InputError.
    """
    with reading(fname), open(fname) as infile:
        # NB: only reading the first pattern, ignoring any others
        line = infile.readline()
    _idx, colon, site_ln = line.partition(":")
    if not colon:
        raise InputError(f"{fname} is not a pattern file: its first line has no ':'")
    site_ln = site_ln.strip()
    if not site_ln:
        logger.warning("Empty pattern file: %s", fname)
        return []
    posns = []
    for site in site_ln.split(","):
        match = PATTERN_SITE.fullmatch(site.strip())
        if match is None:
            raise InputError(f"{fname}: the site {site!r} has no column number")
        posns.append(int(match["posn"]))
    return sorted(posns)


def pattern_marks(patternfname, posns, outfname):
    """An asterisk for each column in posns that the pattern file names.

    Without a pattern file, no column is marked. A site that is not among the
    columns of the '.out' file outfname raises InputError.
    """
    marks = [""] * len(posns)
    if patternfname:
        index = {posn: i for i, posn in enumerate(posns)}
        for posn in parse_pattern(patternfname):
            if posn not in index:
                raise InputError(
                    f"{patternfname} names column {posn}, which {outfname} does "
                    f"not have; its columns are numbered {min(posns)} to "
                    f"{max(posns)}"
                )
            marks[index[posn]] = "*"
    return marks


def p2class(fg, bg, pval):
    """Assign a table cell (td) CSS class given a p-value.

    Heat map:
        higher p-values -> more blue
        p-values near significance -> white
        smaller p-values -> more red

    Foreground deletions -> gray, inserts -> gold.

    """
    if fg == "-":
        return "del"
    if bg == "-":
        return "ins"
    try:
        nlp = -log10(pval)
    except ValueError:
        # pval is 0, but let log10 decide that
        return "r9"
    if nlp < 0.00001:
        return "b5"
    if nlp < 0.0005:
        return "b4"
    if nlp <= 0.005:
        return "b3"
    if nlp <= 0.05:
        return "b2"
    if nlp <= 0.3:
        return "b1"
    if nlp <= 1.0:
        return "wt"
    if nlp >= 9.0:
        return "r9"
    # r1, r2, ..., r8
    return "r%d" % int((nlp + 1) / 2)


# Formatting


def table_single(pvalues, classes, pattern, foreground, background, positions):
    assert len(foreground) == len(background) == len(positions)
    out = "<table>"
    if any(pattern):
        out += tr_plain("", pattern)
    out += tr_class("foreground", foreground, pvalues, classes)
    out += tr_plain("background", background)
    out += tr_pos(positions)
    out += "</table>"
    return out


def table_multi(labels, rowcells, rowclasses, rowpvalues, pattern, positions):
    assert len(labels) == len(rowclasses) == len(rowcells) == len(rowpvalues)
    assert len(positions) == len(rowclasses[0]) == len(rowcells[0])
    out = "<table>"
    if any(pattern):
        out += tr_plain("", pattern)
    for label, cells, classes, pvals in zip(
        labels, rowcells, rowclasses, rowpvalues, strict=True
    ):
        out += tr_class(label, cells, pvals, classes)
    out += tr_pos(positions)
    out += "</table>"
    return out


def tr_plain(label, cells):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % label
    out += "".join(["<td>%s</td>" % c for c in cells])
    out += "\n</tr>"
    return out


def tr_class(label, cells, pvalues, classes):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % html.escape(label)
    out += "".join(
        [
            '<td class="%s" title="p=%s">%s</td>' % (cls, pval, txt)
            for cls, pval, txt in zip(classes, pvalues, cells, strict=True)
        ]
    )
    out += "\n</tr>"
    return out


def tr_pos(cells):
    out = "<tr class='pos'>\n<td></td>\n"
    out += "".join(["<td>%s</td>" % c for c in cells])
    out += "\n</tr>"
    return out


# Logic


def wrap_data(data, width=50):
    """list of atoms -> list of lists, chunked at width"""
    chunks = []
    curr, remain = None, data
    while len(remain) > width:
        curr, remain = remain[:width], remain[width:]
        chunks.append(curr)
    chunks.append(remain)
    return chunks


def do_single(fname, patternfname):
    """Build an HTML document for a single '.out' file and optional pattern."""
    fg_seq, bg_seq, posns, probs, _nstarses = zip(*parse_out(fname), strict=True)
    classes = [
        p2class(fg, bg, p) for fg, bg, p in zip(fg_seq, bg_seq, probs, strict=True)
    ]
    pattern = pattern_marks(patternfname, posns, fname)
    posns = [str(p) if not p % 10 else "" for p in posns]
    # wrap it up
    contents = ""
    chunked_data = map(wrap_data, [probs, classes, pattern, fg_seq, bg_seq, posns])
    for layer in zip(*chunked_data, strict=True):
        contents += table_single(*layer)
        contents += "<p />"
    return fname, contents


def do_multi(fnames, patternfname):
    """Build an HTML document combining multiple '.out' files.

    Each file becomes one row of the table, under one row of column numbers,
    so the files must number the same columns.
    """
    rowcells = []  # dim = #fnames x #cols
    rowclasses = []
    rowpvalues = []
    posns = None
    for fname in fnames:
        fg_seq, bg_seq, file_posns, probs, _nstarses = zip(
            *parse_out(fname), strict=True
        )
        if posns is None:
            posns = file_posns
        elif file_posns != posns:
            raise InputError(
                "%s and %s do not number the same columns (%d and %d columns), "
                "so they cannot be shown in one table"
                % (fnames[0], fname, len(posns), len(file_posns))
            )
        rowcells.append(fg_seq)
        rowclasses.append(
            [
                p2class(fg, bg, p)
                for fg, bg, p in zip(fg_seq, bg_seq, probs, strict=True)
            ]
        )
        rowpvalues.append(probs)
    pattern = pattern_marks(patternfname, posns, fnames[0])
    posns = [str(p) if not p % 10 else "" for p in posns]
    # wrap it up
    contents = ""
    labels = [fn[: -len(".out")] if fn.endswith(".out") else fn for fn in fnames]
    chunked_posns = wrap_data(posns)
    chunked_pattern = wrap_data(pattern)
    wrapped_row_cells = zip(*map(wrap_data, rowcells), strict=True)
    wrapped_row_classes = zip(*map(wrap_data, rowclasses), strict=True)
    wrapped_row_pvalues = zip(*map(wrap_data, rowpvalues), strict=True)
    for w_rowcells, w_rowclasses, w_rowpvalues, ptnrow, positions in zip(
        wrapped_row_cells,
        wrapped_row_classes,
        wrapped_row_pvalues,
        chunked_pattern,
        chunked_posns,
        strict=True,
    ):
        contents += table_multi(
            labels, w_rowcells, w_rowclasses, w_rowpvalues, ptnrow, positions
        )
        contents += "<p />"
    return ", ".join(fnames), contents
