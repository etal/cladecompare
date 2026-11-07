"""Generate an HTML report summarizing CladeCompare output."""

import html
import logging
import re
from bisect import bisect_right
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


# A line of a '.out' file as core.write_values writes it, without the
# surrounding white space, e.g. "Y (Y) 1 : prob=0.00363161\t**". The label
# names the kind of value, a p-value ("prob") or a score ("score"). The
# asterisks follow from the value and alpha, so they are matched but not read.
OUT_LINE = re.compile(
    r"(?P<fg>\S) \((?P<bg>\S)\) (?P<posn>\d+) : "
    r"(?P<label>prob|score)=(?P<value>[^\s*]+)"
    r"(?:\s+\*+)?"
)

# A site of a pattern file: a residue followed by its column number, e.g. D34.
PATTERN_SITE = re.compile(r"\D*(?P<posn>\d+)")


def parse_out(fname):
    """Read the sites of a '.out' file.

    Return a list of (foreground residue, background residue, column number,
    label, value), one for each line, where label names the kind of value,
    "prob" or "score". A file without sites, a line of any other shape, or a
    value that is not a finite number raises InputError.
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
                value = float(match["value"])
            except ValueError:
                value = nan
            if not isfinite(value):
                raise InputError(
                    f"{fname}, line {lineno}: the value {match['value']!r} "
                    "is not a finite number"
                )
            sites.append(
                (match["fg"], match["bg"], int(match["posn"]), match["label"], value)
            )
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


# The upper bounds of the score classes below, in bits of Jensen-Shannon
# divergence, and the class of each interval.
SCORE_CUTS = (0.1, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
SCORE_CLASSES = ("b1", "wt", "r1", "r2", "r3", "r4", "r5", "r6")


def score2class(fg, bg, score):
    """Assign a table cell (td) CSS class given a score.

    Heat map, from the score of identical compositions, 0, to that of
    disjoint ones, 1:
        below 0.1 -> blue
        0.1 to 0.4 -> white
        0.4 and above -> red, one step darker per 0.1

    Under random division of the example pair's sequences into a foreground
    and a background of the same sizes, nearly every column scores below 0.4.
    A score is not a p-value, so a red cell is not significant.

    Foreground deletions -> gray, inserts -> gold.
    """
    if fg == "-":
        return "del"
    if bg == "-":
        return "ins"
    return SCORE_CLASSES[bisect_right(SCORE_CUTS, score)]


# For each label of a '.out' line, the name of the value in a cell's title
# and the function that colors the cell.
VALUE_SCALES = {"prob": ("p", p2class), "score": ("score", score2class)}


# Formatting


def tables(pattern, posns, rows, width=50):
    """Draw rows of cells as HTML tables of at most width columns each.

    pattern holds the mark above each column, and posns the column numbers,
    of which every tenth is shown below the rows. Each row is (label, cells,
    titles, classes); a row whose classes are None is drawn plain.
    """
    shown = [str(p) if not p % 10 else "" for p in posns]
    out = ""
    for start in range(0, len(posns), width):
        cols = slice(start, start + width)
        out += "<table>"
        if any(pattern[cols]):
            out += tr_plain("", pattern[cols])
        for label, cells, titles, classes in rows:
            if classes is None:
                out += tr_plain(label, cells[cols])
            else:
                out += tr_class(label, cells[cols], titles[cols], classes[cols])
        out += tr_pos(shown[cols])
        out += "</table><p />"
    return out


def tr_plain(label, cells):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % html.escape(label)
    out += "".join(["<td>%s</td>" % c for c in cells])
    out += "\n</tr>"
    return out


def tr_class(label, cells, titles, classes):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % html.escape(label)
    out += "".join(
        [
            '<td class="%s" title="%s">%s</td>' % (cls, title, txt)
            for cls, title, txt in zip(classes, titles, cells, strict=True)
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


def read_columns(fname):
    """Read the columns of a '.out' file with the cell title and class of each.

    Return the foreground and background consensus, the column numbers, the
    titles, which name and give each value, and the classes, each as a
    sequence over the columns.
    """
    fg_seq, bg_seq, posns, titles, classes = [], [], [], [], []
    for fg, bg, posn, label, value in parse_out(fname):
        name, value2class = VALUE_SCALES[label]
        fg_seq.append(fg)
        bg_seq.append(bg)
        posns.append(posn)
        titles.append("%s=%s" % (name, value))
        classes.append(value2class(fg, bg, value))
    return fg_seq, bg_seq, posns, titles, classes


def do_single(fname, patternfname):
    """Build an HTML document for a single '.out' file and optional pattern."""
    fg_seq, bg_seq, posns, titles, classes = read_columns(fname)
    rows = [
        ("foreground", fg_seq, titles, classes),
        ("background", bg_seq, None, None),
    ]
    pattern = pattern_marks(patternfname, posns, fname)
    return fname, tables(pattern, posns, rows)


def do_multi(fnames, patternfname):
    """Build an HTML document combining multiple '.out' files.

    Each file becomes one row of the table, under one row of column numbers,
    so the files must number the same columns. Each cell is colored on the
    scale that the label of its line names, a p-value's or a score's.
    """
    rows = []
    posns = None
    for fname in fnames:
        fg_seq, _bg_seq, file_posns, titles, classes = read_columns(fname)
        if posns is None:
            posns = file_posns
        elif file_posns != posns:
            raise InputError(
                "%s and %s do not number the same columns (%d and %d columns), "
                "so they cannot be shown in one table"
                % (fnames[0], fname, len(posns), len(file_posns))
            )
        rows.append((fname.removesuffix(".out"), fg_seq, titles, classes))
    pattern = pattern_marks(patternfname, posns, fnames[0])
    return ", ".join(fnames), tables(pattern, posns, rows)
