"""Generate an HTML report summarizing CladeCompare output."""

import logging
from math import log10


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


def parse_out(fname):
    """Iterate over the data for each site in the '.output' file."""
    with open(fname) as infile:
        for line in infile:
            if not line.strip():
                continue
            tokens = line.split()
            fg_aa = tokens[0]
            bg_aa = tokens[1][1]
            posn = int(tokens[2])
            prob = float(tokens[4].split('=')[1])
            nstars = (len(tokens[5]) if len(tokens) == 6 else 0)
            yield (fg_aa, bg_aa, posn, prob, nstars)


def parse_pattern(fname):
    with open(fname) as infile:
        # NB: only reading the first pattern, ignoring any others
        line = infile.readline()
        _idx, site_ln = map(str.strip, line.split(':'))
        if not site_ln:
            logging.warn("Empty pattern file: %s", fname)
            return []
        sites = site_ln.split(',')
        # Convert to number...
        posns = sorted(int(''.join(filter(str.isdigit, site)))
                       for site in sites)
        return list(posns)


def p2class(fg, bg, pval):
    """Assign a table cell (td) CSS class given a p-value.

    Heat map:
        higher p-values -> more blue
        p-values near significance -> white
        smaller p-values -> more red

    Foreground deletions -> gray, inserts -> gold.

    """
    if fg == '-':
        return 'del'
    if bg == '-':
        return 'ins'
    try:
        nlp = -log10(pval)
    except ValueError:
        # pval is 0, but let log10 decide that
        return 'r9'
    if nlp < 0.00001:  return 'b5'
    if nlp < 0.0005:  return 'b4'
    if nlp <= 0.005: return 'b3'
    if nlp <= 0.05: return 'b2'
    if nlp <= 0.3: return 'b1'
    if nlp <= 1.0: return 'wt'
    if nlp >= 9.0: return 'r9'
    # r1, r2, ..., r8
    return 'r%d' % int((nlp+1)/2)


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
    assert len(labels) == len(rowclasses) == len(rowcells)
    assert len(positions) == len(rowclasses[0]) == len(rowcells[0])
    out = "<table>"
    if any(pattern):
        out += tr_plain("", pattern)
    for label, cells, pvals in zip(labels, rowcells, rowclasses, rowpvalues):
        out += tr_class(label, cells, pvals, classes)
    out += tr_pos(positions)
    out += "</table>"
    return out


def tr_plain(label, cells):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % label
    out += ''.join(['<td>%s</td>' % c for c in cells])
    out += "\n</tr>"
    return out


def tr_class(label, cells, pvalues, classes):
    out = "<tr>\n"
    out += '<td class="label">%s</td>' % label
    out += ''.join(['<td class="%s" title="p=%s">%s</td>' % (cls, pval, txt)
                    for cls, pval, txt in zip(classes, pvalues, cells)])
    out += "\n</tr>"
    return out


def tr_pos(cells):
    out = "<tr class='pos'>\n<td></td>\n"
    out += ''.join(['<td>%s</td>' % c for c in cells])
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
    sites = list(parse_out(fname))
    fg_seq, bg_seq, posns, probs, _nstarses = zip(*sites)
    classes = [p2class(fg, bg, p) for fg, bg, p in zip(fg_seq, bg_seq, probs)]
    posns = [str(p) if not p % 10 else '' for p in posns]
    pattern = [''] * len(posns)
    if patternfname:
        for posn in parse_pattern(patternfname):
            pattern[posn-1] = '*'
    # wrap it up
    contents = ''
    chunked_data = map(wrap_data,
                       [probs, classes, pattern, fg_seq, bg_seq, posns])
    for layer in zip(*chunked_data):
        contents += table_single(*layer)
        contents += "<p />"
    return fname, contents


def do_multi(fnames, patternfname):
    """Build an HTML document combining multiple '.out' files."""
    rowcells = []   # dim = #fnames x #cols
    rowclasses = []
    for fname in fnames:
        sites = list(parse_out(fname))
        fg_seq, bg_seq, posns, probs, _nstarses = zip(*sites)
        rowcells.append(fg_seq)
        rowclasses.append([p2class(fg, bg, p)
                           for fg, bg, p in zip(fg_seq, bg_seq, probs)])
    posns = [str(p) if not p % 10 else '' for p in posns]
    pattern = [''] * len(posns)
    if patternfname:
        for posn in parse_pattern(patternfname):
            pattern[posn-1] = '*'
    # wrap it up
    contents = ''
    labels = [fn[:-len('.out')] if fn.endswith('.out') else fn
              for fn in fnames]
    chunked_posns = wrap_data(posns)
    chunked_pattern = wrap_data(pattern)
    wrapped_row_cells = zip(*map(wrap_data, rowcells))
    wrapped_row_classes = zip(*map(wrap_data, rowclasses))
    for w_rowcells, w_rowclasses, ptnrow, positions in zip(
        wrapped_row_cells, wrapped_row_classes, chunked_pattern, chunked_posns):
        contents += table_multi(labels, w_rowcells, w_rowclasses, ptnrow,
                                positions)
        contents += "<p />"
    return ', '.join(fnames), contents

