#!/usr/bin/env python

"""Visualize noise file(s) in an HTML page.

If given one file:
    Display foreground and background consensus sequences as 2 rows, wrapped.
    Color of foreground cells corresponds to p-value (#stars).

If given multiple files:
    Display all foreground consensus sequences as N rows, wrapped.
    Background is [undecided] either the bottom row, if the same across all
    input seqs, or a separate block, or just not shown.
    Color of foreground cells corresponds to p-value (#stars).

"""

import sys

import logging
from math import log10

from esbglib import sugar
sugar.log_config()


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
    td.b3 { background-color: #99f }
    td.b2 { background-color: #bbf }
    td.b1 { background-color: #ddf }
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


def parse_noise(fname):
    """Iterate over the data for each site in the '.noise' file."""
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
        idx, site_ln = map(str.strip, line.split(':'))
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
    nlp = -log10(pval)
    if nlp < 0.1:  return 'b3'
    if nlp <= 0.3: return 'b2'
    if nlp <= 1.0: return 'b1'
    if nlp <= 2.0: return 'wt'
    if nlp >= 10.0: return 'r9'
    # r1, r2, ..., r8
    return 'r%d' % int(nlp-2)


# Formatting

def table_single(classes, pattern, foreground, background, positions):
    assert len(foreground) == len(background) == len(positions)
    out = "<table>"
    if any(pattern):
        out += tr_plain("", pattern)
    out += tr_class("foreground", foreground, classes)
    out += tr_plain("background", background)
    out += tr_pos(positions)
    out += "</table>"
    return out


def table_multi(labels, rowcells, rowclasses, pattern, positions):
    assert len(labels) == len(rowclasses) == len(rowcells)
    assert len(positions) == len(rowclasses[0]) == len(rowcells[0])
    out = "<table>"
    if any(pattern):
        out += tr_plain("", pattern)
    for label, cells, classes in zip(labels, rowcells, rowclasses):
        out += tr_class(label, cells, classes)
    out += tr_pos(positions)
    out += "</table>"
    return out


def tr_plain(label, cells):
    out = "<tr>\n"
    out += "<td class='label'>%s</td>" % label
    out += ''.join(['<td>%s</td>' % c for c in cells])
    out += "\n</tr>"
    return out


def tr_class(label, cells, classes):
    out = "<tr>\n"
    out += "<td class='label'>%s</td>" % label
    out += ''.join(['<td class="%s">%s</td>' % (cls, txt)
                    for cls, txt in zip(classes, cells)])
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


def do_single(noisefname, patternfname):
    """Build an HTML document for a single 'noise' file and optional pattern."""
    sites = list(parse_noise(noisefname))
    fg_seq, bg_seq, posns, probs, _nstarses = zip(*sites)
    classes = [p2class(fg, bg, p) for fg, bg, p in zip(fg_seq, bg_seq, probs)]
    posns = [str(p) if not p % 10 else '' for p in posns]
    pattern = [''] * len(posns)
    if patternfname:
        for posn in parse_pattern(patternfname):
            pattern[posn-1] = '*'
    # wrap it up
    contents = ''
    chunked_data = map(wrap_data, [classes, pattern, fg_seq, bg_seq, posns])
    for layer in zip(*chunked_data):
        contents += table_single(*layer)
        contents += "<p />"
    return noisefname, contents


def do_multi(fnames, patternfname):
    """Build an HTML document combining multiple 'noise' files."""
    rowcells = []   # dim = #fnames x #cols
    rowclasses = []
    for fname in fnames:
        sites = list(parse_noise(fname))
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
    labels = [fn[:-len('.noise')] if fn.endswith('.noise') else fn
              for fn in fnames]
    chunked_posns = wrap_data(posns)
    chunked_pattern = wrap_data(pattern)
    wrapped_row_cells = zip(*map(wrap_data, rowcells))
    wrapped_row_classes = zip(*map(wrap_data, rowclasses))
    for w_rowcells, w_rowclasses, ptnrow, positions in zip(
        wrapped_row_cells, wrapped_row_classes, chunked_pattern, chunked_posns):
        contents += table_multi(labels, w_rowcells, w_rowclasses, ptnrow, positions)
        contents += "<p />"
    return ', '.join(fnames), contents



if __name__ == '__main__':
    import argparse 
    AP = argparse.ArgumentParser(__doc__)
    AP.add_argument('noisefiles', nargs='+',
                    help="Noise file(s) generated by cladecompare.py")
    AP.add_argument('-p', '--pattern', 
                    help="Pattern file (.pttrn) generated by cladecompare.py")
    args = AP.parse_args()

    if len(args.noisefiles) == 1:
        title, contents = do_single(args.noisefiles[0], args.pattern)
    else:
        title, contents = do_multi(args.noisefiles, args.pattern)
    print html_page_tpl % dict(title=title, contents=contents)

