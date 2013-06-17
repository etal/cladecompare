"""Generate a pair of logos with contrasting sites highlighted."""

import logging
import math
import os.path
from cStringIO import StringIO

from reportlab.pdfgen import canvas
from reportlab.lib import pagesizes, units

from biofrills import logoutils


UNIT = units.inch * 0.5  #= units.inch * 0.1
WUNIT = 0.6 * UNIT  # NB: Seems to be the text width-to-height ratio

# Try to find & register a nice monospace font
MONOSPACE = 'DejaVuSansMono-Bold'
font_loc = '/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansMono-Bold.ttf'
if not os.path.isfile(font_loc):
    # If the font isn't already installed, try to use the bundled copy
    font_loc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                            MONOSPACE) + '.ttf'
if os.path.isfile(font_loc):
    from reportlab.pdfbase import pdfmetrics, ttfonts
    font = ttfonts.TTFont(MONOSPACE, font_loc)
    pdfmetrics.registerFont(font)
    LINE_HT = 0.75
else:
    MONOSPACE = "Courier-Bold"
    logging.warn("Font 'DejaVuSansMono-Bold' unavailable; using lame Courier")
    LINE_HT = 0.6


aacolors = dict(
    A='black',
    C='green',
    D='red',
    E='red',
    F='black',
    G='green',
    H='blue',
    I='black',
    K='blue',
    L='black',
    M='black',
    N='purple',
    P='black',
    Q='purple',
    R='blue',
    S='green',
    T='green',
    V='black',
    W='black',
    Y='green',
)


def make_pairlogos(fg_aln, bg_aln, tophits, name, N):
    """
    """
    # Group together nearby "pattern sites"
    clusters = []
    last_seen = None
    for posn, faa, baa in sorted(tophits, key=lambda pfb: pfb[0]):
        # NB: posns are numbered from 1
        if last_seen and posn - last_seen <= N:
            clusters[-1].append(posn)
        else:
            clusters.append([posn])
        last_seen = posn
    # Make logos
    for sites in clusters:
        draw_pdf(name, fg_aln, bg_aln, sites, N)


def draw_pdf(name, fg_aln, bg_aln, sites, N):
    """Draw a PDF."""
    # Define regions for FG, BG logos
    left_wing = min(N, sites[0] - 1)
    right_wing = N #min(N, len(aln) - sites[-1])    # XXX TODO
    aln_start, aln_end = sites[0] - 1 - left_wing, sites[-1] + right_wing

    # -- Initialize page --
    pdfname = "%s-%s.pdf" % (name, ','.join(map(str, sites)))
    pg_width = (left_wing + right_wing + sites[-1] - sites[0] + 3) * WUNIT
    pg_height = 9 * UNIT
    can = canvas.Canvas(pdfname, (pg_width, pg_height))

    # -- Draw yellow bar(s) at site position(s) --
    site_x_offsets = [site - aln_start for site in sites]
    can.setFillColorRGB(1, 1, .1)
    for x_start in site_x_offsets:
        can.rect(x_start*WUNIT, 0.4*UNIT, WUNIT, 8.4*UNIT, stroke=0, fill=1)

    # -- Draw FG and BG logos --
    fg_ldata = logoutils.aln2logodata(fg_aln[:,aln_start:aln_end])
    bg_ldata = logoutils.aln2logodata(bg_aln[:,aln_start:aln_end])
    draw_logo(fg_ldata, can, WUNIT, UNIT*5)
    draw_logo(bg_ldata, can, WUNIT, UNIT)

    # -- Draw site label(s) at bottom of yellow bar(s) --
    can.setFillColorRGB(0, 0, 0)
    can.setFont("Helvetica", 0.35*UNIT)
    for site, x_start in zip(sites, site_x_offsets):
        can.drawCentredString((x_start + 0.5)*WUNIT, 0.5*UNIT, str(site))

    # -- Save PDF --
    can.showPage()
    can.save()
    logging.info("Wrote %s", pdfname)



def draw_logo(logodata, can, x, y):
    can.setFont(MONOSPACE, UNIT)
    # Calculate character height and width for each col
    x_offset = 0.0
    # TODO - deal w/ Q hanging below line -- raise & rescale height
    for posn, counts, entropy, weight in logodata:
        y_offset = 0.0
        for aa, frac in letter_scales(counts):
            can.saveState()
            scale = 1.6 * frac * math.sqrt(weight) * entropy
            can.scale(1.0, scale)
            can.setFillColor(aacolors[aa])
            can.drawString(x + x_offset, (y + y_offset)/scale, aa)
            can.restoreState()
            y_offset += scale * UNIT * LINE_HT
        x_offset += WUNIT


def letter_scales(counts):
    """Convert letter counts to frequencies, sorted increasing.""" 
    try:
        scale = 1.0 / sum(counts.values())
    except ZeroDivisionError:
        # This logo is all gaps, nothing can be done
        return []
    freqs = [(aa, cnt*scale) for aa, cnt in counts.iteritems() if cnt]
    freqs.sort(key=lambda pair: pair[1])
    return freqs


