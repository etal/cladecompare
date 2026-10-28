"""Generate a pair of logos with contrasting sites highlighted."""

import logging
import math
from importlib import resources

from biofrills import logoutils
from reportlab.lib import units
from reportlab.pdfbase import pdfmetrics, ttfonts
from reportlab.pdfgen import canvas

logger = logging.getLogger(__name__)

UNIT = units.inch * 0.5
WUNIT = 0.6 * UNIT  # NB: Seems to be the text width-to-height ratio

# The font ships as package data. as_file gives reportlab a real file even
# when the package is imported from a zip archive; the font is read in full
# while it is registered.
MONOSPACE = "DejaVuSansMono-Bold"
LINE_HT = 0.75
with resources.as_file(
    resources.files(__package__) / "data" / (MONOSPACE + ".ttf")
) as font_loc:
    pdfmetrics.registerFont(ttfonts.TTFont(MONOSPACE, font_loc))


aacolors = dict(
    A="black",
    C="green",
    D="red",
    E="red",
    F="black",
    G="green",
    H="blue",
    I="black",
    K="blue",
    L="black",
    M="black",
    N="purple",
    P="black",
    Q="purple",
    R="blue",
    S="green",
    T="green",
    V="black",
    W="black",
    Y="green",
)


def make_pairlogos(fg_aln, bg_aln, tophits, name, N):
    """Draw a pair logo for each run of pattern sites at most N columns apart."""
    # Group together nearby "pattern sites"
    clusters = []
    last_seen = None
    for posn in sorted(posn for posn, _faa, _baa in tophits):
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
    right_wing = N
    aln_start, aln_end = sites[0] - 1 - left_wing, sites[-1] + right_wing

    # -- Initialize page --
    pdfname = "%s-%s.pdf" % (name, ",".join(map(str, sites)))
    pg_width = (left_wing + right_wing + sites[-1] - sites[0] + 3) * WUNIT
    pg_height = 9 * UNIT
    can = canvas.Canvas(pdfname, (pg_width, pg_height))

    # -- Draw yellow bar(s) at site position(s) --
    site_x_offsets = [site - aln_start for site in sites]
    can.setFillColorRGB(1, 1, 0.1)
    for x_start in site_x_offsets:
        can.rect(x_start * WUNIT, 0.4 * UNIT, WUNIT, 8.4 * UNIT, stroke=0, fill=1)

    # -- Draw FG and BG logos --
    fg_ldata = logoutils.aln2logodata(fg_aln[:, aln_start:aln_end])
    bg_ldata = logoutils.aln2logodata(bg_aln[:, aln_start:aln_end])
    draw_logo(fg_ldata, can, WUNIT, UNIT * 5)
    draw_logo(bg_ldata, can, WUNIT, UNIT)

    # -- Draw site label(s) at bottom of yellow bar(s) --
    can.setFillColorRGB(0, 0, 0)
    can.setFont("Helvetica", 0.35 * UNIT)
    for site, x_start in zip(sites, site_x_offsets, strict=True):
        can.drawCentredString((x_start + 0.5) * WUNIT, 0.5 * UNIT, str(site))

    # -- Save PDF --
    can.showPage()
    can.save()
    logger.info("Wrote %s", pdfname)


def draw_logo(logodata, can, x, y):
    can.setFont(MONOSPACE, UNIT)
    # Calculate character height and width for each col
    x_offset = 0.0
    # TODO - deal w/ Q hanging below line -- raise & rescale height
    for _posn, counts, entropy, weight in logodata:
        y_offset = 0.0
        for aa, frac in letter_scales(counts):
            can.saveState()
            scale = 1.6 * frac * math.sqrt(weight) * entropy
            can.scale(1.0, scale)
            can.setFillColor(aacolors[aa])
            can.drawString(x + x_offset, (y + y_offset) / scale, aa)
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
    freqs = [(aa, cnt * scale) for aa, cnt in counts.items() if cnt]
    freqs.sort(key=lambda pair: pair[1])
    return freqs
