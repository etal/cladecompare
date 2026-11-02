"""Compare two alignments for potential phosphorylation sites (Ser/Thr/Tyr).

Based on the "urn" model, but testing for conservation of Ser/Thr/Tyr instead of
the foreground consensus residue type. As in urn, the foreground count has no
pseudocounts; they enter only the background counts of the pooled frequency.

Under -w, the weighting follows CHAIN, as in urn: the foreground is not
weighted, and the background is weighted with biofrills' "andy" scaling. With
a single alignment, -w changes nothing.

"""

import math

from scipy.stats import binom

from .shared import count_col

# Foreground and background scalings of the sequence weights under -w.
WEIGHT_SCALINGS = (None, "andy")


def sty_total(counts):
    """The summed counts or frequencies of Ser, Thr and Tyr; an absent type
    counts as 0."""
    return counts.get("S", 0.0) + counts.get("T", 0.0) + counts.get("Y", 0.0)


def sty_tail(count, size, p_j):
    """P(X >= count) for X ~ Binomial(size, p_j), with count and size rounded up.

    A column without Ser/Thr/Tyr gives exactly 1.
    """
    return binom.sf(math.ceil(count) - 1, math.ceil(size), p_j)


def compare_cols(
    fg_col,
    fg_cons,
    fg_size,
    fg_weights,
    bg_col,
    bg_cons,
    bg_size,
    bg_weights,
    aa_freqs,
    pseudo_size,
):
    "Compare Ser/Thr/Tyr conservation using the ball-in-urn model."
    fg_tot = sty_total(count_col(fg_col, fg_weights))
    # Ser/Thr/Tyr frequency in the combined foreground and background column,
    # with the background's pseudocounts
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    p_j = (sty_total(bg_counts) + fg_tot) / (bg_size + fg_size + pseudo_size)
    return sty_tail(fg_tot, fg_size, p_j)


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Column Ser/Thr/Tyr conservation probability using the ball-in-urn model.

    A type absent from the alignment has frequency 0; an alignment without any
    of Ser, Thr and Tyr therefore gives 1 in every column.
    """
    col_tot = sty_total(count_col(col, weights))
    p_j = sty_total(aa_freqs)
    return sty_tail(col_tot, aln_size, p_j)
