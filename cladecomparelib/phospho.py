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
    fg_counts = count_col(fg_col, fg_weights)
    fg_tot = fg_counts["S"] + fg_counts["T"] + fg_counts["Y"]
    # Ser/Thr/Tyr frequency in the combined foreground and background column,
    # with the background's pseudocounts
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    p_j = (bg_counts["S"] + bg_counts["T"] + bg_counts["Y"] + fg_tot) / (
        bg_size + fg_size + pseudo_size
    )
    return sty_tail(fg_tot, fg_size, p_j)


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    "Column Ser/Thr/Tyr conservation probability using the ball-in-urn model."
    col_counts = count_col(col, weights)
    col_tot = col_counts["S"] + col_counts["T"] + col_counts["Y"]
    p_j = aa_freqs["S"] + aa_freqs["T"] + aa_freqs["Y"]
    return sty_tail(col_tot, aln_size, p_j)
