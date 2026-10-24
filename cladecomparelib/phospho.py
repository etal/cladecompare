"""Compare two alignments for potential phosphorylation sites (Ser/Thr/Tyr).

Based on the "urn" model, but testing for conservation of Ser/Thr/Tyr instead of
the foreground consensus residue type.

"""
from __future__ import division

import logging
import math

from scipy.stats import binom

from .shared import count_col


def compare_cols(fg_col, fg_cons, fg_size, fg_weights,
                 bg_col, bg_cons, bg_size, bg_weights,
                 aa_freqs, pseudo_size):
    "Compare Ser/Thr/Tyr conservation using the ball-in-urn model."
    # Weighted count of Ser/Thr/Tyr in the foreground column, with pseudocounts
    fg_counts = count_col(fg_col, fg_weights, aa_freqs, pseudo_size)
    fg_tot = fg_counts['S'] + fg_counts['T'] + fg_counts['Y']
    # Ser/Thr/Tyr frequency in the combined foreground and background column
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    p_j = (bg_counts['S'] + bg_counts['T'] + bg_counts['Y'] + fg_tot
          ) / (bg_size + fg_size + 2.0) # pseudocount size = 1.0
    # The count carries fractional pseudocounts and the size sums weights;
    # range() and the binomial trial count need integers
    fg_tot_i = int(math.ceil(fg_tot))
    fg_size_i = int(math.ceil(fg_size))

    # Probability of fg col conservation vs. the combined/main set
    # (P_j_LB in the publication)
    # NB: Some tweaks for pseudocounts
    pvalue = binom.pmf(range(fg_tot_i, fg_size_i + 2),
                       fg_size_i + 1, p_j).sum()
    if pvalue == 1.0:
        logging.info("Meaningless p-value: p_j=%s, fg=%s vs. bg=%s",
                        p_j, fg_tot, bg_counts)

    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    "Column Ser/Thr/Tyr conservation probability using the ball-in-urn model."
    col_counts = count_col(col, weights, aa_freqs, pseudo_size)
    col_tot = col_counts['S'] + col_counts['T'] + col_counts['Y']
    p_j = (aa_freqs['S'] + aa_freqs['T'] + aa_freqs['Y'])
    # The count carries fractional pseudocounts and the size sums weights;
    # range() and the binomial trial count need integers
    col_tot_i = int(math.ceil(col_tot))
    size_i = int(math.ceil(aln_size))
    pvalue = binom.pmf(range(col_tot_i, size_i + 2), size_i + 1, p_j).sum()
    return pvalue
