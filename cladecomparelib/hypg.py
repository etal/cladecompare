"""Hypergeometric distribution.

Thought experiment:

    The foreground and background sequence sets are pre-defined.
    Given N foreground sequences and M - N background sequences, select N
    sequences at random from the M. Take the consensus residue of the
    foreground column as type I, and ask for the probability of observing at
    least as many type I sequences in the selection as the foreground has.
"""

from math import ceil

from scipy.stats import hypergeom

from .shared import count_col


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
    "Compare alignments using the hypergeometric model"
    # Number of consensus-type residues in the foreground column
    fg_cons_count = count_col(fg_col, fg_weights)[fg_cons]
    # Consensus residue count in the combined alignment column
    p_j = count_col(bg_col, bg_weights)[fg_cons] + fg_cons_count
    # The hypergeometric distribution is defined on integer counts, and
    # sequence weights make these counts non-integral, so every count and
    # size is rounded up. ceil is monotone, so the foreground count stays at
    # most the combined count, and ceil(a + b) <= ceil(a) + ceil(b), so the
    # combined count stays at most the population, the sum of the rounded
    # sizes. Rounding to nearest can exceed it: two sizes of 1.4 give a
    # population of 2, but a combined count of 2.8 rounds to 3, and SciPy
    # returns nan for more successes than the population holds.
    fg_cons_count_i = max(1, int(ceil(fg_cons_count)))
    p_j_i = int(ceil(p_j))
    fg_size_i = int(ceil(fg_size))
    bg_size_i = int(ceil(bg_size))
    # Probability of fg col conservation vs. the combined/main set: the upper
    # tail P(X >= k), computed directly, since 1 - cdf cancels for small tails
    pvalue = hypergeom.sf(fg_cons_count_i - 1, fg_size_i + bg_size_i, p_j_i, fg_size_i)
    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """The fraction of the column that matches the consensus residue.

    The weighted count of the consensus residue is rounded up and divided by
    the number of sequences, gaps included. A single alignment gives no
    second set to draw from, so this is not a hypergeometric probability, and
    not a p-value.
    """
    cons_count = count_col(col, weights)[cons_aa]
    cons_count_i = int(ceil(cons_count))
    pvalue = float(cons_count_i) / len(col)
    return pvalue
