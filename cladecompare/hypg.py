"""Hypergeometric distribution.

Thought experiment:

    The foreground and background sequence sets are pre-defined.
    Given N foreground sequences and M - N background sequences, select N
    sequences at random from the M. For a residue type of the foreground
    column, ask for the probability of observing at least as many sequences
    of that type in the selection as the foreground has: the one-sided
    Fisher exact test of the 2 x 2 table (foreground, background) x (type,
    other). The type tested is the one with the smallest probability among
    the types that the foreground column holds, and that probability is
    multiplied by the number of those types, since the type is chosen from
    the column it is tested on.
"""

from math import ceil

from scipy.stats import hypergeom

from .shared import count_col, most_significant_type


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
    """Compare alignments using the hypergeometric model.

    Return the residue type tested and its p-value.
    """
    fg_counts = count_col(fg_col, fg_weights)
    bg_counts = count_col(bg_col, bg_weights)
    # The hypergeometric distribution is defined on integer counts, and
    # sequence weights make these counts non-integral, so every count and
    # size is rounded up. ceil is monotone, so the foreground count stays at
    # most the combined count, and ceil(a + b) <= ceil(a) + ceil(b), so the
    # combined count stays at most the population, the sum of the rounded
    # sizes. Rounding to nearest can exceed it: two sizes of 1.4 give a
    # population of 2, but a combined count of 2.8 rounds to 3, and SciPy
    # returns nan for more successes than the population holds.
    fg_size_i = int(ceil(fg_size))
    population = fg_size_i + int(ceil(bg_size))

    def tail(aa):
        # The upper tail P(X >= k), computed directly, since 1 - cdf cancels
        # for small tails
        fg_count_i = int(ceil(fg_counts[aa]))
        successes = int(ceil(bg_counts[aa] + fg_counts[aa]))
        return hypergeom.sf(fg_count_i - 1, population, successes, fg_size_i)

    return most_significant_type(fg_counts, fg_cons, tail)
