"""Numerical functions shared by multiple CladeCompare modules."""

import collections
from math import log

from biofrills import alnutils


MAX_ENTROPY = log(20, 2)


def count_col(col, weights, pseudo=None, pseudo_size=1.0):
    """Get a weighted count of the characters in a column, skipping gaps.

    weights:
        Sequence weight coefficients. List output of alnutils.sequence_weights.
    pseudo:
        Chacter pseudocounts. Dict output of combined_frequencies.
    """
    counts = collections.defaultdict(float)
    for char, wt in zip(col, weights):
        if char in '-.X':
            continue
        counts[char] += wt
    if pseudo:
        # Apply pseudocounts
        for char, pcount in pseudo.iteritems():
            counts[char] += pcount * pseudo_size
    return counts


def entropy(counts):
    """Entropy of a set of observations, in the information-theoretic sense."""
    scale = 1.0 / sum(counts)
    probs = [c * scale for c in counts if c]
    # try:
    return -sum(p * log(p, 2) for p in probs)
    # except ValueError:
    #     print locals()
    #     return 0.


def combined_frequencies(fg_aln, fg_weights, bg_aln, bg_weights):
    """Get overall frequencies from FG + BG alignments (for pseudocounts).

    Returns a dict of characters and their overall frequencies.
    """
    freqs = alnutils.aa_frequencies(list(fg_aln) + list(bg_aln),
                                    fg_weights + bg_weights,
                                    gap_chars='-.X')
    assert abs(sum(freqs.values()) - 1.0) < 1e-6, \
            "Freqs sum to %s:\n%s" % (sum(freqs.values()), freqs)
    return freqs


def scale_counts(counts):
    """Scale a dict of keys-to-integer counts so values sum to 1."""
    scale = 1.0 / sum(counts.values())
    return dict((key, cnt * scale) for key, cnt in counts.iteritems())

