"""Numerical functions shared by multiple CladeCompare modules."""

import collections
from math import fsum, isnan, log, nan

from biofrills import alnutils

MAX_ENTROPY = log(20, 2)

# Characters that the counts and the alignment-wide frequencies skip: the gaps
# and X, an unknown residue. A column whose consensus is one of them has no
# residue type to test.
NON_RESIDUES = "-.X"


def count_col(col, weights, pseudo=None, pseudo_size=1.0):
    """Get a weighted count of the characters in a column, skipping gaps.

    weights:
        Sequence weight coefficients. List output of alnutils.sequence_weights.
    pseudo:
        Character pseudocounts. Dict output of combined_frequencies.
    """
    counts = collections.defaultdict(float)
    for char, wt in zip(col, weights, strict=True):
        if char in NON_RESIDUES:
            continue
        counts[char] += wt
    if pseudo:
        # Apply pseudocounts
        for char, pcount in pseudo.items():
            counts[char] += pcount * pseudo_size
    return counts


def most_significant_type(counts, consensus, tail):
    """The residue type of a column with the smallest tail probability, and
    that probability multiplied by the number of types, at most 1.

    counts:
        Counts of the column, as count_col returns them; its keys, the types
        that the column holds, are the candidates.
    consensus:
        The column's consensus residue, which wins a tie.
    tail:
        Function of a type that returns its tail probability.

    The consensus residue of a column is itself chosen from the types that the
    column holds, so a test of a type chosen from the column is one of K tests,
    K being the number of types; multiplying the smallest of their values by K
    (Bonferroni) keeps the result a p-value. Ties go to the consensus, then to
    the alphabetically first type. A tail that is nan makes the result nan, so
    that it reaches the output rather than losing the comparison.
    """
    candidates = [(tail(aa), aa != consensus, aa) for aa in counts]
    if any(isnan(pvalue) for pvalue, _, _ in candidates):
        return consensus, nan
    pvalue, _, residue = min(candidates)
    return residue, min(1.0, len(counts) * pvalue)


def entropy(counts):
    """Entropy of a set of observations, in the information-theoretic sense."""
    scale = 1.0 / fsum(counts)
    probs = [c * scale for c in counts if c]
    return -fsum(p * log(p, 2) for p in probs)


def combined_frequencies(fg_aln, fg_weights, bg_aln, bg_weights):
    """Get overall frequencies from FG + BG alignments (for pseudocounts).

    Returns a dict of characters and their overall frequencies.
    """
    freqs = alnutils.aa_frequencies(
        list(fg_aln) + list(bg_aln), fg_weights + bg_weights, gap_chars=NON_RESIDUES
    )
    assert abs(fsum(freqs.values()) - 1.0) < 1e-6, "Freqs sum to %s:\n%s" % (
        fsum(freqs.values()),
        freqs,
    )
    return freqs


def scale_counts(counts):
    """Scale a dict of keys-to-integer counts so values sum to 1."""
    scale = 1.0 / fsum(counts.values())
    return dict((key, cnt * scale) for key, cnt in counts.items())
