"""Tests of each strategy's column statistics against hand computations.

The expected values are computed here from the definitions, in exact rational
arithmetic where the statistic is a sum of probabilities, and for the G-test
with the closed form of the chi-squared tail at odd degrees of freedom, so that
a test does not share its numerical method with the code it checks.
"""

import math
from fractions import Fraction
from math import comb, log2

import pytest

from cladecomparelib import gtest, hypg, jsd, phospho, urn

AMINO = "ACDEFGHIKLMNPQRSTVWY"
PSEUDO = 0.5
# Background amino-acid frequencies for the pseudocounts.
FREQS = {"A": 0.2, "D": 0.1, "E": 0.4, "S": 0.2, "T": 0.05, "Y": 0.05}


def ones(col):
    return [1] * len(col)


def binom_tail(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    return math.fsum(comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1))


def hypergeom_tail(k, M, n, N):
    """P(X >= k) for N draws from M items of which n are successes."""
    return Fraction(
        sum(comb(n, i) * comb(M - n, N - i) for i in range(k, min(n, N) + 1)),
        comb(M, N),
    )


def chisq_tail(x, df):
    """P(X >= x) for X ~ chi-squared with odd df >= 3 (Abramowitz & Stegun 26.4.5)."""
    chi = math.sqrt(x)
    term = chi
    terms = [term]
    for r in range(2, (df - 1) // 2 + 1):
        term *= x / (2 * r - 1)
        terms.append(term)
    return math.erfc(chi / math.sqrt(2)) + 2 * math.exp(-x / 2) / math.sqrt(
        2 * math.pi
    ) * math.fsum(terms)


def pair(module, fg_col, bg_col, fg_weights=None, bg_weights=None):
    fg_weights = fg_weights or ones(fg_col)
    bg_weights = bg_weights or ones(bg_col)
    return module.compare_cols(
        fg_col,
        "D",
        math.fsum(fg_weights),
        fg_weights,
        bg_col,
        "E",
        math.fsum(bg_weights),
        bg_weights,
        FREQS,
        PSEUDO,
    )


# --- hypg ---


def test_hypg_is_the_hypergeometric_tail():
    # 4 of 5 foreground D; 1 more D among 8 background sequences.
    expected = hypergeom_tail(4, 13, 5, 5)
    assert pair(hypg, "DDDDE", "DEEEEEAA") == pytest.approx(float(expected), rel=1e-12)


def test_hypg_keeps_precision_for_small_tails():
    # 20 of 20 foreground D, none among 180 background: p = 1 / C(200, 20).
    expected = hypergeom_tail(20, 200, 20, 20)
    assert float(expected) < 1e-25
    assert pair(hypg, "D" * 20, "E" * 180) == pytest.approx(
        float(expected), rel=1e-9, abs=0
    )


def test_hypg_rounds_weighted_counts_up():
    # Weighted: foreground D count 2.5 of size 3.5; background D count 0.9
    # of size 7.9. Rounded up: k = 3, N = 4, M = 4 + 8, and n = ceil(3.4)
    # = 4 successes, where rounding to nearest would give 3.
    fg_weights = [0.5, 0.5, 0.5, 1, 1]
    bg_weights = [0.9, 1, 1, 1, 1, 1, 1, 1]
    expected = hypergeom_tail(3, 12, 4, 4)
    assert expected == Fraction(33, 495)
    assert pair(hypg, "DDDDE", "DEEEEEAA", fg_weights, bg_weights) == pytest.approx(
        float(expected), rel=1e-12
    )


def test_hypg_single_rejects_a_consensus_it_cannot_count():
    """A consensus of X raises instead of scoring 0, the most significant value.

    The alignment-wide frequencies skip X, as count_col does.
    """
    col = "XXXAD"
    with pytest.raises(KeyError):
        hypg.compare_one(col, "X", 5.0, ones(col), FREQS, PSEUDO)


# --- urn ---


def test_urn_pair_is_the_binomial_tail_of_the_pooled_frequency():
    # Pooled D frequency: background count with pseudocount, plus the
    # foreground count, over both sizes plus the pseudocount size.
    p_j = (1 + FREQS["D"] * PSEUDO + 4) / (8 + 5 + PSEUDO)
    assert pair(urn, "DDDDE", "DEEEEEAA") == pytest.approx(
        binom_tail(4, 5, p_j), rel=1e-12
    )


def test_urn_single_is_the_binomial_tail_of_the_overall_frequency():
    col = "DDDEE"
    assert urn.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == pytest.approx(
        binom_tail(3, 5, FREQS["D"]), rel=1e-12
    )


# --- phospho ---


def test_phospho_pair_is_the_binomial_tail_of_the_pooled_ser_thr_tyr():
    # Foreground S+T+Y count, without pseudocounts: 4. Background with
    # pseudocounts: 1 + (0.2 + 0.05 + 0.05) * 0.5 = 1.15.
    p_j = (1.15 + 4) / (8 + 5 + PSEUDO)
    assert pair(phospho, "SSTYE", "SEEEEEAA") == pytest.approx(
        binom_tail(4, 5, p_j), rel=1e-12
    )


def test_phospho_single_is_the_binomial_tail_of_the_overall_ser_thr_tyr():
    col = "SSTEE"
    p_j = FREQS["S"] + FREQS["T"] + FREQS["Y"]
    assert phospho.compare_one(
        col, "S", 5.0, ones(col), FREQS, PSEUDO
    ) == pytest.approx(binom_tail(3, 5, p_j), rel=1e-12)


def test_phospho_column_without_ser_thr_tyr_is_not_significant():
    """No Ser/Thr/Tyr in the foreground gives 1, not the chance of one."""
    assert pair(phospho, "DDDEE", "SSSSSSAA") == 1.0
    col = "DDDEE"
    assert phospho.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == 1.0


# --- jsd ---


def entropy(dist):
    return -math.fsum(p * log2(p) for p in dist.values() if p)


def normalized(counts):
    total = math.fsum(counts.values())
    return {aa: counts.get(aa, 0.0) / total for aa in AMINO}


def similarity(p, q):
    """(1 - JSD)^2 with equal mixture weights."""
    mix = {aa: 0.5 * p[aa] + 0.5 * q[aa] for aa in AMINO}
    return (1 - (entropy(mix) - 0.5 * (entropy(p) + entropy(q)))) ** 2


def test_jsd_pair_compares_with_the_pseudocounted_background():
    p = normalized({"D": 4, "E": 1})
    q = normalized(
        {
            aa: {"D": 1, "E": 5, "A": 2}.get(aa, 0) + f * PSEUDO
            for aa, f in FREQS.items()
        }
    )
    assert pair(jsd, "DDDDE", "DEEEEEAA") == pytest.approx(similarity(p, q), rel=1e-12)


def test_jsd_single_compares_with_uniform_frequencies():
    col = "DDDEA"
    p = normalized({"D": 3, "E": 1, "A": 1})
    q = {aa: 1 / 20 for aa in AMINO}
    assert jsd.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == pytest.approx(
        similarity(p, q), rel=1e-12
    )


# --- gtest ---


def test_gtest_pair_compares_with_the_scaled_background():
    # Expected counts: the pseudocounted background frequencies, scaled to
    # the foreground size; 19 degrees of freedom.
    bg = {
        aa: n + FREQS.get(aa, 0.0) * PSEUDO
        for aa, n in {"D": 1, "E": 5, "A": 2}.items()
    }
    expected = {aa: 5 * bg[aa] / (8 + PSEUDO) for aa in bg}
    G = 2 * math.fsum(
        o * math.log(o / expected[aa]) for aa, o in {"D": 4, "E": 1}.items()
    )
    assert pair(gtest, "DDDDE", "DEEEEEAA") == pytest.approx(
        chisq_tail(G, 19), rel=1e-12
    )


def test_gtest_single_compares_with_overall_frequencies_at_the_column_total():
    # Observed: the column's counts plus pseudocounts, 5 + 0.5 in total;
    # expected: the overall frequencies scaled to that total.
    col = "DDDEA"
    observed = {
        aa: {"D": 3, "E": 1, "A": 1}.get(aa, 0) + f * PSEUDO for aa, f in FREQS.items()
    }
    G = 2 * math.fsum(o * math.log(o / (5.5 * FREQS[aa])) for aa, o in observed.items())
    assert gtest.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == pytest.approx(
        chisq_tail(G, 19), rel=1e-12
    )


def test_gtest_single_column_with_the_overall_composition_is_not_significant():
    """A column whose composition equals the alignment's gives G = 0, p = 1,
    however many sequences it holds."""
    col = "AAAADDEEEEEEEESSSSTY"
    assert gtest.compare_one(col, "E", 20.0, ones(col), FREQS, PSEUDO) == pytest.approx(
        1.0, rel=1e-12
    )
