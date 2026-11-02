"""Tests of each strategy's column statistics against hand computations.

The expected values are computed here from the definitions, in exact rational
arithmetic where the statistic is a sum of probabilities, and for the G-test
with the closed forms of the chi-squared tail, so that a test does not share
its numerical method with the code it checks.
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
    """P(X >= x) for X ~ chi-squared, from the closed forms for even and odd df
    (Abramowitz & Stegun, section 26.4)."""
    if df % 2 == 0:
        term = 1.0
        terms = []
        for j in range(1, df // 2 + 1):
            terms.append(term)
            term *= x / (2 * j)
        return math.exp(-x / 2) * math.fsum(terms)
    chi = math.sqrt(x)
    term = chi
    terms = []
    for r in range(1, (df - 1) // 2 + 1):
        terms.append(term)
        term *= x / (2 * r + 1)
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


def test_phospho_single_counts_an_absent_ser_thr_tyr_type_as_frequency_zero():
    """The alignment-wide frequencies have no entry for an absent type."""
    col = "DDDEE"
    freqs = {"A": 0.5, "D": 0.3, "E": 0.2}
    assert phospho.compare_one(col, "D", 5.0, ones(col), freqs, PSEUDO) == 1.0
    col = "TTEEE"
    freqs = {"E": 0.5, "T": 0.3, "Y": 0.2}
    assert phospho.compare_one(
        col, "T", 5.0, ones(col), freqs, PSEUDO
    ) == pytest.approx(binom_tail(2, 5, 0.5), rel=1e-12)


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


def homogeneity_g(fg, bg):
    """G of the 2 x K table of two dicts of residue counts, and its df.

    Computed as 2 (sum O ln O - sum R ln R - sum C ln C + N ln N) over the
    cells O, the row totals R, the column totals C and the table total N.
    """
    rows = [math.fsum(fg.values()), math.fsum(bg.values())]
    cols = [fg.get(aa, 0) + bg.get(aa, 0) for aa in set(fg) | set(bg)]

    def xlogx(values):
        return math.fsum(v * math.log(v) for v in values if v)

    total = math.fsum(rows)
    G = 2 * (
        xlogx([*fg.values(), *bg.values()])
        - xlogx(rows)
        - xlogx(cols)
        + total * math.log(total)
    )
    return G, len(cols) - 1


def test_gtest_pair_is_the_g_test_of_homogeneity():
    # The 2 x 3 table D, E, A: foreground 4, 1, 0; background 1, 5, 2.
    G, df = homogeneity_g({"D": 4, "E": 1}, {"D": 1, "E": 5, "A": 2})
    assert df == 2
    assert pair(gtest, "DDDDE", "DEEEEEAA") == pytest.approx(
        chisq_tail(G, df), rel=1e-12
    )
    assert chisq_tail(G, df) == pytest.approx(0.03155, abs=1e-5)


def test_gtest_pair_two_types_have_no_continuity_correction():
    """A 2 x 2 table gives the uncorrected G on 1 degree of freedom, where
    Yates' correction would give about 0.021."""
    G, df = homogeneity_g({"D": 8, "E": 3}, {"D": 10, "E": 26})
    assert df == 1
    assert pair(gtest, "D" * 8 + "E" * 3, "D" * 10 + "E" * 26) == pytest.approx(
        chisq_tail(G, df), rel=1e-12
    )
    assert chisq_tail(G, df) == pytest.approx(0.0076, abs=1e-4)


def test_gtest_pair_ignores_gaps():
    """Gaps in either column change neither the table nor the p-value."""
    assert pair(gtest, "DD--DDE-", "DEEEE-EAA") == pytest.approx(
        pair(gtest, "DDDDE", "DEEEEEAA"), rel=1e-12
    )


def test_gtest_pair_single_residue_type_is_not_significant():
    """One residue type in both columns leaves no degrees of freedom: p = 1."""
    assert pair(gtest, "DDD", "DDDDD-") == 1.0


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
