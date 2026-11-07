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

from cladecompare import gtest, hypg, jsd, phospho, urn
from cladecompare.shared import most_significant_type

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


def pair(
    module, fg_col, bg_col, fg_weights=None, bg_weights=None, fg_cons="D", freqs=FREQS
):
    fg_weights = fg_weights or ones(fg_col)
    bg_weights = bg_weights or ones(bg_col)
    return module.compare_cols(
        fg_col,
        fg_cons,
        math.fsum(fg_weights),
        fg_weights,
        bg_col,
        "E",
        math.fsum(bg_weights),
        bg_weights,
        freqs,
        PSEUDO,
    )


# --- hypg ---


def test_hypg_is_the_hypergeometric_tail_times_the_number_of_types():
    # 4 of 5 foreground D; 1 more D among 8 background sequences. The
    # foreground column holds D and E, so the smaller tail, D's, is doubled.
    expected = hypergeom_tail(4, 13, 5, 5)
    assert hypergeom_tail(1, 13, 6, 5) > expected  # E's tail
    assert pair(hypg, "DDDDE", "DEEEEEAA") == (
        "D",
        pytest.approx(2 * float(expected), rel=1e-12),
    )


def test_hypg_keeps_precision_for_small_tails():
    # 20 of 20 foreground D, none among 180 background: p = 1 / C(200, 20),
    # and a column of one type is not multiplied.
    expected = hypergeom_tail(20, 200, 20, 20)
    assert float(expected) < 1e-25
    assert pair(hypg, "D" * 20, "E" * 180) == (
        "D",
        pytest.approx(float(expected), rel=1e-9, abs=0),
    )


def test_hypg_rounds_weighted_counts_up():
    # Weighted: foreground D count 2.5 of size 3.5; background D count 0.9
    # of size 7.9. Rounded up: k = 3, N = 4, M = 4 + 8, and n = ceil(3.4)
    # = 4 successes, where rounding to nearest would give 3.
    fg_weights = [0.5, 0.5, 0.5, 1, 1]
    bg_weights = [0.9, 1, 1, 1, 1, 1, 1, 1]
    expected = hypergeom_tail(3, 12, 4, 4)
    assert expected == Fraction(33, 495)
    assert pair(hypg, "DDDDE", "DEEEEEAA", fg_weights, bg_weights) == (
        "D",
        pytest.approx(2 * float(expected), rel=1e-12),
    )


def test_hypg_tests_the_type_with_the_smallest_tail():
    """As in column 180 of the example pair: 6 A and 5 C in the foreground,
    where the consensus is A, against 24 A and no C among 36 background
    sequences. C is the type that sets the foreground apart."""
    fg = "A" * 6 + "C" * 5
    bg = "A" * 24 + "YYYSSVVFHI" + "--"
    expected = hypergeom_tail(5, 47, 5, 11)
    assert hypergeom_tail(6, 47, 30, 11) > expected  # A's tail
    assert pair(hypg, fg, bg, fg_cons="A") == (
        "C",
        pytest.approx(2 * float(expected), rel=1e-12),
    )


def test_hypg_tie_goes_to_the_consensus():
    # D and E have the same counts in both sets, so the same tail.
    residue, pvalue = pair(hypg, "DDEE", "AAAAAAAA", fg_cons="E")
    assert residue == "E"
    assert pvalue == pytest.approx(2 * float(hypergeom_tail(2, 12, 2, 4)), rel=1e-12)


def test_hypg_tie_without_the_consensus_goes_to_the_first_type():
    # D and E tie with the smallest tail; the consensus A is not among them.
    residue, pvalue = pair(hypg, "EEDDA", "AAAAAAAA", fg_cons="A")
    assert residue == "D"
    assert pvalue == pytest.approx(3 * float(hypergeom_tail(2, 13, 2, 5)), rel=1e-12)


def test_hypg_multiplied_tail_is_at_most_one():
    assert pair(hypg, "DE", "DDEE") == ("D", 1.0)


# --- most_significant_type ---


def test_most_significant_type_makes_a_nan_tail_nan():
    """A tail that cannot be computed is not hidden behind another type's."""
    for tails in [{"A": math.nan, "C": 0.1}, {"A": 0.1, "C": math.nan}]:
        residue, pvalue = most_significant_type({"A": 1, "C": 1}, "A", tails.get)
        assert residue == "A"
        assert math.isnan(pvalue)


# --- urn ---


def test_urn_pair_is_the_binomial_tail_of_the_pooled_frequency():
    # Pooled D frequency: background count with pseudocount, plus the
    # foreground count, over both sizes plus the pseudocount size. The
    # foreground column holds D and E, so the smaller tail, D's, is doubled.
    p_j = (1 + FREQS["D"] * PSEUDO + 4) / (8 + 5 + PSEUDO)
    p_e = (5 + FREQS["E"] * PSEUDO + 1) / (8 + 5 + PSEUDO)
    assert binom_tail(1, 5, p_e) > binom_tail(4, 5, p_j)
    assert pair(urn, "DDDDE", "DEEEEEAA") == (
        "D",
        pytest.approx(2 * binom_tail(4, 5, p_j), rel=1e-12),
    )


def test_urn_single_is_the_binomial_tail_of_the_overall_frequency():
    col = "DDDEE"
    assert binom_tail(2, 5, FREQS["E"]) > binom_tail(3, 5, FREQS["D"])
    assert urn.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == (
        "D",
        pytest.approx(2 * binom_tail(3, 5, FREQS["D"]), rel=1e-12),
    )


def test_urn_single_tests_the_type_with_the_smallest_tail():
    # A is the consensus, but 3 D of 7 is less likely than 4 A of 7, since D
    # is the rarer type overall.
    col = "AAAADDD"
    expected = binom_tail(3, 7, FREQS["D"])
    assert binom_tail(4, 7, FREQS["A"]) > expected
    assert urn.compare_one(col, "A", 7.0, ones(col), FREQS, PSEUDO) == (
        "D",
        pytest.approx(2 * expected, rel=1e-12),
    )


# --- phospho ---


def test_phospho_pair_is_the_binomial_tail_of_the_pooled_ser_thr_tyr():
    # Foreground S+T+Y count, without pseudocounts: 4. Background with
    # pseudocounts: 1 + (0.2 + 0.05 + 0.05) * 0.5 = 1.15. The set of types is
    # fixed in advance, so the tail is not multiplied.
    p_j = (1.15 + 4) / (8 + 5 + PSEUDO)
    assert pair(phospho, "SSTYE", "SEEEEEAA") == (
        "S",
        pytest.approx(binom_tail(4, 5, p_j), rel=1e-12),
    )


def test_phospho_single_is_the_binomial_tail_of_the_overall_ser_thr_tyr():
    col = "SSTEE"
    p_j = FREQS["S"] + FREQS["T"] + FREQS["Y"]
    assert phospho.compare_one(col, "S", 5.0, ones(col), FREQS, PSEUDO) == (
        "S",
        pytest.approx(binom_tail(3, 5, p_j), rel=1e-12),
    )


def test_phospho_column_without_ser_thr_tyr_is_not_significant():
    """No Ser/Thr/Tyr in the foreground gives 1, not the chance of one."""
    assert pair(phospho, "DDDEE", "SSSSSSAA") == ("D", 1.0)
    col = "DDDEE"
    assert phospho.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == ("D", 1.0)


def test_phospho_single_counts_an_absent_ser_thr_tyr_type_as_frequency_zero():
    """The alignment-wide frequencies have no entry for an absent type."""
    col = "DDDEE"
    freqs = {"A": 0.5, "D": 0.3, "E": 0.2}
    assert phospho.compare_one(col, "D", 5.0, ones(col), freqs, PSEUDO) == ("D", 1.0)
    col = "TTEEE"
    freqs = {"E": 0.5, "T": 0.3, "Y": 0.2}
    assert phospho.compare_one(col, "E", 5.0, ones(col), freqs, PSEUDO) == (
        "T",
        pytest.approx(binom_tail(2, 5, 0.5), rel=1e-12),
    )


@pytest.mark.parametrize(
    "fg_col, residue",
    [
        ("DDTTY", "T"),  # the consensus D is not what phospho tests
        ("DSSTT", "S"),  # a tie goes to the alphabetically first
        ("DYYYS", "Y"),
    ],
)
def test_phospho_names_the_site_after_its_most_frequent_ser_thr_tyr(fg_col, residue):
    assert pair(phospho, fg_col, "DEEEEEAA")[0] == residue


# --- jsd ---


def entropy(dist):
    return -math.fsum(p * log2(p) for p in dist.values() if p)


def normalized(counts, letters=AMINO):
    total = math.fsum(counts.values())
    return {aa: counts.get(aa, 0.0) / total for aa in letters}


def divergence(p, q, letters=AMINO):
    """JSD in bits with equal mixture weights."""
    mix = {aa: 0.5 * p[aa] + 0.5 * q[aa] for aa in letters}
    return entropy(mix) - 0.5 * (entropy(p) + entropy(q))


def test_jsd_pair_compares_with_the_pseudocounted_background():
    p = normalized({"D": 4, "E": 1})
    q = normalized(
        {
            aa: {"D": 1, "E": 5, "A": 2}.get(aa, 0) + f * PSEUDO
            for aa, f in FREQS.items()
        }
    )
    assert pair(jsd, "DDDDE", "DEEEEEAA") == (
        "D",
        pytest.approx(divergence(p, q), rel=1e-12),
    )


def test_jsd_single_compares_with_uniform_frequencies():
    col = "DDDEA"
    p = normalized({"D": 3, "E": 1, "A": 1})
    q = {aa: 1 / 20 for aa in AMINO}
    assert jsd.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == (
        "D",
        pytest.approx(divergence(p, q), rel=1e-12),
    )


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="the two entropies of a uniform column differ by a rounding error, "
    "which gives a JSD of -8.9e-16",
)
def test_jsd_of_a_uniform_column_is_not_negative():
    # Equal weights of 1.0 cancel exactly; 0.3 leaves the two entropies
    # differing by a rounding error.
    assert jsd.compare_one(AMINO, "A", 6.0, [0.3] * 20, FREQS, PSEUDO)[1] >= 0.0


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="the counts include B, which the mixture of the 20 amino acids "
    "omits, so the JSD of 'BBDD' against 'DDDD' is -0.35",
)
def test_jsd_of_a_pair_with_an_ambiguity_code_lies_between_0_and_1():
    assert 0.0 <= pair(jsd, "BBDD", "DDDD")[1] <= 1.0


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="the counts include B, which the mixture of the 20 amino acids "
    "omits, so the JSD of 'BBBB' against uniform frequencies is 2.16",
)
def test_jsd_of_a_single_column_with_an_ambiguity_code_lies_between_0_and_1():
    col = "BBBB"
    assert 0.0 <= jsd.compare_one(col, "B", 4.0, ones(col), FREQS, PSEUDO)[1] <= 1.0


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="a B in the pseudocounts enters the background composition but not "
    "the mixture of the 20 amino acids, so a column without B is not scored as "
    "the divergence of its compositions",
)
def test_jsd_pair_with_an_ambiguity_code_in_the_pseudocounts_is_a_divergence():
    freqs = {aa: 0.9 * f for aa, f in FREQS.items()} | {"B": 0.1}
    letters = AMINO + "B"
    p = normalized({"D": 4, "E": 1}, letters)
    q = normalized(
        {
            aa: {"D": 1, "E": 5, "A": 2}.get(aa, 0) + f * PSEUDO
            for aa, f in freqs.items()
        },
        letters,
    )
    assert pair(jsd, "DDDDE", "DEEEEEAA", freqs=freqs)[1] == pytest.approx(
        divergence(p, q, letters), rel=1e-12
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
    assert pair(gtest, "DDDDE", "DEEEEEAA") == (
        "D",
        pytest.approx(chisq_tail(G, df), rel=1e-12),
    )
    assert chisq_tail(G, df) == pytest.approx(0.03155, abs=1e-5)


def test_gtest_pair_two_types_have_no_continuity_correction():
    """A 2 x 2 table gives the uncorrected G on 1 degree of freedom, where
    Yates' correction would give about 0.021."""
    G, df = homogeneity_g({"D": 8, "E": 3}, {"D": 10, "E": 26})
    assert df == 1
    assert pair(gtest, "D" * 8 + "E" * 3, "D" * 10 + "E" * 26) == (
        "D",
        pytest.approx(chisq_tail(G, df), rel=1e-12),
    )
    assert chisq_tail(G, df) == pytest.approx(0.0076, abs=1e-4)


def test_gtest_pair_ignores_gaps():
    """Gaps in either column change neither the table nor the p-value."""
    residue, pvalue = pair(gtest, "DDDDE", "DEEEEEAA")
    assert pair(gtest, "DD--DDE-", "DEEEE-EAA") == (
        residue,
        pytest.approx(pvalue, rel=1e-12),
    )


def test_gtest_pair_single_residue_type_is_not_significant():
    """One residue type in both columns leaves no degrees of freedom: p = 1."""
    assert pair(gtest, "DDD", "DDDDD-") == ("D", 1.0)


def test_gtest_single_compares_with_overall_frequencies_at_the_column_total():
    # Observed: the column's counts plus pseudocounts, 5 + 0.5 in total;
    # expected: the overall frequencies scaled to that total.
    col = "DDDEA"
    observed = {
        aa: {"D": 3, "E": 1, "A": 1}.get(aa, 0) + f * PSEUDO for aa, f in FREQS.items()
    }
    G = 2 * math.fsum(o * math.log(o / (5.5 * FREQS[aa])) for aa, o in observed.items())
    assert gtest.compare_one(col, "D", 5.0, ones(col), FREQS, PSEUDO) == (
        "D",
        pytest.approx(chisq_tail(G, 19), rel=1e-12),
    )


def test_gtest_single_column_with_the_overall_composition_is_not_significant():
    """A column whose composition equals the alignment's gives G = 0, p = 1,
    however many sequences it holds."""
    col = "AAAADDEEEEEEEESSSSTY"
    assert gtest.compare_one(col, "E", 20.0, ones(col), FREQS, PSEUDO) == (
        "E",
        pytest.approx(1.0, rel=1e-12),
    )
