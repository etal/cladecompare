r"""Compare two alignments via CHAIN's ball-in-urn stats.

The model is the selective-constraint statistic of CHAIN: Neuwald AF, Kannan
N, Poleksic A, Hata N, Liu JS (2003) Genome Research 13(4):673-692,
doi:10.1101/gr.862303, Methods, "CHAIN Analysis". For c conserved residues
among the N residues of column j of a subalignment, and p_j the frequency of
the conserved residue in column j of the superalignment,

    P_j = sum_{i=c..N} C(N, i) p_j^i (1 - p_j)^(N - i)

The section Methods of README.rst gives the formula as implemented here, with
its choice of subalignment, superalignment and conserved residue, and lists
where it departs from CHAIN.

CHAIN reports K_j = 1 / P_j, the expected number of random trials needed to
observe the event, and displays it as a histogram bar whose height scales
the number of trials t as

    h = (t^{1-sigma}) / (1 - sigma)

where sigma in [0, 1) moves the scaling from linear (sigma = 0) toward
logarithmic (sigma -> 1) and is chosen by the display program. When the
relative bar height doubles, t increases by log_10(2^{1/(1-sigma)}) orders
of magnitude. CladeCompare reports P_j and does not draw these bars.
"""

# ENH: Dirichlet mixture priors; Brown et al. 1993

from math import ceil

from scipy.stats import binom

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
    "Compare alignments using the ball-in-urn model (cumulative binomial test)"
    # Number of consensus-type residues in the foreground column
    fg_cons_count = count_col(fg_col, fg_weights)[fg_cons]
    # Consensus residue frequency in the combined alignment column
    p_j = (
        count_col(bg_col, bg_weights, aa_freqs, pseudo_size)[fg_cons] + fg_cons_count
    ) / (bg_size + fg_size + pseudo_size)
    # Round the fg count and size up to integers for the binomial test
    fg_cons_count_i = max(1, int(ceil(fg_cons_count)))
    fg_size_i = int(ceil(fg_size))
    # Probability of fg col conservation vs. the combined/main set
    # (P_j_LB in the CHAIN publication)
    pvalue = binom.pmf(range(fg_cons_count_i, fg_size_i + 1), fg_size_i, p_j).sum()
    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    "Column probability using the ball-in-urn model."
    cons_count = count_col(col, weights)[cons_aa]
    p_j = aa_freqs[cons_aa]
    cons_count_i = int(ceil(cons_count))
    size_i = int(ceil(aln_size))
    pvalue = binom.pmf(range(cons_count_i, size_i + 1), size_i, p_j).sum()
    return pvalue
