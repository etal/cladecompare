r"""Compare two alignments via CHAIN's ball-in-urn stats.

The model is the selective-constraint statistic of CHAIN: Neuwald AF, Kannan
N, Poleksic A, Hata N, Liu JS (2003) Genome Research 13(4):673-692,
doi:10.1101/gr.862303, Methods, "CHAIN Analysis". For c conserved residues
among the N residues of column j of a subalignment, and p_j the frequency of
the conserved residue in column j of the superalignment,

    P_j = sum_{i=c..N} C(N, i) p_j^i (1 - p_j)^(N - i)

The section Methods of README.rst gives the formula as implemented here, with
its choice of subalignment, superalignment and conserved residue, and lists
where it departs from CHAIN. The conserved residue is the type with the
smallest P_j among the types that the subalignment's column holds, and that
P_j is multiplied by the number of those types, since the type is chosen from
the column it is tested on.

Under -w, the weighting follows CHAIN: the foreground, CHAIN's query family,
is not weighted, because CHAIN takes its sequences to be independent, and the
background is weighted with CHAIN's robust weights, normalized to a mean of
one-half and truncated at one (biofrills' "andy" scaling). With a single
alignment, which plays the part of the query family, -w changes nothing.

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

from .shared import count_col, most_significant_type

# Foreground and background scalings of the sequence weights under -w.
WEIGHT_SCALINGS = (None, "andy")


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
    """Compare alignments using the ball-in-urn model (cumulative binomial test).

    Return the residue type tested and its p-value.
    """
    fg_counts = count_col(fg_col, fg_weights)
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    total = bg_size + fg_size + pseudo_size
    # Round the fg counts and size up to integers for the binomial test
    fg_size_i = int(ceil(fg_size))

    def tail(aa):
        # The type's frequency in the combined alignment column, and the
        # probability of fg col conservation vs. the combined/main set (P_j_LB
        # in the CHAIN publication): the upper tail P(X >= k), computed
        # directly, since a sum of point probabilities can exceed 1 by rounding
        p_j = (bg_counts[aa] + fg_counts[aa]) / total
        return binom.sf(int(ceil(fg_counts[aa])) - 1, fg_size_i, p_j)

    return most_significant_type(fg_counts, fg_cons, tail)


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Column probability using the ball-in-urn model.

    Return the residue type tested and its p-value.
    """
    counts = count_col(col, weights)
    size_i = int(ceil(aln_size))

    def tail(aa):
        return binom.sf(int(ceil(counts[aa])) - 1, size_i, aa_freqs[aa])

    return most_significant_type(counts, cons_aa, tail)
