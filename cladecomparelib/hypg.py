"""Write up for hypergeometric.

CHAIN formula:

The selective constraint acting at position j in a subalignment is expressed in
terms of the number of random trials needed to draw from among the residues in
a superalignment (with replacement) at least as many conserved residues as are
observed in the subalignment at that position.

    P_j^{(L,B} = \sum_{i=c_j^{(L)}^{N_j^{(L)}
                    \binom{N_j}{i}
                    (p_j^{(B)})^i 
                    (1-p_j^{(B)})^{N_j^{(L)}-i}

where c_j^{(L)} and N_j^{(L)} are the number of conserved residues and total
number of residues, respectively, in the j^th column of subalignment L, and
p_j^{(B)} is the frequency of the conserved residues observed at that position
for superalignment B, which serves as the background model.

Neuwald: "Note that weights are not computed for the query family alignment,
because these sequences are selected from distinct phyla or kingdoms and,
therefore, are treated as statistically independent."

The corresponding selective constraint acting on subalignment L is then defined
as

    K_j^{(L,B)} = 1 / P_j^{(L,B)}

the expected number of random trials needed to observed this event.
(e.g. P = 0.01 => 100 trials)

Histogram bar height ~ number of random trials implied by K (i.e., K).

A hack of logarithmic scaling:

    h = (t^{1-sigma}) / (1 - sigma)

where
    t = number of random trials
    sigma \in [0,1) is a scaling parameter for adjusting the relative bar
    heights so as to converge to linear scaling at sigma=0 and logarithmic
    scaling as sigma->1. (Automatically determined by the display program)

The order-of-magnitude increase in t as a function of sigma, when he relative
bar height increases by twofold, is given by
    log_10 (t_2h / t_h) = log_10 (2^{1/(1-sigma)})

"""
from math import ceil

from scipy.stats import hypergeom
from biofrills import consensus, alnutils

from .shared import count_col, combined_frequencies


def compare_cols(fg_col, fg_cons, fg_size, fg_weights,
                 bg_col, bg_cons, bg_size, bg_weights,
                 aa_freqs, pseudo_size):
    "Compare alignments using the hypergeometric model"
    # Number of consensus-type residues in the foreground column
    fg_cons_count = count_col(fg_col, fg_weights)[fg_cons]
    # Consensus residue frequency in the combined alignment column
    p_j = count_col(bg_col, bg_weights)[fg_cons] + fg_cons_count
    # Round fg counts & size to nearest integer for hypergeometric test
    fg_cons_count_i = max(1, int(ceil(fg_cons_count)))
    fg_size_i = int(ceil(fg_size))
    bg_size_i = int(ceil(bg_size))
    # Probability of fg col conservation vs. the combined/main set
    pvalue = 1-hypergeom.cdf(fg_cons_count_i-1,fg_size_i+bg_size_i,
                        p_j, fg_size_i)
    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    "Column probability using the hypergeometric model."
    # cons_count = col.count(cons_aa)
    cons_count = count_col(col, weights)[cons_aa]
    cons_count_i = int(ceil(cons_count))
    size_i = int(ceil(aln_size))
    pvalue = 1-hypergeom.cdf(cons_count_i-1,size_i,
                        cons_count_i, size_i)
    return pvalue

