"""Compare two alignments via CHAIN's ball-in-urn stats.

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
# ENH: Dirichlet mixture priors; Brown et al. 1993

from scipy.stats import binom
from biofrills import consensus, alnutils

from .shared import count_col, combined_frequencies


def compare_aln(fg_aln, bg_aln):
    """Compare alignments using the ball-in-urn model, like CHAIN does.
    """
    # BG seqs are weighted, FG seqs are not
    bg_weights = alnutils.sequence_weights(bg_aln, 'none')
    bg_size = sum(bg_weights)
    bg_cons = consensus.consensus(bg_aln, weights=bg_weights)
    # Height of the foreground alignment column
    fg_size = len(fg_aln)
    fg_cons = consensus.consensus(fg_aln)
    fg_cols = zip(*fg_aln)
    bg_cols = zip(*bg_aln)
    pseudocounts = combined_frequencies(fg_aln, [1]*fg_size,
                                        bg_aln, bg_weights)
    hits = []
    for faa, baa, fg_col, bg_col in zip(fg_cons, bg_cons, fg_cols, bg_cols):
        if faa == '-' or baa == '-':
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.
        else:
            # Cumulative binomial test
            # Number of consensus-type residues in the foreground column
            fg_cons_count = fg_col.count(faa)
            # Consensus residue frequency in the combined alignment column
            p_j = (count_col(bg_col, bg_weights, pseudocounts)[faa]
                   + fg_cons_count
                ) / (bg_size + fg_size + 1.0)
            # Probability of fg col conservation vs. the combined/main set
            # (P_j_LB in the publication)
            pvalue = binom.pmf(range(fg_cons_count, fg_size+1),
                               fg_size, p_j).sum()
        hits.append((faa, baa, pvalue))
    return hits

