"""Compare two alignments for potential phosphorylation sites (Ser/Thr/Tyr).

Based on the "urn" model, but testing for conservation of Ser/Thr/Tyr instead of
the foreground consensus residue type.

"""

import logging
import math

from scipy.stats import binom
from biofrills import consensus, alnutils

from .shared import count_col, combined_frequencies


def compare_cols(fg_col, fg_cons, fg_size, fg_weights,
                 bg_col, bg_cons, bg_size, bg_weights,
                 aa_freqs, pseudo_size):
    "Compare alignments using the ball-in-urn model (cumulative binomial test)"
    # Number of consensus-type residues in the foreground column
    fg_counts = count_col(fg_col, fg_weights, aa_freqs, pseudo_size)
    fg_tot = fg_counts['S'] + fg_counts['T'] + fg_counts['Y']
    # Consensus residue frequency in the combined alignment column
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    p_j = (bg_counts['S'] + bg_counts['T'] + bg_counts['Y'] + fg_tot
          ) / (bg_size + fg_size + 2.0) # pseudocount size = 1.0

    # Probability of fg col conservation vs. the combined/main set
    # (P_j_LB in the publication)
    # NB: Some tweaks for pseudocounts
    pvalue = binom.pmf(range(int(math.ceil(fg_tot)), fg_size + 2),
                       fg_size + 1, p_j).sum()
    if pvalue == 1.0:
        logging.info("Meaningless p-value: p_j=%s, fg=%s vs. bg=%s",
                        p_j, fg_tot, bg_counts)

    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    "Column probability using the ball-in-urn model."
    col_counts = count_col(col, weights, aa_freqs, pseudo_size)
    col_tot = col_counts['S'] + col_counts['T'] + col_counts['Y']
    p_j = (aa_freqs['S'] + aa_freqs['T'] + aa_freqs['Y'])
    # pvalue = binom.pmf(range(cons_count, aln_size + 1),
                       # aln_size, p_j).sum()
    pvalue = binom.pmf(range(int(math.ceil(col_tot)), aln_size + 2),
                       aln_size + 1, p_j).sum()
    return pvalue



# --- old ---

def compare_aln(fg_aln, bg_aln):
    """Compare alignments using the ball-in-urn model.

    Like CHAIN does.
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
    fg_weights = [1]*fg_size
    pseudocounts = combined_frequencies(fg_aln, fg_weights, bg_aln, bg_weights)
    hits = []
    for faa, baa, fg_col, bg_col in zip(fg_cons, bg_cons, fg_cols, bg_cols):
        if faa == '-' or baa == '-':
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.0
        else:
            # Cumulative binomial test
            # Number of consensus-type residues in the foreground column
            fg_counts = count_col(fg_col, fg_weights, pseudocounts)
            fg_tot = fg_counts['S'] + fg_counts['T'] + fg_counts['Y']
            # Consensus residue frequency in the combined alignment column
            bg_counts = count_col(bg_col, bg_weights, pseudocounts)
            p_j = (bg_counts['S'] + bg_counts['T'] + bg_counts['Y'] + fg_tot
                  ) / (bg_size + fg_size + 2.0) # pseudocount size = 1.0

            # Probability of fg col conservation vs. the combined/main set
            # (P_j_LB in the publication)
            # NB: Some tweaks for pseudocounts
            pvalue = binom.pmf(range(int(math.ceil(fg_tot)), fg_size+2),
                               fg_size+1, p_j).sum()
            if pvalue == 1.0:
                logging.info("Meaningless p-value: p_j=%s, fg=%s vs. bg=%s",
                             p_j, fg_tot, bg_counts)
        hits.append((faa, baa, pvalue))
    return hits

