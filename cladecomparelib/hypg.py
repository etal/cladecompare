"""Hypergeometric Distribution

    Gendankenexperiment:
    
        Foreground and background sequence sets are pre-defined.

        Given N foreground sequences and M-N background sequences,
    we randomly select N sequences from M.  We consider the consensus
    residue in the foreground as being type I and ask what is the probability
    of observing at least as many type I sequences in our selection as we see
    in the foreground. 

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
    p_j = int(ceil(aa_freqs[cons_aa]*aln_size))
    size_i = int(ceil(aln_size))
    pvalue = float(cons_count_i)/len(col)
    #pvalue = hypergeom.cdf(cons_count_i-1,size_i,
                        #max(cons_count_i,p_j), len(col))
    return pvalue

