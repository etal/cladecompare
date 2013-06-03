"""Compare two alignments via G-test of column amino acid frequencies.

The contingency table looks like::

                Ala Cys Asp ...
    observed    OA  OC  OD  ...
    expected    EA  EC  ED  ...

where O<resname> is the "observed" residue count in the foreground (using
weighted sequences, so usually not an integer) and E<resname> is the "expected"
residue count, calculated from the background column residue frequencies and
scaled to match the size of the foreground. A pseudocount size of 1 is used in
the background frequencies to avoid divide-by-zero issues and, I suppose,
account for possible sampling and alignment errors/issues.

This is a goodness-of-fit test where the null hypothesis is that the foreground
column residues were drawn from the the distribution represented by the
background set. The test assumes the foreground is a sub-clade within the
background clade.

Formula::

        G = 2 * sum_i( O_i * ln(Oi/Ei) )

The test statistic G follows the chi-squared distribution with 19 degrees of
freedom (20 amino acids, 2 compared sets -> (20-1)x(2-1) = 19).

"""

import math

from biofrills import consensus, alnutils
from biofrills.stats.chisq import chisqprob

from .shared import count_col, combined_frequencies


def compare_cols(fg_col, fg_cons, fg_size, fg_weights,
                 bg_col, bg_cons, bg_size, bg_weights,
                 aa_freqs, pseudo_size):
    """Compare amino acid frequencies between aligned columns via G-test."""
    # Calculate the "expected" aa frequencies
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    expected = {}
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
            # Scale to same size as foreground
        expected[aa] = fg_size * (bg_counts[aa] / (bg_size + pseudo_size)) 
    # Calculate the G-value of observed vs. expected
    observed = count_col(fg_col, fg_weights)
    G = 2 * sum(obsv * math.log(obsv/expected[aa])
                for aa, obsv in observed.iteritems())
    # 4. Calculate the Chi-squared p-value of G
    pvalue = chisqprob(G, 19)
    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Compare column amino acid frequencies to overall via G-test."""
    observed = count_col(col, weights, aa_freqs, pseudo_size)
    G = 2 * sum(obsv * math.log(obsv / aa_freqs.get(aa, 0.0))
                for aa, obsv in observed.iteritems())
    pvalue = chisqprob(G, 19)
    return pvalue

