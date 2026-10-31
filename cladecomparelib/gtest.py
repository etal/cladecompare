"""Compare two alignments via G-test of column amino acid frequencies.

The contingency table looks like::

                Ala Cys Asp ...
    observed    OA  OC  OD  ...
    expected    EA  EC  ED  ...

where O<resname> is the "observed" residue count in the foreground (using
weighted sequences, so usually not an integer) and E<resname> is the "expected"
residue count, calculated from the background column residue frequencies and
scaled to match the size of the foreground. Pseudocounts of total size 0.5
(core.PSEUDO_SIZE), in proportion to the amino acid frequencies of the whole
alignment, are added to the background counts so that no expected count is
zero.

This is a goodness-of-fit test where the null hypothesis is that the foreground
column residues were drawn from the distribution represented by the
background set. The test assumes the foreground is a sub-clade within the
background clade.

Formula::

        G = 2 * sum_i( O_i * ln(Oi/Ei) )

The test statistic G follows the chi-squared distribution with 19 degrees of
freedom (20 amino acids, 2 compared sets -> (20-1)x(2-1) = 19), and the p-value
is its upper tail, scipy.stats.chi2.sf(G, 19).

For a single alignment, the observed counts are each column's weighted residue
counts plus the same pseudocounts, and the expected counts are the amino acid
frequencies of the whole alignment scaled to the column's observed total, so
that the test compares the column's composition with the alignment's on 19
degrees of freedom.

"""

import math

from scipy.stats import chi2

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
    """Compare amino acid frequencies between aligned columns via G-test."""
    # Calculate the "expected" aa frequencies
    bg_counts = count_col(bg_col, bg_weights, aa_freqs, pseudo_size)
    expected = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        # Scale to same size as foreground
        expected[aa] = fg_size * (bg_counts[aa] / (bg_size + pseudo_size))
    # Calculate the G-value of observed vs. expected
    observed = count_col(fg_col, fg_weights)
    G = 2 * sum(obsv * math.log(obsv / expected[aa]) for aa, obsv in observed.items())
    # 4. Calculate the Chi-squared p-value of G
    pvalue = chi2.sf(G, 19)
    return pvalue


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Compare column amino acid frequencies to overall via G-test."""
    observed = count_col(col, weights, aa_freqs, pseudo_size)
    # Scale the expected counts to the observed total, pseudocounts included
    total = math.fsum(observed.values())
    G = 2 * sum(
        obsv * math.log(obsv / (total * aa_freqs.get(aa, 0.0)))
        for aa, obsv in observed.items()
    )
    pvalue = chi2.sf(G, 19)
    return pvalue
