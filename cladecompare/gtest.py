"""Compare two alignments via G-test of column amino acid frequencies.

In paired mode, each column gives a contingency table of weighted residue
counts (usually not integers when the sequences are weighted)::

                    Ala Cys Asp ...
    foreground      FA  FC  FD  ...
    background      BA  BC  BD  ...

with one column of the table for each residue type found in either set. The
test is the G-test of homogeneity (Dunning 1993, section 5.3): the null
hypothesis is that the two sets draw the column's residues from one
distribution, whose frequencies are estimated from both sets together. The
expected count of each cell is its row total times its column total divided by
the table total, and::

        G = 2 * sum_ij( O_ij * ln(O_ij / E_ij) )

over the cells with O_ij > 0. With K residue types in the table, G follows the
chi-squared distribution with (K-1)x(2-1) = K-1 degrees of freedom, and the
p-value is its upper tail, scipy.stats.chi2.sf(G, K - 1). A column with fewer
than two types gives 1. No pseudocounts are added, so the p-value of a column
depends on that column alone.

For a single alignment, the observed counts are each column's weighted residue
counts plus pseudocounts of total size 0.5 (core.PSEUDO_SIZE), in proportion to
the amino acid frequencies of the whole alignment, and the expected counts are
those frequencies scaled to the column's observed total, so that the test
compares the column's composition with the alignment's on 19 degrees of
freedom.

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
    """Compare amino acid frequencies between aligned columns via G-test.

    The test covers the whole composition, so the residue returned with the
    p-value, which names the site, is the foreground consensus.
    """
    fg_counts = count_col(fg_col, fg_weights)
    bg_counts = count_col(bg_col, bg_weights)
    # Column totals of the table, over the residue types present in either set
    col_totals = {
        aa: both
        for aa in set(fg_counts) | set(bg_counts)
        if (both := fg_counts.get(aa, 0.0) + bg_counts.get(aa, 0.0))
    }
    if len(col_totals) < 2:
        # The compositions cannot differ, and the test has no degrees of freedom
        return fg_cons, 1.0
    fg_total = math.fsum(fg_counts.values())
    bg_total = math.fsum(bg_counts.values())
    total = fg_total + bg_total
    G = 2 * math.fsum(
        count * math.log(count * total / (row_total * col_totals[aa]))
        for counts, row_total in ((fg_counts, fg_total), (bg_counts, bg_total))
        for aa, count in counts.items()
        if count
    )
    return fg_cons, chi2.sf(G, len(col_totals) - 1)


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Compare column amino acid frequencies to overall via G-test.

    Return the consensus residue, which names the site, and the p-value.
    """
    observed = count_col(col, weights, aa_freqs, pseudo_size)
    # Scale the expected counts to the observed total, pseudocounts included
    total = math.fsum(observed.values())
    G = 2 * sum(
        obsv * math.log(obsv / (total * aa_freqs.get(aa, 0.0)))
        for aa, obsv in observed.items()
    )
    return cons_aa, chi2.sf(G, 19)
