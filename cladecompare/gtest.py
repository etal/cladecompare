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

from scipy.stats import chisqprob

from biofrills import consensus, alnutils

from shared import count_col


def compare_aln(fg_aln, bg_aln):
    fg_weights = alnutils.sequence_weights(fg_aln, 'none')
    fg_size = sum(fg_weights)
    bg_weights = alnutils.sequence_weights(bg_aln, 'none')
    bg_size = sum(bg_weights)
    pseudo_size = 1.0 # math.sqrt(bg_size)
    # Overall aa freqs for pseudocounts
    aa_freqs = alnutils.aa_frequencies(fg_aln._records + bg_aln._records,
                                       gap_chars='-.X',
                                       weights=fg_weights + bg_weights)
    fg_cons = consensus.consensus(fg_aln, weights=fg_weights, trim_ends=False,
                                  gap_threshold=0.8)
    bg_cons = consensus.consensus(bg_aln, weights=bg_weights, trim_ends=False,
                                  gap_threshold=0.8)
    fg_cols = zip(*fg_aln)
    bg_cols = zip(*bg_aln)
    hits = []
    for faa, baa, fg_col, bg_col in zip(fg_cons, bg_cons, fg_cols, bg_cols):
        if faa == '-' or baa == '-':
            # Ignore indel columns -- there are better ways to look at these
            pvalue = 1.
        else:
            # Calculate the "expected" aa frequencies
            bg_counts = count_col(bg_col, bg_weights)
            expected = {}
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                expected[aa] = ((bg_counts[aa] + pseudo_size*aa_freqs[aa])
                                / (bg_size + pseudo_size)
                            ) * fg_size  # Scale to same size as foreground

            # Calculate the G-value of observed vs. expected
            observed = count_col(fg_col, fg_weights)
            G = 2 * sum(obsv * math.log(obsv/expected[aa])
                        for aa, obsv in observed.iteritems())
            # 4. Calculate the Chi-squared p-value of G
            pvalue = chisqprob(G, 19)
        hits.append((faa, baa, pvalue))

    return hits

