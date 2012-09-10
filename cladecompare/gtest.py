"""Compare two alignments via G-test of column amino acid frequencies.

Column composition:
    - do the whole vector = characters in the column?
        - null = main set, alt = foreground, 

                Ala Cys Asp ...
    observed    OA  OC  OD  ...
    expected    EA  EC  ED  ...

    G-test for independence/fit:
        * need pseudocounts in BG to avoid errors
        * loses info abt. set size
        * gets that back from df & individual observations?

        G = 2 * sum_ij( O_ij * ln(Oij/Eij) )
        sum is taken over all non-empty cells
        -> chisq value

        df = (N-k)
"""

import math

from scipy.stats import chisqprob

from biofrills import consensus, alnutils

from shared import count_col, standard_aa


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
    fg_cons = consensus.consensus(fg_aln, weights=fg_weights)
    bg_cons = consensus.consensus(bg_aln, weights=bg_weights)
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
            for aa in standard_aa:
                expected[aa] = ((bg_counts[aa] + pseudo_size*aa_freqs[aa])
                                / (bg_size + pseudo_size)
                            ) * fg_size  # Scale to same size as foreground

            # Calculate the G-value of observed vs. expected
            observed = count_col(fg_col, fg_weights)
            G = 2 * sum(obsv * math.log(obsv/expected[aa])
                        for aa, obsv in observed.iteritems())
            # 4. Calculate the Chi-squared p-value of G
            pvalue = chisqprob(G, bg_size)  # df = N - k
        hits.append((faa, baa, pvalue))

    return hits

