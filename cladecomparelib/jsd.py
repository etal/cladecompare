"""Jensen-Shannon divergence (JSD).

Also known as information radius (IRad) or total divergence to the average.

JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))

Single mode: let Q be the distribution of equal amino acid frequencies.

See:
    https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
"""

from biofrills import consensus, alnutils

from .shared import (count_col, combined_frequencies, entropy, MAX_ENTROPY, 
                     scale_counts)


def compare_cols(fg_col, fg_cons, fg_size, fg_weights,
                 bg_col, bg_cons, bg_size, bg_weights,
                 aa_freqs, pseudo_size):
    """1 - Jensen-Shannon divergence between the column.

    JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))
    """
    # ENH: choose \pi weights according to size of each set (after weighting)?
    fg_counts = scale_counts(count_col(fg_col, fg_weights))
    bg_counts = scale_counts(count_col(bg_col, bg_weights, aa_freqs, pseudo_size))
    # bg_counts = scale_counts(count_col(bg_col, bg_weights))
    mix = {}
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        mix[aa] = 0.5*fg_counts.get(aa, 0.) + 0.5*bg_counts.get(aa, 0)
    jsd = entropy(mix.values()) - 0.5*(entropy(fg_counts.values()) +
                                        entropy(bg_counts.values()))
    return (1 - jsd)**2


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """1 - Jensen-Shannon divergence of column from equal aa frequencies."""
    # XXX use equal freqs or aa_freqs as background?
    counts = scale_counts(count_col(col, weights))
    mix = {}
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        mix[aa] = 0.5*counts.get(aa, 0.) + 0.5*0.05
    jsd = entropy(mix.values()) - 0.5*(entropy(counts.values()) +
                                        MAX_ENTROPY)
    return (1 - jsd)**2

