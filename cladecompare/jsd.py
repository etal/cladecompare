"""Jensen-Shannon divergence (JSD).

Also known as information radius (IRad) or total divergence to the average.

JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))

Single mode: let Q be the distribution of equal amino acid frequencies.

See:
    https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
"""
# XXX choose \pi weights according to size of each set (after weighting)?

from biofrills import consensus, alnutils

from shared import count_col, entropy, MAX_ENTROPY


def compare_one(aln):
    """Jensen-Shannon divergence within each column of the alignment.

    JSD(P) = H(0.5*P + 0.5*[0.05...]) - (0.5*H(P) + 0.5*H_max)
    """
    weights = alnutils.sequence_weights(aln, 'sum1')
    cons = consensus.consensus(aln, weights=weights, trim_ends=False,
                               gap_threshold=0.8)
    cols = zip(*aln)
    hits = []
    for aa, col in zip(cons, cols):
        if aa == '-':
            # Ignore indel columns -- there are better ways to look at these
            jsd = 1.
        else:
            mix = {}
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                mix[aa] = 0.5*counts.get(aa, 0.) + 0.5*0.05
            jsd = entropy(mix.values()) - 0.5*(entropy(fg_counts.values()) +
                                               MAX_ENTROPY)
        hits.append((aa, jsd))

    return hits


def compare_aln(fg_aln, bg_aln):
    """Jensen-Shannon divergence between each column of the alignments.

    JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))
    """
    fg_weights = alnutils.sequence_weights(fg_aln, 'sum1')
    bg_weights = alnutils.sequence_weights(bg_aln, 'sum1')
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
            jsd = 1.
        else:
            bg_counts = count_col(bg_col, bg_weights)
            fg_counts = count_col(fg_col, fg_weights)
            mix = {}
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                mix[aa] = 0.5*fg_counts.get(aa, 0.) + 0.5*bg_counts.get(aa, 0)
            jsd = entropy(mix.values()) - 0.5*(entropy(fg_counts.values()) +
                                               entropy(bg_counts.values()))
        hits.append((faa, baa, 1 - jsd))

    return hits


