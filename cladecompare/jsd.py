"""Jensen-Shannon divergence (JSD).

Also known as information radius (IRad) or total divergence to the average.

JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))

Single mode: let Q be the distribution of equal amino acid frequencies.

The value is JSD itself, in bits, larger for more divergent columns. It lies
between 0 and 1 for columns of the 20 amino acids, X and gaps. Any other
letter, such as B or Z, is counted by count_col but absent from the
mixture of the 20, so the value of a column holding one, and in paired mode
of every column of a pair holding one, through the background's
pseudocounts, is not a divergence: it can fall below 0, exceed 1, or lie
within [0, 1] and still be wrong. It is a score, not a p-value; VALUE_LABEL
tells core, which labels the '.out' lines accordingly.

See:
    https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
"""

from .shared import MAX_ENTROPY, count_col, entropy, scale_counts

VALUE_LABEL = "score"

# Two entropies of a few bits each can differ by a rounding error, which
# leaves the JSD of nearly equal compositions up to a few 1e-15 below 0. A
# letter outside the 20 shifts the JSD by far more: less than 1e-12 only for
# a frequency of that letter under about 5e-14, which no real alignment has.
ROUNDING = 1e-12


def nonnegative(jsd):
    """Take a value less than ROUNDING below 0 as 0.

    Any other value is returned as it is: a larger negative, such as a letter
    outside the 20 gives, is not a rounding error, and a nan stays nan.
    """
    return 0.0 if -ROUNDING < jsd < 0.0 else jsd


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
    """Jensen-Shannon divergence between the columns.

    JSD(P,Q) = H(0.5*P + 0.5*Q) - (0.5*H(P) + 0.5*H(Q))

    Return the foreground consensus, which names the site, and the score.
    """
    # ENH: choose \pi weights according to size of each set (after weighting)?
    fg_counts = scale_counts(count_col(fg_col, fg_weights))
    bg_counts = scale_counts(count_col(bg_col, bg_weights, aa_freqs, pseudo_size))
    mix = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        mix[aa] = 0.5 * fg_counts.get(aa, 0.0) + 0.5 * bg_counts.get(aa, 0)
    jsd = entropy(mix.values()) - 0.5 * (
        entropy(fg_counts.values()) + entropy(bg_counts.values())
    )
    return fg_cons, nonnegative(jsd)


def compare_one(col, cons_aa, aln_size, weights, aa_freqs, pseudo_size):
    """Jensen-Shannon divergence of column from equal aa frequencies.

    Return the consensus residue, which names the site, and the score.
    """
    # XXX use equal freqs or aa_freqs as background?
    counts = scale_counts(count_col(col, weights))
    mix = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        mix[aa] = 0.5 * counts.get(aa, 0.0) + 0.5 * 0.05
    jsd = entropy(mix.values()) - 0.5 * (entropy(counts.values()) + MAX_ENTROPY)
    return cons_aa, nonnegative(jsd)
