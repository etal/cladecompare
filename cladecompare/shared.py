import collections
from math import log

MAX_ENTROPY = log(20, 2)

def count_col(col, weights):
    counts = collections.defaultdict(float)
    for char, wt in zip(col, weights):
        if char in '-.X':
            continue
        counts[char] += wt
    return counts


def entropy(counts):
    scale = 1.0 / sum(counts)
    probs = [c * scale for c in counts if c]
    # try:
    return -sum(p * log(p, 2) for p in probs)
    # except ValueError:
    #     print locals()
    #     return 0.
