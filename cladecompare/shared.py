import collections

standard_aa = 'ACDEFGHIKLMNPQRSTVWY'

def count_col(col, weights):
    counts = collections.defaultdict(float)
    for char, wt in zip(col, weights):
        if char in '-.X':
            continue
        counts[char] += wt
    return counts
