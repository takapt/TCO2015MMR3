#! /usr/bin/env python

import sys

def get_scores(filename):
    scores = {}
    for line in open(filename):
        s = line.split()
        seed = int(s[0])
        score = int(s[1])
        scores[seed] = score
    return scores

a, b = sys.argv[1], sys.argv[2]
p, q = get_scores(a), get_scores(b)
seeds = list(set(p) & set(q))
sum_ratio = 0
for seed in seeds:
    x, y = p[seed], q[seed]
    ratio = y / x
    sum_ratio += ratio

    print('{:>5} {:>8} {:>8} {:>7.3f}'.format(seed, x, y, ratio))

total_ratio = sum_ratio / len(seeds)
print('total_ratio: {}'.format(total_ratio))

