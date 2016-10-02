#!/usr/bin/env python2
import argparse
import collections
import os
import pyximport;
import sys

pyximport.install()
from model import BPGraph, CType
import file_ops
import numpy as np
from operator import mul
import itertools
from decimal import Decimal
import random
import matplotlib
matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt


def expected_dcj_distance(genome1, genome2, n=0):
    a = len(genome1.adjacency_set())
    a2 = len(genome2.adjacency_set())
    BP = a - len(genome1.common_adjacencies(genome2))
    g = len(genome1.gene_set())
    # n = np.math.sqrt(g)
    n = g
    if BP == n:
        BP = n-1
    # import ipdb;ipdb.set_trace()
    return np.math.log(1.0 - (BP * (2.0*n - 1)) / (a * (2.0*n - 2))) / np.math.log(1 - (1.0 / (n - 1)) - 1.0 / n)


def probability(n, cycle_dist, st):
    return Decimal(n_scenarios(cycle_dist, st)) / (n * (n - 1)) ** (st)

def cycle_splits(size):
    # total: s(s-1)/2
    # s of (1,s-1)
    # s of (2,s-2)
    # s of (3,s-3)
    # ...
    # and s/2 of (s/2,s/2) if s even
    for i in range(1, (size - 1) / 2 + 1):  # python 2: size/2 is rounded down, that is like a floor;
        yield size, (i, size - i)
    if size % 2 == 0:
        yield size / 2, (size / 2, size / 2)


def memoize(f):
    cache = {}
    return lambda *args: cache[args] if args in cache else cache.update({args: f(*args)}) or cache[args]


@memoize
def n_scenarios(cycle_dist, steps):
    n = sum(cycle_dist)
    c = len(cycle_dist)
    dist = n - c
    if steps < dist:
        return 0
    # d+1 I know:
    elif steps == dist + 1:
        l = [(cycle - 1) for cycle in cycle_dist if cycle > 1]
        m = reduce(mul, [(l_i + 1) ** (l_i - 1) for l_i in l], 1)
        s = 0
        for l_p in l:
            f = np.math.factorial(l_p)
            s1 = sum([f*((l_p+1)**i)/np.math.factorial(i) for i in range(l_p)])
            s1 *= m
            s1 /= (l_p+1)**(l_p-1)
            s += s1
        p1 = np.math.factorial(dist + 1)*s/2
        p1 /= reduce(mul, [np.math.factorial(l_i) for l_i in l], 1)
        return p1
    # d is simple:
    elif steps == dist:
        l = [(cycle - 1) for cycle in cycle_dist if cycle > 1]
        p1 = np.math.factorial(dist) / reduce(mul, [np.math.factorial(l_i) for l_i in l], 1)
        p2 = reduce(mul, [(l_i + 1) ** (l_i - 1) for l_i in l], 1)
        return p1 * p2
    else:  # more steps than distance; recursive:
        # generate all possible cycle distributions from the current one:
        cycle_dist_count = collections.defaultdict(lambda: 0)

        cycle_dist_l = list(cycle_dist)
        # find all cycle splits:
        for idx, size_i in enumerate(cycle_dist_l):
            for qty, (size_1, size_2) in cycle_splits(size_i):
                new_dist = tuple(sorted(cycle_dist_l[:idx] + [size_1, size_2] + cycle_dist_l[(idx + 1):]))
                cycle_dist_count[new_dist] += qty

        # cycle freezes:
        # freezes: C(s_i,2) for each cycle;

        n_freezes = sum([l_i * (l_i - 1) / 2 for l_i in cycle_dist])
        cycle_dist_count[cycle_dist] += n_freezes

        # cycle merges:
        # s_i x s_j of (s_i+s_j) for each pair
        for i, j in itertools.combinations(range(len(cycle_dist)), 2):
            l_i, l_j = cycle_dist[i], cycle_dist[j]
            new_dist = tuple(sorted(cycle_dist_l[:i] + cycle_dist_l[(i + 1):j] + [l_i + l_j] + cycle_dist_l[(j + 1):]))
            cycle_dist_count[new_dist] += 2 * l_i * l_j

        # print cycle_dist_count
        return sum(
            [count_i * n_scenarios(cycle_dist_i, steps - 1) for cycle_dist_i, count_i in cycle_dist_count.iteritems()])


def random_walk(g1, g2, steps, n_walks=100000):
    adj_2 = sorted(g2.adjacency_set())
    hit = 0
    for i in range(n_walks):
        adj_1 = [[a, b] for a, b in g1.adjacency_set()]
        for j in range(steps):
            p, q = random.sample(range(len(adj_1)), 2)
            if p < q:
                adj_1[p][0], adj_1[q][0] = adj_1[q][0], adj_1[p][0]
            else:
                adj_1[p][0], adj_1[q][1] = adj_1[q][1], adj_1[p][0]
        adj_1 = sorted([tuple(sorted((a, b))) for a, b in adj_1])
        if adj_1 == adj_2:
            hit += 1
    print "hits: %e" % (float(hit) / n_walks)


# def sort_cycle_with_one_freeze(n):
#     return sum([reduce(mul, [x for x in range(k+1, n+1)]) * (n ** k)  for k in range(n - 1)])/2

def sort_cycle_with_one_freeze(n):
    f = np.math.factorial(n)
    return sum([f * n ** k / np.math.factorial(k) for k in range(n - 1)])/2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Finds the ML estimate for the DCJ distance between 2 genomes.")
    parser.add_argument("-g", type=str, nargs='+', help="Genomes file(s)")
    parser.add_argument("-i", type=int, nargs=2, default=[0, 1], help="Idx of the genomes")
    param = parser.parse_args()


    # n = 100
    # print n_scenarios((n,), n)
    # print sort_cycle_with_one_freeze(n)
    # print sort_cycle_with_one_freeze2(n)
    # print ",".join(map(str, [(i,n_scenarios((i,), i)) for i in range(30, 31)]))
    # print ",".join(map(str, [ (n ,sort_cycle_with_one_freeze(n)) for n in range(30, 31)]))
    # sys.exit()
    n1, n2 = param.i
    for filename in param.g:
        genomes = file_ops.open_genome_file(filename, as_list=True)
        g1 = genomes[int(n1)]
        g2 = genomes[int(n2)]

        bp = BPGraph(g1, g2)
        n = len(bp.common_AB)
        c = len(bp.type_dict[CType.CYCLE])
        cycle_distribution = tuple(sorted([len(x) / 2 for x in bp.type_dict[CType.CYCLE]]))
        # cycle_distribution = tuple([len(x) / 2 for x in bp.type_dict[CType.CYCLE]])
        d = n - c
        x = []
        y = []
        last_step = 0
        down = 0
        max_p = 0
        max_k = 0
        # DCJ estimate:
        est_DCJ = expected_dcj_distance(g1,g2)
        print "Distance:%d" % d,
        print " Estimate: %.1f" % est_DCJ

        # if there is no common adjacency, estimate goes to infinity, also in the DCJ estimate;
        if all([c > 1 for c in cycle_distribution]):
            # cheat and build a new one, by randomly picking an element and then removing a cycle from it;
            cycle_distribution = list(cycle_distribution)
            random.shuffle(cycle_distribution)
            cycle = cycle_distribution.pop()
            cycle_distribution = tuple(sorted([1, cycle - 1] + cycle_distribution))
        for i in range(2*n):
            prob = probability(n, cycle_distribution, d + i)
            print >> sys.stderr, "Steps:%d P:%e" % (d + i, prob)
            x.append(d + i)
            y.append(prob)
            if prob < last_step:
                down += 1
                if down == 2:
                    break
            else:
                down = 0
            if max_p < prob:
                max_p = prob
                max_k = i + d
            last_step = prob

        plt.plot(x, y, 'o-')
        plt.savefig(os.path.join(os.path.dirname(filename), 'ml.pdf'), bbox_inches='tight')
        print "Max:", max_k
        # save results:
        with open(filename+".ml", "w") as f:
            print >> f, "DCJ\tML\tEDCJ"
            print >> f, "%d\t%d\t%.1f" % (d, max_k, est_DCJ)
