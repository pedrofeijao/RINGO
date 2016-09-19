#!/usr/bin/env python2
import argparse
import os

import pyximport;
pyximport.install()
from model import Genome
import file_ops
import random


def random_walk(n_genes, n_chromosomes, n_steps):
    id = Genome.identity(n_genes, n_chromosomes, circular=True, name="G_0")
    adjacencies = list([ [a,b] for a,b in id.adjacency_set()])
    for j in range(n_steps):
        p, q = random.sample(range(len(adjacencies)), 2)
        if p < q:
            adjacencies[p][0], adjacencies[q][0] = adjacencies[q][0], adjacencies[p][0]
        else:
            adjacencies[p][0], adjacencies[q][1] = adjacencies[q][1], adjacencies[p][0]
    return [id, Genome.from_adjacency_list("G_%d" % n_steps, adjacencies)]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Does a DCJ random walk.")
    parser.add_argument("-n", type=int, help="Number of genes.")
    parser.add_argument("-c", type=int, help="Number of chromosomes.")
    parser.add_argument("-s", type=int, help="Number of steps.")
    parser.add_argument("-o", type=str, help="Output folder.")
    param = parser.parse_args()

    genomes = random_walk(param.n, param.c, param.s)
    if not os.path.exists(param.o):
        os.mkdir(param.o)
    file_ops.write_genomes_to_file(genomes, os.path.join(param.o, "genomes.txt"))



