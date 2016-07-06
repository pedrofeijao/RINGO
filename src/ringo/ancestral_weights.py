#!/usr/bin/env python2
import pyximport; pyximport.install()
import argparse
from dendropy import Tree
import algorithms
import file_ops

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates an adjacency weight file for internal node adjacencies, given a newick tree and the leaf genomes file.")
    parser.add_argument("-i", "--input_genomes", required=True, type=str, help="Leaf genomes file.")
    parser.add_argument("-t","--tree", type=str, required=True, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file.")
    param = parser.parse_args()

    # read input:
    leaf_genomes = file_ops.open_genome_file(param.input_genomes)
    tree = file_ops.open_newick_tree(param.tree, label_internal_nodes=True)

    # calc weight:
    internalAdjWeight = algorithms.ancestral_adjacency_weights(Tree(tree), leaf_genomes)

    # write output:
    file_ops.write_ancestral_weights(internalAdjWeight, param.output)
