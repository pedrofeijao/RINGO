#!/usr/bin/env python2
import pyximport; pyximport.install()
import argparse
import os
from dendropy import Tree
import file_ops
import algorithms
import plot_bp
from model import BPGraph, CType, Genome, Chromosome
from pyx import canvas, style, color
from ringo_config import RingoConfig

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solves the Small Phylogeny problem with I.G. InDel")
    # parser.add_argument('-r', '--random', action="store_true", default=False,
    #                     help="Choose random adjacencies in the ambiguous components.")
    # parser.add_argument('-v', '--verbose', action="store_true", default=False)
    # parser.add_argument('-rep', '--repeat', type=int, default=1, help="Number of repeated runs.")
    parser.add_argument("-i", "--input_genomes", required=True, type=str, help="Leaf genomes file.")
    parser.add_argument("-t", "--tree", type=str, required=True, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str, help="Output folder. If not given, output is written to the same location of the genomes file.")

    parser.add_argument("-w", "--adj_weights_file", default="custom_weight", type=str, help="Ancestral adjacency weights file.")

    parser.add_argument("-f", "--weight_filter", type=float, default=0,
                        help="Filter cutoff for adjacency weights, smaller weights are removed.")

    parser.add_argument("-bp", action="store_true", default=False, help="Writes PDFs files each with a plot of the Breakpoint Graph between two siblings that was used to reconstruct the parent node.")
    param = parser.parse_args()

    leaf_genomes = file_ops.open_genome_file(param.input_genomes)
    tree = file_ops.open_newick_tree(param.tree, label_internal_nodes=True)

    # if custom, use my weighting scheme:
    if param.adj_weights_file == "custom_weight":
        internalAdjWeight = algorithms.ancestral_adjacency_weights(Tree(tree), leaf_genomes)
    else:
        # if weights are given, use: (usually DeClone weights):
        internalAdjWeight = file_ops.open_ancestral_weights(param.adj_weights_file, cutoff=param.weight_filter)

    filename = os.path.basename(param.adj_weights_file)

    reconstructed = algorithms.ig_indel_small_phylogeny(leaf_genomes, tree, internalAdjWeight)

    # output:
    cfg = RingoConfig()
    folder = param.output if param.output is not None else os.path.dirname(param.input_genomes)
    out_filename = os.path.join(folder, cfg.ringo_output_genomes())
    tree_out_filename = os.path.join(folder, cfg.ringo_output_tree())
    if not os.path.exists(folder):
        os.mkdir(folder)
    file_ops.write_genomes_to_file(reconstructed, out_filename)
    file_ops.write_newick_tree(tree, tree_out_filename)

    # BP graph plot:
    if param.bp:
        reconstructed.update(leaf_genomes)
        plot_bp.draw_all_bp(reconstructed, tree, folder, internalAdjWeight, ancestral)
