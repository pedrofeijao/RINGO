#!/usr/bin/env python2
import ringo_config
cfg = ringo_config.RingoConfig()
import pyximport; pyximport.install(build_dir=cfg.pyximport_build())
import argparse
import os
from dendropy import Tree
import file_ops
import algorithms
import plot_bp
from model import BPGraph, CType, Genome, Chromosome
from pyx import canvas, style, color
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solves the Small Phylogeny problem with I.G. InDel")
    # parser.add_argument('-r', '--random', action="store_true", default=False,
    #                     help="Choose random adjacencies in the ambiguous components.")
    # parser.add_argument('-v', '--verbose', action="store_true", default=False)
    # parser.add_argument('-rep', '--repeat', type=int, default=1, help="Number of repeated runs.")
    parser.add_argument("-i", "--input_genomes", required=True, type=str, help="Leaf genomes file.")
    parser.add_argument("-t", "--tree", type=str, required=True, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str, help="Output folder. If not given, output is written to the same location of the genomes file.")

    parser.add_argument("-w", "--adj_weights_file", type=str, help="Ancestral adjacency weights file.")

    parser.add_argument("-f", "--weight_filter", type=float, default=0,
                        help="Filter cutoff for adjacency weights, smaller weights are removed.")
    parser.add_argument("-p", "--perfect", action="store_true", default=False, help="Force perfect matching for the maximum weight matching of open components.")
    parser.add_argument("-bl", "--estimate_lenghts", type=str, choices=["lp", "least_squares"], help="Estimate tree branch lenghts, greatly improves RINGO results if the input tree does not have branch lengths. Choose one of the available methods.")
    parser.add_argument("-bp", action="store_true", default=False, help="Writes PDFs files each with a plot of the Breakpoint Graph between two siblings that was used to reconstruct the parent node.")
    parser.add_argument("-a", "--ancestral", type=str, help="Ancestral genomes file for a simulation. Can be used to draw the ancestors in the BP plot.")

    parser.add_argument("-r", "--random_repeat", type=int, default=0, help="Number of repeats of randomly fill adjacencies, reNumber of randomly repeats of  ")

    # DEBUG:
    parser.add_argument("--add_open_2_cycles", action="store_true", default=False, help="adds adjacencies from 2-cycles that are AA- or BB-open directly.")

    # PARSE params:
    param = parser.parse_args()


    # test if blossom5 is needed:
    if param.perfect:
        if not file_ops.blossom5_is_available():
            print >> sys.stderr, "Blossom5 not found, it is needed for perfect matching. Either install blossom5 or do not use the '-p' perfect matching option."
            sys.exit(-1)

    extant_genomes = file_ops.open_genome_file(param.input_genomes)
    tree = file_ops.open_newick_tree(param.tree, label_internal_nodes=True)

    # estimate lenghts:
    if param.estimate_lenghts is not None:
        if param.estimate_lenghts == "lp":
            algorithms.estimate_branch_lengths_lp(tree, extant_genomes)
        elif param.estimate_lenghts == "least_squares":
            algorithms.estimate_branch_lengths_least_squares(tree, extant_genomes)
    else:
        # default is to run it anyway, if the tree has no branch lengths:
        if any([node.edge.length is None for node in tree.leaf_nodes()]):
            algorithms.estimate_branch_lengths_lp(tree, extant_genomes)

    # if no weights file is given, use default weighting scheme:
    if param.adj_weights_file is None:
        internalAdjWeight = algorithms.ancestral_adjacency_weights(Tree(tree), extant_genomes)


    else:
        # if weights are given, use: (usually DeClone weights):
        internalAdjWeight = file_ops.open_ancestral_weights(param.adj_weights_file, cutoff=param.weight_filter)

    reconstructed = algorithms.ig_indel_small_phylogeny(extant_genomes, tree, internalAdjWeight,
                        perfect_matching=param.perfect, random_repeat=param.random_repeat, add_open_2_cycles=param.add_open_2_cycles)

    # output:
    if param.output is None:
        param.output = os.getcwd()

    folder = param.output
    out_filename = os.path.join(folder, cfg.ringo_output_genomes())
    tree_out_filename = os.path.join(folder, cfg.ringo_output_tree())
    if folder != "" and not os.path.exists(folder):
        os.mkdir(folder)
    file_ops.write_genomes_to_file(reconstructed, out_filename)
    file_ops.write_newick_tree(tree, tree_out_filename)
    # Save parameters:
    file_ops.write_ringo_parameters(param, folder)

    # ancestral genomes:
    ancestral = file_ops.open_genome_file(param.ancestral) if param.ancestral is not None else None
    # BP graph plot:
    if param.bp:
        reconstructed.update(extant_genomes)
        plot_bp.draw_all_bp(reconstructed, tree, folder, internalAdjWeight, ancestral)
