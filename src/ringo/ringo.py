#!/usr/bin/env python2
import ringo_config

cfg = ringo_config.RingoConfig()
import pyximport;

pyximport.install(build_dir=cfg.pyximport_build())
import argparse
import os
from dendropy import Tree
import file_ops
import algorithms
import plot_bp
import sys
import numpy as np
import dcj


def find_min_q(q, d):
    min_pair = min_val = (0, 0)
    for i in range(1, len(q)):
        for j in range(i):
            if (q[i][j], d[i][j]) < min_val:
                min_val = (q[i][j], d[i][j])
                min_pair = i, j
    return min_pair


def get_dist(d_matrix, i, j):
    if i < j:
        i, j = j, i
    return d_matrix[i][j]


def calc_q_matrix(d_matrix):
    q = np.zeros((len(d_matrix), len(d_matrix)), int)
    # vector of row sums:
    r = [sum([d_matrix[i][k] for k in range(i)]) + sum([d_matrix[k][i] for k in range(i + 1, len(d_matrix))]) for i in
         range(len(d_matrix))]
    # fill lower diagonal of q matrix:
    for i in range(1, len(d_matrix)):
        for j in range(i):
            q[i][j] = (len(d_matrix) - 2) * d_matrix[i][j] - r[i] - r[j]
    return r, q


def resolve_tree(tree):
    changed = False
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            continue
        else:
            children = node.child_nodes()
            # resolve to binary if needed:
            if len(children) > 2:
                changed = True
                # NJ for these children:
                taxa = [n for n in node.child_nodes()]

                d_matrix = np.zeros((len(taxa), len(taxa)), float)
                for i in range(len(taxa) - 1):
                    for j in range(i+1, len(taxa)):
                        d_matrix[j][i] = taxa[i].edge.length + taxa[j].edge.length

                while True:
                    r, q = calc_q_matrix(d_matrix)
                    # find pair to join:
                    idx_a, idx_b = find_min_q(q, d_matrix)
                    if idx_a > idx_b:
                        idx_a, idx_b = idx_b, idx_a
                    # if >=3 children, create new node:
                    if len(taxa) >= 3:
                        # TODO: review naming:
                        new_node = node.new_child()
                        new_node.label = "%s_%s" % (taxa[idx_a].label[:4], taxa[idx_b].label[:4])
                        # join the neighbours in the new node:
                        taxa[idx_a].parent_node = new_node
                        taxa[idx_b].parent_node = new_node

                    # adjust edge lengths
                        edge_a_new = d_matrix[idx_b][idx_a] * 0.5 + (r[idx_a] - r[idx_b]) / (
                            2 * len(d_matrix) - 4)
                        edge_b_new = d_matrix[idx_b][idx_a] - edge_a_new
                    else:
                        edge_a_new = edge_b_new = d_matrix[idx_b][idx_a]/2
                    taxa[idx_a].edge.length = edge_a_new
                    taxa[idx_b].edge.length = edge_b_new

                    # test end:
                    if len(taxa) == 2:

                        break

                    # redo matrices:
                    new_taxa = [new_node] + taxa[:idx_a] + taxa[(idx_a + 1):idx_b] + taxa[(idx_b + 1):]
                    new_d = np.delete(np.delete(d_matrix, (idx_a, idx_b), 0), (idx_a, idx_b), 1)
                    new_d = np.r_[[np.zeros(len(new_taxa))], np.c_[np.zeros(len(new_taxa) - 1), new_d]]
                    new_idx = 1
                    for i in range(0, len(taxa)):
                        if i == idx_a or i == idx_b:
                            continue
                        new_d[new_idx][0] = 0.5 * (
                        get_dist(d_matrix, idx_a, i) + get_dist(d_matrix, idx_b, i) - d_matrix[idx_b][idx_a])
                        new_idx += 1
                    # update:
                    taxa = new_taxa
                    d_matrix = new_d

    return changed

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solves the Small Phylogeny problem with I.G. InDel")
    # parser.add_argument('-r', '--random', action="store_true", default=False,
    #                     help="Choose random adjacencies in the ambiguous components.")
    # parser.add_argument('-v', '--verbose', action="store_true", default=False)
    # parser.add_argument('-rep', '--repeat', type=int, default=1, help="Number of repeated runs.")
    parser.add_argument("-i", "--input_genomes", required=True, type=str, help="Leaf genomes file.")
    parser.add_argument("-t", "--tree", type=str, required=True, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str,
                        help="Output folder. If not given, output is written to the same location of the genomes file.")

    parser.add_argument("-w", "--adj_weights_file", type=str, help="Ancestral adjacency weights file.")

    parser.add_argument("-f", "--weight_filter", type=float, default=0,
                        help="Filter cutoff for adjacency weights, smaller weights are removed.")
    parser.add_argument("-p", "--perfect", action="store_true", default=False,
                        help="Force perfect matching for the maximum weight matching of open components.")
    parser.add_argument("-bl", "--estimate_lenghts", type=str, choices=["lp", "least_squares"],
                        help="Estimate tree branch lenghts, greatly improves RINGO results if the input tree does not have branch lengths. Choose one of the available methods.")
    parser.add_argument("-bp", action="store_true", default=False,
                        help="Writes PDFs files each with a plot of the Breakpoint Graph between two siblings that was used to reconstruct the parent node.")
    parser.add_argument("-a", "--ancestral", type=str,
                        help="Ancestral genomes file for a simulation. Can be used to draw the ancestors in the BP plot.")

    parser.add_argument("-r", "--random_repeat", type=int, default=0,
                        help="Number of repeats of randomly fill adjacencies, reNumber of randomly repeats of  ")

    # DEBUG:
    parser.add_argument("--add_open_2_cycles", action="store_true", default=False,
                        help="adds adjacencies from 2-cycles that are AA- or BB-open directly.")

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

    # Fix multibranching trees:
    changed = resolve_tree(tree)
    # if the tree changed, it might be a good idea to rerun the branch length detection?
    if changed and param.estimate_lenghts is not None:
        if param.estimate_lenghts == "lp":
            algorithms.estimate_branch_lengths_lp(tree, extant_genomes)
        elif param.estimate_lenghts == "least_squares":
            algorithms.estimate_branch_lengths_least_squares(tree, extant_genomes)

    # if no weights file is given, use default weighting scheme:
    if param.adj_weights_file is None:
        internalAdjWeight = algorithms.ancestral_adjacency_weights(Tree(tree), extant_genomes)
    else:
        # if weights are given, use: (usually DeClone weights):
        internalAdjWeight = file_ops.open_ancestral_weights(param.adj_weights_file, cutoff=param.weight_filter)

    # main alg:
    reconstructed = algorithms.ig_indel_small_phylogeny(extant_genomes, tree, internalAdjWeight,
                                                        perfect_matching=param.perfect,
                                                        random_repeat=param.random_repeat,
                                                        add_open_2_cycles=param.add_open_2_cycles)

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
