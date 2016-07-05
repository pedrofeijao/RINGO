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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solves the Small Phylogeny problem with I.G. InDel")
    # parser.add_argument('-r', '--random', action="store_true", default=False,
    #                     help="Choose random adjacencies in the components without guide.")
    # parser.add_argument('-v', '--verbose', action="store_true", default=False)
    parser.add_argument('-rep', '--repeat', type=int, default=1, help="Number of repeated runs.")
    parser.add_argument("-i", "--input_genomes", type=str, help="Leaf genomes file.")
    parser.add_argument("-tree", type=str, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output folder.")

    parser.add_argument("-w", "--adj_weights_file", default="custom", type=str, help="internal weights file.")

    parser.add_argument("-f", "--weight_filter", type=float, default=0,
                        help="Filter cutoff for adjacency weights, smaller weights are removed.")


    # DEBUG:
    parser.add_argument("-a", "--ancestral", type=str, help="Ancestral genomes file, just to test.")
    param = parser.parse_args()

    leaf_genomes = file_ops.open_genome_file(param.input_genomes)
    tree = file_ops.open_newick_tree(param.tree, label_internal_nodes=True)

    # if custom, use my weighting scheme:
    if param.adj_weights_file == "custom":
        internalAdjWeight = algorithms.ancestral_adjacency_weights(Tree(tree), leaf_genomes)
    else:
        # if weights are given, use: (usually DeClone weights):
        internalAdjWeight = file_ops.open_ancestral_weights(param.adj_weights_file, cutoff=param.weight_filter)

    folder = param.output
    filename = os.path.basename(param.adj_weights_file) if param.adj_weights_file is not None else "custom"
    out_filename = os.path.join(folder, "ig_indel_genomes_%s.txt" % filename)
    ancestral = None
    if param.ancestral is not None:
        ancestral = file_ops.open_genome_file(param.ancestral)

    reconstructed = algorithms.ig_indel_small_phylogeny(leaf_genomes, tree, internalAdjWeight)

    # output:
    if not os.path.exists(folder):
        os.mkdir(folder)
    file_ops.write_genomes_to_file(reconstructed, out_filename)

    # quick test:
    if param.ancestral is not None:
        ancestral = file_ops.open_genome_file(param.ancestral)
        for label, g in ancestral.iteritems():
            if label not in reconstructed:
                continue
            r = reconstructed[label]
            print label
            print "ADJ:", len(g.adjacency_set()), len(r.adjacency_set())
            print "COMMON:", len(g.common_adjacencies(r))
            print

    # BP graph:
    # BP plot:
    plot_ancestral_bp = False
    if plot_ancestral_bp:
        for node in tree.preorder_internal_node_iter():
            if node == tree.seed_node:
                continue
            g1, g2 = reconstructed[node.child_nodes()[0].label], reconstructed[node.child_nodes()[1].label]
            bp = BPGraph(g1, g2)
            # comp_list = bp.type_dict.values()
            comp_list = [{'c': val, 'type': key} for key, c_list in bp.type_dict.iteritems() for val in c_list]
            c = canvas.canvas()
            pos = plot_bp.draw_bp_graph(c, bp, g1, g2)
            # ancestor:
            plot_bp.draw_genome(c, ancestral[node.label].adjacency_set(), pos,
                                [style.linewidth(0.2), color.cmyk.Orange, style.linecap.round],
                                draw_only_in_given_components=comp_list)
            # reconstructed:
            plot_bp.draw_genome(c, reconstructed[node.label].adjacency_set(), pos,
                                [style.linewidth(0.1), color.cmyk.Red, style.linecap.round],
                                draw_only_in_given_components=comp_list)
            # weights:
            for max_w in [x / 10.0 for x in range(1, 11)]:
                plot_bp.draw_genome(c, [adj for adj, w in internalAdjWeight[node.label].iteritems() if w > max_w], pos,
                                    [style.linewidth(max_w / 8), color.cmyk.Sepia],
                                    draw_only_in_given_components=comp_list)

            c.writePDFfile(os.path.join(folder, "bp_%s_%s.pdf" % (g1.name, g2.name)))
