#!/usr/bin/env python
# -*- coding: utf-8
#
import argparse

from dendropy import Taxon, Tree
import file_ops
from model import Genome

__author__ = 'pedro'


def scj_bottom_up_fitch(tree, genome_adj_set, a):
    for node in tree.postorder_node_iter():
        if node.is_internal():
            g1, g2 = node.child_nodes()
            if g1.set.isdisjoint(g2.set):
                node.set = g1.set.union(g2.set)
            else:
                node.set = g1.set.intersection(g2.set)

        else:  # leaf
            node.set = {a in genome_adj_set[node.label]}
    return tree


def scj_top_down_fitch(tree, adjacencies, adj, root_value=False):
    # if a == [0,2]:
    #     import pdb; pdb.set_trace()
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue

        # Root: choose False if we have both:
        if node.parent_node is None:
            if len(node.set) > 1:  # Has False and True
                node.set = root_value  # might be True, leading to conflict.
            else:
                node.set = node.set.pop()

        # Internal: choose own if unitary set, otherwise choose parent.
        else:
            if len(node.set) == 1:
                node.set = node.set.pop()
            else:
                node.set = node.parent_node.set

        # If it is True, add adjacency to genome:
        if node.set:
            adjacencies[node.label].append(adj)
    return tree


def scj_small_phylogeny_adjacencies(tree, genome_list, root_value=False):
    return _scj_small_phylogeny(tree, genome_list, genes=False, root_value=root_value)


def scj_small_phylogeny_genes(tree, genome_list, root_value=True):
    return _scj_small_phylogeny(tree, genome_list, genes=True, root_value=root_value)


def _scj_small_phylogeny(tree, genome_list, genes=False, root_value=False):
    # find all adjacencies in leaves:
    all_adjacencies = set()
    genome_adj_set = dict()
    for name, genome in genome_list.iteritems():
        if genes:
            genome_adj_set[name] = genome.gene_set()
        else:
            genome_adj_set[name] = genome.adjacency_set()
        all_adjacencies = all_adjacencies.union(genome_adj_set[genome.name])

    adj_set = {node.label: [] for node in tree.internal_nodes()}
    for adj in all_adjacencies:
        scj_bottom_up_fitch(tree, genome_adj_set, adj)
        scj_top_down_fitch(tree, adj_set, adj, root_value)

    # Filter conflicts: TODO: maybe a smarter way to filter during the process?
    if not genes and root_value:
        new_adj_set = {}
        for key, adjacencies in adj_set.iteritems():
            new_adj = []
            seen = set()
            for a in adjacencies:
                if a[0] in seen or a[1] in seen:
                    continue
                seen.add(a[0])
                seen.add(a[1])
                new_adj.append(a)
            new_adj_set[key] = new_adj
        adj_set = new_adj_set

    return adj_set

#TODO: write a "solve SCJ script."
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="SCJ algorithms.")
    parser.add_argument("sim_name", type=str, help="Name of simulation to solve.")
    param = parser.parse_args()

    # # open sim:
    # sim = simulations.Simulation.open_folder(param.sim_name)
    # tree, adjacencies = scj_small_phylogeny(sim.n_genes, sim.sim_tree, sim.leaf_genomes)
    #
    # file_ops.write_genomes_to_file(sim.n_genes, {l:Genome(2*sim.n_genes,adj) for l,adj in adjacencies.iteritems()}, simulations.SCJ_OUTPUT % sim.folder )
