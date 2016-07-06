#!/usr/bin/env python2
import argparse

import sys

import os
import random

import math

from dendropy import Tree
from dendropy.simulate import treesim
import file_ops
import numpy as np
import model

__author__ = 'pfeijao'

# Config names:
TREE_FILE = "%s/input_tree.nwk"
EVOLVED_TREE_FILE = "%s/evolved_tree.nwk"
EVOLVED_BASIC_TREE_FILE = "%s/evolved_basic.nwk"
LEAF_GENOMES_FILE = "%s/leaf_genomes.txt"
ANCESTRAL_GENOMES_FILE = "%s/ancestral_genomes.txt"
LOG_FILE = "%s/simulation.log"
MGRA2_CONFIG = "%s/mgra2.cfg"

# noinspection PyClassHasNoInit
class EventType:
    REARRANGEMENT, DELETION, INSERTION = range(3)

# noinspection PyClassHasNoInit
class RearrangementType:
    REVERSAL, TRANSLOCATION, TRANSPOSITION = range(3)


class Simulation:
    def __init__(self, folder, set_distances=False):
        self.leaf_genomes = file_ops.open_genome_file(LEAF_GENOMES_FILE % folder)
        self.ancestral_genomes = file_ops.open_genome_file(ANCESTRAL_GENOMES_FILE % folder)
        self.evolved_tree = file_ops.open_newick_tree(EVOLVED_TREE_FILE % folder)
        self.input_tree = file_ops.open_newick_tree(TREE_FILE % folder)
        self.folder = folder

## Main: Generate simulation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="Simulates rearrangement evolution on a given newick tree")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the root genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the root genome.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-i",    "--indel", type=float, default=0.0, help="Percentage of indels")
    parser.add_argument("--indel_length", type=int, default=5, help="Maximum size of indel event in genes.")
    scaling = parser.add_mutually_exclusive_group(required=False)
    scaling.add_argument("-r", "--rate", type=float, default=1, help="Multiplier on the input tree number of events.")
    scaling.add_argument("-sc", "--scale", type=float,
                         help="Scales the tree so each leaf has on average scale/2 *n_genes number of events. Almost the same as diameter, if tree is ultrametric. ")
    # TODO: parameter to overwrite stuff; if not present, do not rewrite on same folder
    scaling.add_argument("-e", "--events_per_edge", type=int, help="Events on edge: random [e/2, e]")
    tree_input = parser.add_mutually_exclusive_group(required=True)
    tree_input.add_argument("-f", "--file", help="Input a Newick tree")
    tree_input.add_argument("-s", "--sim", type=int, help="Simulate a new birth_death with SIM species")

    parser.add_argument("-d", "--disturb", type=float, default=0,
                        help="Disturb branch lengths multiplying each by e^r, where r in [-d,+d]. ")
    param = parser.parse_args()

    # Rearrangement, Insertion and Deletion prob:
    rearrangement_p = 1 - param.indel
    insertion_p = param.indel/2
    deletion_p = param.indel/2

    # Rearrangement probabilities: # TODO: include as parameter; as of now, only reversals and transloc.
    if param.num_chr > 1:
        reversal_p = 0.9
        transposition_p = 0
        translocation_p = 0.1
    else:
        reversal_p = 1
        transposition_p = 0
        translocation_p = 0

    # insertion and deletions parameters:
    insertion_length_range = xrange(1, param.indel_length+1)
    deletion_length_range = xrange(1, param.indel_length+1)

    # current insertion genes: (new genes)
    insertion_gene = param.num_genes + 1


    def apply_random_reversal(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = sorted(np.random.choice(len(chromosome.gene_order)+1, 2))
        chromosome.gene_order[bp[0]:bp[1]] = reversed([-x for x in chromosome.gene_order[bp[0]:bp[1]]])


    def apply_random_transposition(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = sorted(np.random.choice(len(chromosome.gene_order)+1, 3))
        chromosome.gene_order[bp[0]:bp[2]] = chromosome.gene_order[bp[1]:bp[2]] + chromosome.gene_order[bp[0]:bp[1]]


    def apply_random_translocation(genome):
        chromosomes = np.random.choice(genome.chromosomes, 2, replace=False)
        bp1 = np.random.choice(len(chromosomes[0].gene_order))
        bp2 = np.random.choice(len(chromosomes[1].gene_order))
        chromosomes[0].gene_order[bp1:], chromosomes[1].gene_order[bp2:] = \
            chromosomes[1].gene_order[bp2:], chromosomes[0].gene_order[bp1:]


    def apply_random_deletion(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(deletion_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        chromosome.gene_order[bp:bp + length] = []
        ## remove chromosome if empty:
        if len(chromosome.gene_order) == 0:
            genome.chromosomes.remove(chromosome)


    def apply_random_insertion(genome, gene):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(insertion_length_range)
        chromosome.gene_order[bp:bp] = range(gene, gene + length)
        return gene + length


    def apply_random_events(genome, n, current_insertion_gene):
        rearrangement_count = 0
        insertion_count = 0
        deletion_count = 0
        # choose events and apply:
        events = np.random.choice([EventType.REARRANGEMENT, EventType.INSERTION, EventType.DELETION], n,
                                  p=[rearrangement_p, insertion_p, deletion_p])

        for event in events:  # number of events, can be weighted by 'scaling' parameters
            if event == EventType.REARRANGEMENT:
                rearrangement = np.random.choice([RearrangementType.REVERSAL, RearrangementType.TRANSLOCATION,
                                                  RearrangementType.TRANSPOSITION], 1,
                                                 p=[reversal_p, translocation_p, transposition_p])
                if rearrangement == RearrangementType.REVERSAL:
                    apply_random_reversal(genome)
                elif rearrangement == RearrangementType.TRANSLOCATION:
                    apply_random_translocation(genome)
                elif rearrangement == RearrangementType.TRANSPOSITION:
                    apply_random_transposition(genome)
                else:
                    raise RuntimeError("Unknown rearrangement type.")
                rearrangement_count += 1

            elif event == EventType.DELETION:
                apply_random_deletion(genome)
                deletion_count += 1
            elif event == EventType.INSERTION:
                current_insertion_gene = apply_random_insertion(genome, current_insertion_gene)
                insertion_count += 1
            else:
                raise RuntimeError("Unknown evolutionary event.")
        return rearrangement_count, insertion_count, deletion_count, current_insertion_gene


    def evolve_tree(ev_tree, root_genome, current_insertion_gene):
        idx = 1
        for ev_node in ev_tree.preorder_node_iter():
            if ev_node.parent_node is None:
                ev_node.value = root_genome.clone("Root")
                ev_node.rearrangements = 0
                ev_node.insertions = 0
                ev_node.deletions = 0
                ev_node.label = "Root"
            else:
                # evolve genome:
                if ev_node.is_internal():
                    ev_node.label = "M%02d" % idx
                    idx += 1
                else:  # complete labelling for leaves
                    ev_node.label = ev_node.taxon.label
                current_genome = ev_node.parent_node.value.clone(ev_node.label)
                weight = ev_node.edge.length

                # evolution
                rearrangements, insertions, deletions, current_insertion_gene = \
                    apply_random_events(current_genome, int(weight), current_insertion_gene)

                ev_node.rearrangements = ev_node.parent_node.rearrangements + rearrangements
                ev_node.deletions = ev_node.parent_node.deletions + deletions
                ev_node.insertions = ev_node.parent_node.insertions + insertions
                ev_node.value = current_genome

        return ev_tree

    num_genes = param.num_genes
    num_chr = param.num_chr
    if param.file is not None:
        tree = Tree.get_from_path(param.file, schema="newick")  # , as_rooted=True)
        tree.reroot_at_midpoint()
    else:
        tree = treesim.birth_death_tree(birth_rate=0.001, death_rate=0, ntax=param.sim)
        tree.seed_node.edge.length = 0
        if param.disturb > 0:
            d = param.disturb
            for edge in tree.postorder_edge_iter():
                r = random.random() * 2 * d - d
                edge.length *= math.exp(r)

    # Scaling

    if param.rate != 1:
        tree.scale_edges(param.rate)
        for edge in tree.postorder_edge_iter():
            edge.length = round(edge.length, 0) if edge.length is not None else 0
    elif param.scale is not None:
        # UPDATE: different scale, like MGRA;
        tree.scale_edges(
                0.5 * param.scale * num_genes / np.mean([node.distance_from_root() for node in tree.leaf_nodes()]))
        for edge in tree.postorder_edge_iter():
            edge.length = round(edge.length, 0) if edge.length is not None else 0
    elif param.events_per_edge is not None:
        for edge in tree.postorder_edge_iter():
            edge.length = random.randint(param.events_per_edge/2, param.events_per_edge)

    output = param.output
    if not os.path.exists(output):
        os.makedirs(output)

    # Output starting tree: (before the sim, because it might change it?)
    tree.write_to_path(TREE_FILE % output, schema='newick', suppress_rooting=True)

    # Simulate!
    evolved = evolve_tree(Tree(tree), model.Genome.identity(num_genes, num_chr), insertion_gene)

    # Output evolved tree:
    evolved.write_to_path(EVOLVED_TREE_FILE % output, schema='newick')
    evolved.write_to_path(EVOLVED_BASIC_TREE_FILE % output, schema='newick',
                          suppress_rooting=True, suppress_edge_lengths=True)

    # create genomes:
    leaf_genomes = {node.taxon.label: node.value for node in evolved.leaf_nodes()}
    ancestral_genomes = {node.label: node.value for node in evolved.internal_nodes()}

    # Output sim result:

    # Genomes:
    file_ops.write_genomes_to_file(leaf_genomes, LEAF_GENOMES_FILE % output)
    file_ops.write_genomes_to_file(ancestral_genomes, ANCESTRAL_GENOMES_FILE % output)

    # Software-specific files:
    # TODO: Other outputs; some are ready from the previous version; I have to implement for the new methods with indel

    # MGRA2:
    file_ops.write_mgra2_config(leaf_genomes, evolved, MGRA2_CONFIG % output)

    #
    # # == Procar format:
    #
    # # Genomes:
    # procar.write_procar_genomes(num_genes, leaf_genomes, PROCAR_GENOMES % output)
    # # Trees
    # procar.write_all_procar_trees(evolved, ancestral_genomes.keys(), PROCAR_TREE % output)
    #
    # # == PATHGROUPs format:
    # ids = pathgroups.genomes_to_pathgroups(leaf_genomes, PATHGROUPS_GENOMES % output)
    # pathgroups.tree_to_pathgroups(evolved, PATHGROUPS_TREE % output, ids)
    #
    # # == GASTS format:
    # gasts.genomes_to_gasts(leaf_genomes, GASTS_GENOMES % output)

    # Log:
    with open(LOG_FILE % output, "w") as f:
        f.write("Num_genes(total)\t%d\n" % num_genes)
        f.write("Num_chromosomes\t%d\n" % num_chr)
        f.write("Num_species\t%d\n" % len(evolved.leaf_nodes()))
        f.write("Input Tree\t%s\n" % (param.file if param.file is not None else "BD tree"))
        f.write("Evol.Rate\t%d\n" % param.rate)
        f.write("Num_events\t%d\n" % int(evolved.length()))
        f.write("Avg events per branch\t%.2f\n" % (evolved.length() / len([x for x in evolved.postorder_edge_iter()])))
        d = {node.taxon.label: node.distance_from_root() for node in tree.leaf_nodes()}
        f.write("Distances to root. Average: %.1f\n" % np.mean(d.values()))

        f.write(",".join(["(%s:%d)" % (node, distance) for node, distance in d.iteritems()]))
        f.write("\n\n")
        f.write(evolved.as_ascii_plot(show_internal_node_labels=True, plot_metric='length') + "\n")
