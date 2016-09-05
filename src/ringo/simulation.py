#!/usr/bin/env python2
import ringo_config

cfg = ringo_config.RingoConfig()
import pyximport;pyximport.install(build_dir=cfg.pyximport_build())
import argparse
import os
import random
import math
from dendropy import Tree
from dendropy.simulate import treesim
import numpy as np
import algorithms
import file_ops
import model

__author__ = 'pfeijao'


# Config names:

# noinspection PyClassHasNoInit
class EventType:
    all = ["rearrangement", "deletion", "insertion", "duplication"]
    REARRANGEMENT, DELETION, INSERTION, DUPLICATION = all



# noinspection PyClassHasNoInit
class RearrangementType:
    REVERSAL, TRANSLOCATION, TRANSPOSITION = range(3)


class SimParameters:
    def __init__(self, num_genes=0, num_chr=0, del_p=0, ins_p=0, indel_length=0, duplication_p=0,
                 duplication_length=0, rate=0, scale=0, disturb=0, pre_dup=0):
        self.num_genes = num_genes
        self.num_chr = num_chr
        self.deletion_p = del_p
        self.insertion_p = ins_p
        self.indel_length = indel_length
        self.duplication_p = duplication_p
        self.duplication_length = duplication_length
        self.rate = rate
        self.scale = scale
        self.disturb = disturb
        self.pre_duplications = pre_dup

        # Rearrangement, Insertion and Deletion prob:
        self.rearrangement_p = 1 - del_p - ins_p - duplication_p

        assert self.rearrangement_p + self.insertion_p + self.deletion_p + self.duplication_p == 1
        # Rearrangement probabilities: # TODO: include as parameter; as of now, only reversals and transloc.
        if num_chr > 1:
            self.reversal_p = 0.9
            self.transposition_p = 0
            self.translocation_p = 0.1
        else:
            self.reversal_p = 1
            self.transposition_p = 0
            self.translocation_p = 0


class Simulation:
    def __init__(self, folder, sim_parameters=None):
        self.sim_parameters = sim_parameters
        self.sim_tree = None
        self.folder = folder
        self.extant_genomes = None
        self.ancestral_genomes = None

    @staticmethod
    def open_folder(folder):
        sim = Simulation(folder)

        sim.extant_genomes = file_ops.open_genome_file(os.path.join(folder, cfg.sim_extant_genomes()))
        sim.ancestral_genomes = file_ops.open_genome_file(os.path.join(folder, cfg.sim_ancestral_genomes()))
        sim.sim_tree = file_ops.open_newick_tree(os.path.join(folder, cfg.sim_tree()))
        sim.folder = folder

        sim.sim_parameters = SimParameters()
        # Reading data back
        sim.sim_parameters.__dict__ = file_ops.read_simulation_parameters(folder)
        return sim

    @staticmethod
    def apply_random_reversal(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = sorted(np.random.choice(len(chromosome.gene_order) + 1, 2))
        chromosome.gene_order[bp[0]:bp[1]] = reversed([-x for x in chromosome.gene_order[bp[0]:bp[1]]])
        if chromosome.copy_number is not None:
            chromosome.copy_number[bp[0]:bp[1]] = reversed(chromosome.copy_number[bp[0]:bp[1]])

    @staticmethod
    def apply_random_transposition(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = sorted(np.random.choice(len(chromosome.gene_order) + 1, 3))
        chromosome.gene_order[bp[0]:bp[2]] = chromosome.gene_order[bp[1]:bp[2]] + chromosome.gene_order[bp[0]:bp[1]]
        if chromosome.copy_number is not None:
            chromosome.copy_number[bp[0]:bp[2]] = chromosome.copy_number[bp[1]:bp[2]] + chromosome.copy_number[
                                                                                        bp[0]:bp[1]]

    @staticmethod
    def apply_random_translocation(genome):
        chromosomes = np.random.choice(genome.chromosomes, 2, replace=False)
        bp1 = np.random.choice(len(chromosomes[0].gene_order))
        bp2 = np.random.choice(len(chromosomes[1].gene_order))
        chromosomes[0].gene_order[bp1:], chromosomes[1].gene_order[bp2:] = \
            chromosomes[1].gene_order[bp2:], chromosomes[0].gene_order[bp1:]
        if chromosomes[0].copy_number is not None and chromosomes[1].copy_number is not None:
            chromosomes[0].copy_number[bp1:], chromosomes[1].copy_number[bp2:] = \
                chromosomes[1].copy_number[bp2:], chromosomes[0].copy_number[bp1:]

    @staticmethod
    def apply_random_tandem_duplication(genome, duplication_length_range):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(duplication_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        chromosome.gene_order[bp:bp] = chromosome.gene_order[bp:bp + length]

    @staticmethod
    def apply_random_segmental_duplication(genome, duplication_length_range, current_copy_number):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        # length = np.random.choice(duplication_length_range)
        # update to match COSER sims; lenght is fixed at max;
        length = max(duplication_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        # position:
        position = np.random.choice(range(bp) + range(bp + length, chromosome.length()))
        block = chromosome.gene_order[bp:bp + length]
        # update gene copy number
        for gene in block:
            current_copy_number[abs(gene)] += 1
        # apply dup:
        chromosome.copy_number[position:position] = [current_copy_number[abs(x)] for x in block]
        chromosome.gene_order[position:position] = block

    @staticmethod
    def apply_random_deletion(genome, deletion_length_range):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(deletion_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        chromosome.gene_order[bp:bp + length] = []
        ## remove chromosome if empty:
        if len(chromosome.gene_order) == 0:
            genome.chromosomes.remove(chromosome)

    @staticmethod
    def apply_random_insertion(genome, gene, insertion_length_range, current_copy_number):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(insertion_length_range)
        block = range(gene, gene + length)
        chromosome.gene_order[bp:bp] = block
        if chromosome.copy_number is not None:
            for gene in block:
                current_copy_number[abs(gene)] = 1
            chromosome.copy_number[bp:bp] = [1] * length
        return gene + length

    @staticmethod
    def apply_random_rearrangement(param, genome):
        rearrangement = np.random.choice([RearrangementType.REVERSAL, RearrangementType.TRANSLOCATION,
                                          RearrangementType.TRANSPOSITION], 1,
                                         p=[param.reversal_p, param.translocation_p, param.transposition_p])
        if rearrangement == RearrangementType.REVERSAL:
            Simulation.apply_random_reversal(genome)
        elif rearrangement == RearrangementType.TRANSLOCATION:
            Simulation.apply_random_translocation(genome)
        elif rearrangement == RearrangementType.TRANSPOSITION:
            Simulation.apply_random_transposition(genome)
        else:
            raise RuntimeError("Unknown rearrangement type.")

    @staticmethod
    def apply_random_events(param, genome, n, current_insertion_gene, current_copy_number=None):
        rearrangement_count = 0
        insertion_count = 0
        deletion_count = 0
        duplication_count = 0
        # param = self.sim_parameters

        insertion_length_range = xrange(1, param.indel_length + 1)
        deletion_length_range = xrange(1, param.indel_length + 1)
        duplication_length_range = xrange(1, param.duplication_length + 1)

        # choose events and apply:
        event_count = {event:0 for event in EventType.all}
        events = np.random.choice(
            [EventType.REARRANGEMENT, EventType.INSERTION, EventType.DELETION, EventType.DUPLICATION], n,
            p=[param.rearrangement_p, param.insertion_p, param.deletion_p, param.duplication_p])

        for event in events:  # number of events, can be weighted by 'scaling' parameters
            if event == EventType.REARRANGEMENT:
                Simulation.apply_random_rearrangement(param, genome)

            elif event == EventType.DELETION:
                Simulation.apply_random_deletion(genome, deletion_length_range)
            elif event == EventType.INSERTION:
                current_insertion_gene = Simulation.apply_random_insertion(genome, current_insertion_gene,
                                                                           insertion_length_range, current_copy_number)
            elif event == EventType.DUPLICATION:
                # Simulation.apply_random_tandem_duplication(genome, duplication_length_range)
                Simulation.apply_random_segmental_duplication(genome, duplication_length_range, current_copy_number)

            else:
                raise RuntimeError("Unknown evolutionary event.")
            event_count[event] += 1
        return event_count, current_insertion_gene

    def run_simulation(self):

        param = self.sim_parameters
        # insertion and deletions parameters:

        # current insertion genes: (new genes)
        current_insertion_gene = param.num_genes + 1
        current_copy_number = None  # will init at root
        idx = 1
        ev_tree = self.sim_tree
        for ev_node in ev_tree.preorder_node_iter():
            if ev_node.parent_node is None:
                # identity genome:
                ev_node.value = current_genome = model.Genome.identity(param.num_genes, param.num_chr)
                ev_node.events = {ev:0 for ev in EventType.all}

                # add copy number information to track orthologous/paralogous, when duplications are present:
                for chromosome in current_genome.chromosomes:
                    chromosome.copy_number = [1] * len(chromosome.gene_order)
                current_copy_number = current_genome.gene_count()
                # pre-duplications:
                for i in range(param.pre_duplications):
                    Simulation.apply_random_segmental_duplication(current_genome,
                                                                  range(1, param.duplication_length + 1),
                                                                  current_copy_number)
                ev_node.events[EventType.DUPLICATION] = param.pre_duplications

                if ev_node.label is None:
                    ev_node.label = "Root"
            else:
                # evolve genome:
                if ev_node.is_internal():
                    if ev_node.label is None:
                        ev_node.label = "M%02d" % idx
                        idx += 1
                else:  # complete labelling for leaves
                    ev_node.label = ev_node.taxon.label

                current_genome = ev_node.parent_node.value.clone(ev_node.label)
                ev_node.value = current_genome
                weight = ev_node.edge.length

                # evolution
                ev_count, current_insertion_gene = \
                    Simulation.apply_random_events(param, current_genome, int(weight), current_insertion_gene,
                                                   current_copy_number)

                # update count of events at node (from the root) and at each edge
                ev_node.events = {ev: ev_node.parent_node.events[ev] + count for ev, count in ev_count.iteritems()}
                ev_node.edge.events = ev_count



    # Custom simulation on the model of COSER paper, with D DCJs,
    def open_tree(self, treefile):
        self.sim_tree = Tree.get_from_path(treefile, schema="newick")  # , as_rooted=True)
        self.sim_tree.reroot_at_midpoint()

    def simulate_tree(self, ntax):
        tree = treesim.birth_death_tree(birth_rate=0.001, death_rate=0, ntax=ntax)
        tree.seed_node.edge.length = 0
        # disturb and scaling:
        try:
            if param.disturb > 0:
                d = param.disturb
                for edge in tree.postorder_edge_iter():
                    r = random.random() * 2 * d - d
                    edge.length *= math.exp(r)

            # Scaling
            if self.sim_parameters.rate != 1:
                tree.scale_edges(param.rate)
                for edge in tree.postorder_edge_iter():
                    edge.length = round(edge.length, 0) if edge.length is not None else 0
            elif self.sim_parameters.scale is not None:
                diameter = algorithms.tree_diameter(tree)
                tree.scale_edges(self.sim_parameters.scale * self.sim_parameters.num_genes / diameter)

                # round to integer
                for edge in tree.postorder_edge_iter():
                    edge.length = round(edge.length, 0) if edge.length is not None else 0
        except NameError:
            pass
        self.sim_tree = tree

    def save_simulation(self, save_copies=False):
        # Output sim result:
        output = self.folder
        if not os.path.exists(output):
            os.makedirs(output)

        # Output simulated tree:
        tree = self.sim_tree
        tree.write_to_path(os.path.join(output, cfg.sim_tree()), suppress_rooting=True, schema='newick')
        tree.write_to_path(os.path.join(output, cfg.sim_tree_no_lengths()), suppress_rooting=True,
                           suppress_edge_lengths=True, schema='newick')

        # create genomes:
        self.extant_genomes = {node.taxon.label: node.value for node in tree.leaf_nodes()}
        self.ancestral_genomes = {node.label: node.value for node in tree.internal_nodes()}

        # Genomes:
        file_ops.write_genomes_to_file(self.extant_genomes, os.path.join(output, cfg.sim_extant_genomes()))
        file_ops.write_genomes_to_file(self.ancestral_genomes, os.path.join(output, cfg.sim_ancestral_genomes()))

        # Copy number, if duplications present:
        if save_copies:
            leaf_cn = cfg.sim_extant_genomes() + cfg.copy_number_file_extension()
            file_ops.write_genomes_copy_number_to_file(self.extant_genomes, os.path.join(output, leaf_cn))
            ancestral_cn = cfg.sim_ancestral_genomes() + cfg.copy_number_file_extension()
            file_ops.write_genomes_to_file(self.ancestral_genomes, os.path.join(output, ancestral_cn))

        # Software-specific files:
        # MLGO, and maybe others, doesn't like the "#" lines; output a 'simple' version:
        file_ops.write_genomes_to_file(self.extant_genomes, os.path.join(output, cfg.sim_leaf_simple_genomes()),
                                       write_chr_line=False)

        #COSER:
        file_ops.write_genomes_coser_format(sim.extant_genomes, output)

        # MGRA2:
        file_ops.write_mgra2_config(self.extant_genomes, tree, os.path.join(output, cfg.sim_mgra_config()))

        # TODO: If no indel, also write output for Procar, Pathgroups and GASTS. Other softs?
        # # == Procar format:
        #
        # # Genomes:
        # procar.write_procar_genomes(num_genes, extant_genomes, PROCAR_GENOMES % output)
        # # Trees
        # procar.write_all_procar_trees(evolved, ancestral_genomes.keys(), PROCAR_TREE % output)
        #
        # # == PATHGROUPs format:
        # ids = pathgroups.genomes_to_pathgroups(extant_genomes, PATHGROUPS_GENOMES % output)
        # pathgroups.tree_to_pathgroups(evolved, PATHGROUPS_TREE % output, ids)
        #
        # # == GASTS format:
        # gasts.genomes_to_gasts(extant_genomes, GASTS_GENOMES % output)

        # Log:
        param = self.sim_parameters
        with open(os.path.join(output, cfg.sim_logfile()), "w") as f:
            f.write("Num_genes(at the root genome)\t%d\n" % param.num_genes)
            f.write("Num_chromosomes\t%d\n" % param.num_chr)
            f.write("Num_species\t%d\n" % len(tree.leaf_nodes()))
            # f.write("Input Tree\t%s\n" % (param.file if param.file is not None else "BD tree"))
            f.write("Evol.Rate\t%d\n" % param.rate)
            f.write("Num_events\t%d\n" % int(tree.length()))
            f.write("Avg events per branch\t%.2f\n" % (tree.length() / (len(list(tree.postorder_edge_iter())) - 1)))
            d = {node.taxon.label: node.distance_from_root() for node in tree.leaf_nodes()}
            f.write("Diameter: %.1f\n" % algorithms.tree_diameter(tree))
            f.write(",".join(["(%s:%d)" % (node, distance) for node, distance in d.iteritems()]))
            f.write("\n\n")
            f.write(tree.as_ascii_plot(show_internal_node_labels=True, plot_metric='length') + "\n")

        # Events per edge:
        with open(os.path.join(sim.folder, cfg.sim_events_file()), "w") as f:
            f.write("Edge\t%s\n" % "\t".join(EventType.all))
            for node in tree.postorder_node_iter():
                if node == tree.seed_node:
                    continue
                f.write("%s\t%s\n" % (node.label, "\t".join(str(node.edge.events[ev]) for ev in EventType.all)))

        # Save parameters:
        file_ops.write_simulation_parameters(param, output)

## Main: Generate simulation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulates rearrangement evolution on a given newick tree")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the root genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the root genome.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0, help="Percentage of deletions, from 0 to 1.0")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0, help="Percentage of insertions, from 0 to 1.0")
    parser.add_argument("-dup_p", "--duplication_p", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0")
    parser.add_argument("-dl", "--duplication_length", type=int, default=5, help="Maximum length of duplication event.")
    parser.add_argument("-il", "--indel_length", type=int, default=5, help="Maximum size of indel event in genes.")
    scaling = parser.add_mutually_exclusive_group(required=False)
    scaling.add_argument("-r", "--rate", type=float, default=1, help="Multiplier on the input tree number of events.")
    scaling.add_argument("-sc", "--scale", type=float,
                         help="Scales the tree so each leaf has on average scale/2 *n_genes number of events. Almost the same as diameter, if tree is ultrametric. ")
    # TODO: parameter to overwrite stuff; if not present, do not rewrite on same folder
    parser.add_argument("--predup", type=int, default=0,
                        help="Apply a number of duplications in the initial genome, before the random walk.")

    tree_input = parser.add_mutually_exclusive_group(required=True)
    tree_input.add_argument("-f", "--file", help="Input a Newick tree")
    tree_input.add_argument("-s", "--sim", type=int, help="Simulate a new birth_death with SIM species")

    parser.add_argument("-d", "--disturb", type=float, default=0,
                        help="Disturb branch lengths multiplying each by e^r, where r in [-d,+d]. ")
    param = parser.parse_args()

    # Simulation parameters:
    sim_par = SimParameters(num_genes=param.num_genes, num_chr=param.num_chr,
                            del_p=param.deletion_p, ins_p=param.insertion_p, indel_length=param.indel_length,
                            duplication_p=param.duplication_p, duplication_length=param.duplication_length,
                            rate=param.rate, scale=param.scale, disturb=param.disturb, pre_dup=param.predup)
    # start sim object;
    sim = Simulation(param.output, sim_par)

    # sim tree:
    if param.file is not None:
        sim.open_tree(param.file)
    else:
        sim.simulate_tree(param.sim)

    sim.run_simulation()

    sim.save_simulation(save_copies=sim_par.pre_duplications > 0 or sim_par.duplication_p > 0)
