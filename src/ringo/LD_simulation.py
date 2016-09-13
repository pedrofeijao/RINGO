#!/usr/bin/env python2
import ringo_config

cfg = ringo_config.RingoConfig()
import pyximport;pyximport.install(build_dir=cfg.pyximport_build())
import argparse
import random
import numpy as np
import model
from simulation import Simulation, SimParameters, EventType, RearrangementType


def run_L_D_simulation(self, L, D):
    # L = duplication length
    # D = number of DCJs in each branch.
    #
    param = self.sim_parameters
    # pre_dups (at root) and post_dups (at branches) to achieve 1.5 genes/family in average.
    pre_duplications = int(0.45 * param.num_genes / L)
    post_duplications = int(0.05 * param.num_genes / L)
    post_duplications = [int(0.6 * post_duplications), int(1.4 * post_duplications)]
    param.pre_duplications = pre_duplications
    current_copy_number = None  # will init at root
    deletion_length_range = xrange(1, param.indel_length + 1)
    duplication_length_range = xrange(1, L + 1)
    idx = 1
    ev_tree = self.sim_tree
    for ev_node in ev_tree.preorder_node_iter():
        if ev_node.parent_node is None:
            # identity genome:
            ev_node.value = current_genome = model.Genome.identity(param.num_genes, param.num_chr)
            ev_node.events = {ev: 0 for ev in EventType.all}

            # add copy number information to track orthologous/paralogous, when duplications are present:
            for chromosome in current_genome.chromosomes:
                chromosome.copy_number = [1] * len(chromosome.gene_order)
            current_copy_number = current_genome.gene_count()
            # pre-duplications:
            for i in range(pre_duplications):
                Simulation.apply_random_segmental_duplication(current_genome,
                                                              range(1, param.duplication_length + 1),
                                                              current_copy_number)

            ev_node.events[EventType.DUPLICATION] = pre_duplications
            # ev_node.edge.length = pre_duplications

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
            pd = post_duplications.pop()
            ev_node.edge.length = D + pd

            # events
            events = [EventType.DUPLICATION] * pd + [EventType.REARRANGEMENT] * D

            ev_node.edge.events = {ev: 0 for ev in EventType.all}
            random.shuffle(events)
            for event in events:
                if event == EventType.DUPLICATION:
                    Simulation.apply_random_segmental_duplication(current_genome, duplication_length_range, current_copy_number)
                    ev_node.edge.events[event] += 1
                elif event == EventType.REARRANGEMENT:
                    # here, I can also have deletions:
                    ev = np.random.choice([RearrangementType.REVERSAL, EventType.DELETION], 1,
                                          p=[param.rearrangement_p, param.deletion_p])[0]
                    if ev == RearrangementType.REVERSAL:
                        Simulation.apply_random_reversal(current_genome)
                        ev_node.edge.events[event] += 1
                    else:
                        Simulation.apply_random_deletion(current_genome, deletion_length_range)
                        ev_node.edge.events[EventType.DELETION] += 1

            ev_node.events = {ev: ev_node.parent_node.events[ev] + count for ev, count in
                              ev_node.edge.events.iteritems()}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulates rearrangement evolution on a given newick tree")
    parser.add_argument("-s", "--sim", type=int, help="Simulate a new birth_death with SIM species")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the root genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the root genome.")
    parser.add_argument("-L", "-dl", "--duplication_length", type=int, default=5, help="Maximum length of duplication event.")
    parser.add_argument("-D", "--rearrangements", type=int, default=5, help="Number of rearrangements.")

    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0, help="Percentage of deletions, from 0 to 1.0")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0, help="Percentage of insertions, from 0 to 1.0")
    parser.add_argument("-il", "--indel_length", type=int, default=5, help="Maximum size of indel event in genes.")

    parser.add_argument("-d", "--disturb", type=float, default=0,
                        help="Disturb branch lengths multiplying each by e^r, where r in [-d,+d]. ")
    param = parser.parse_args()

    # Simulation parameters:
    sim_par = SimParameters(num_genes=param.num_genes, num_chr=param.num_chr,
                            del_p=param.deletion_p, ins_p=param.insertion_p, indel_length=param.indel_length,
                            duplication_length=param.duplication_length)
    # start sim object;
    sim = Simulation(param.output, sim_par)

    sim.simulate_tree(param.sim)

    run_L_D_simulation(sim, param.duplication_length, param.rearrangements)

    sim.save_simulation(save_copies=True)

