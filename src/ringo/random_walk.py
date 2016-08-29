#!/usr/bin/env python2
import os
import sys

import ringo_config;

cfg = ringo_config.RingoConfig()
import argparse
from simulation import SimParameters, Simulation
from model import Genome
import file_ops


class RandomWalkParams(SimParameters):
    def __init__(self, pre_dup=0, *args, **kwargs):
        SimParameters.__init__(self, *args, **kwargs)
        self.pre_duplications = pre_dup


class RandomWalk:
    number_of_events = ["total", "rearrangement", "insertion", "deletion", "duplication"]

    def __init__(self, folder, length, steps, sim_parameters=None):
        self.sim_parameters = sim_parameters
        self.folder = folder
        self.length = length
        self.steps = steps
        self.genomes = []
        self.events = {event: [] for event in RandomWalk.number_of_events}

    def do_the_walk(self):
        sim_param = self.sim_parameters

        # current insertion genes: (new genes)
        current_insertion_gene = sim_param.num_genes + 1

        # start genome:
        current_genome = Genome.identity(sim_param.num_genes, sim_param.num_chr, name="G_0")

        # do some pre-dups if necessary:
        if sim_param.pre_duplications > 0:
            for i in range(sim_param.pre_duplications):
                Simulation.apply_random_segmental_duplication(current_genome, range(1, param.duplication_length+1))

        self.genomes.append(current_genome)

        for key in RandomWalk.number_of_events:
            self.events[key].append(0)

        # walk:
        for step in range(self.length):
            # apply random event on current;
            current_genome = current_genome.clone("G_%d" % (step+1))
            n_rearrangements, n_insertions, n_deletions, n_duplications, current_insertion_gene = \
                Simulation.apply_random_events(sim_param, current_genome, self.steps, current_insertion_gene)
            for key, value in zip(RandomWalk.number_of_events,
                                  [self.steps, n_rearrangements, n_insertions, n_deletions, n_duplications]):
                self.events[key].append(value)
            self.genomes.append(current_genome)

    def save_genomes(self, repeat):
        # create folder:
        output = os.path.join(self.folder, "rep_%0d" % (repeat + 1))
        if not os.path.exists(output):
            os.makedirs(output)

        # save genomes:
        file_ops.write_genomes_to_file(self.genomes, os.path.join(output, "genomes.txt"))

        # save in COSER format:
        file_ops.write_genomes_coser_format(self.genomes, output)

        # save distances:
        with open(os.path.join(output, "distances.txt"), "w") as f:
            f.write("%s\n" % ("\t".join(RandomWalk.number_of_events)))
            for idx, g in enumerate(self.genomes):
                f.write("%s\t%s\n" % (
                g.name, "\t".join([str(sum(self.events[k][:(idx + 1)])) for k in RandomWalk.number_of_events])))

        # Save parameters:
        file_ops.write_simulation_parameters(param, output)


## Main: Generate random walk
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulates a random walk with rearrangement, indel and duplication evolution.")
    parser.add_argument("-l", "--length", type=int, required=True, help="Length of the random walk.")
    parser.add_argument("-s", "--steps", type=int, required=True,
                        help="Number of mutation in each step of the random walk.")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the start genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the start genome.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0, help="Percentage of indels, from 0 to 1.0")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0, help="Percentage of indels, from 0 to 1.0")
    parser.add_argument("-dup_p", "--duplication_p", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0")
    parser.add_argument("-dl", "--duplication_length", type=int, default=5, help="Maximum length of duplication event.")
    parser.add_argument("-il", "--indel_length", type=int, default=5, help="Maximum size of indel event in genes.")
    parser.add_argument("-r", "--repeats", type=int, default=1, help="How many repetitions.")
    parser.add_argument("--predup", type=int, default=0,
                        help="Apply a number of duplications in the initial genome, before the random walk.")
    param = parser.parse_args()

    sim_par = RandomWalkParams(pre_dup=param.predup, num_genes=param.num_genes, num_chr=param.num_chr,
                            del_p=param.deletion_p, ins_p=param.insertion_p, indel_length=param.indel_length,
                            duplication_p=param.duplication_p, duplication_length=param.duplication_length)

    for i in range(param.repeats):
        if (i + 1) % 1 == 0:
            print >> sys.stderr, "Generating walk %d ... " % (i + 1)
        random_walk = RandomWalk(param.output, param.length, param.steps, sim_parameters=sim_par)
        random_walk.do_the_walk()
        random_walk.save_genomes(i)
    print >> sys.stderr, "Finished!"
