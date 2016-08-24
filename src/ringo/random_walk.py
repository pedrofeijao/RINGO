#!/usr/bin/env python2
import os

import sys

import ringo_config
cfg = ringo_config.RingoConfig()
import argparse
from simulation import SimParameters, Simulation
from model import Genome
import file_ops

class RandomWalk:
    def __init__(self, folder, length, steps, sim_parameters=None):
        self.sim_parameters = sim_parameters
        self.folder = folder
        self.length = length
        self.steps = steps
        self.genomes = []
        self.distances = []

    def do_the_walk(self):
        sim_param = self.sim_parameters

        # current insertion genes: (new genes)
        current_insertion_gene = sim_param.num_genes + 1

        # start genome:
        current_genome = Genome.identity(sim_param.num_genes, sim_param.num_chr, name="G_0")
        distance = 0
        self.genomes.append(current_genome)
        self.distances.append(distance)
        # walk:

        for step in range(self.length):
            distance += self.steps
            # apply random event on current;
            current_genome = current_genome.clone("G_%d" % distance)
            rearrangements, insertions, deletions, duplications, current_insertion_gene = \
                Simulation.apply_random_events(sim_param, current_genome, self.steps, current_insertion_gene)
            self.genomes.append(current_genome)
            self.distances.append(distance)
    def save_genomes(self, repeat):
        # create folder:
        output = os.path.join(self.folder, "rep_%0d" % (repeat+1))
        if not os.path.exists(output):
            os.makedirs(output)

        # save genomes:
        file_ops.write_genomes_to_file(self.genomes, os.path.join(output, "genomes.txt"))

        # save in COSER format:
        file_ops.write_genomes_coser_format(self.genomes, output)

        # save distances:
        with open(os.path.join(output, "distances.txt"), "w") as f:
            for g, dist in zip(self.genomes, self.distances):
                f.write("%s\t%d\n" % (g.name, dist))

        # Save parameters:
        file_ops.write_simulation_parameters(param, output)


## Main: Generate random walk
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="Simulates a random walk with rearrangement, indel and duplication evolution.")
    parser.add_argument("-l", "--length", type=int, required=True, help="Length of the random walk.")
    parser.add_argument("-s", "--steps", type=int, required=True, help="Number of mutation in each step of the random walk.")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the start genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the start genome.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0, help="Percentage of indels, from 0 to 1.0")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0, help="Percentage of indels, from 0 to 1.0")
    parser.add_argument("-dup_p", "--duplication_p", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0")
    parser.add_argument("-dl",  "--duplication_length", type=int, default=5, help="Maximum length of duplication event.")
    parser.add_argument("-il", "--indel_length", type=int, default=5, help="Maximum size of indel event in genes.")
    parser.add_argument("-r", "--repeats", type=int, default=1, help="How many repetitions.")
    param = parser.parse_args()

    sim_par = SimParameters(num_genes=param.num_genes, num_chr=param.num_chr,
                            del_p=param.deletion_p, ins_p=param.insertion_p, indel_length=param.indel_length,
                            duplication_p=param.duplication_p, duplication_length=param.duplication_length)

    for i in range(param.repeats):
        if (i+1) % 1 == 0:
            print >> sys.stderr, "Generating walk %d ... " % (i+1)
        random_walk = RandomWalk(param.output, param.length, param.steps, sim_parameters=sim_par)
        random_walk.do_the_walk()
        random_walk.save_genomes(i)
    print >> sys.stderr, "Finished!"

