#!/usr/bin/env python2
import os

import pyximport;

pyximport.install()
import file_ops
import argparse
import re

matching_regexp = re.compile("<variable name=\"x_A(\d+)_(\d+)h,B(\d+)_(\d+)h.*value=\"(.+)\"")
#   <variable name="x_A5_1t,B5_2t" index="150" value="1"/>

def parse_assignment_quality(sol_file, correct_assignment):
    matching = {}
    # open genomes to get the correct matching:
    correct = wrong = 0
    with open(sol_file) as f:
        for l in f:
            m = matching_regexp.match(l.strip())
            if m is not None:
                gene_a, copy_a, gene_b, copy_b, val = m.groups()
                if float(val) >= 0.9 and int(gene_a) > 0:
                    matching[(int(gene_a), int(copy_a))] = int(copy_b)
    for gene_i, copy_i in correct_assignment.iteritems():
        if matching[gene_i] == copy_i:
            correct += 1
        else:
            wrong += 1
    return correct, wrong


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a orthology detetion quality summary file from a list of CPLEX solution files to the DCJdupindel ILP.")
    parser.add_argument("genome_file", type=str, help="Genome file")
    parser.add_argument("copy_number_file", type=str, help="Copy number file")
    # parser.add_argument("files", type=str, nargs="+", help="Solution file(s).")
    param = parser.parse_args()
    # open genome file:
    genomes = file_ops.open_genomes_with_copy_number(param.genome_file, param.copy_number_file)

    for i in range(1, len(genomes)):
        # find the copy_to_ILP-idx assignment in the 2nd genome; the assignment for ILP is sequential
        copy_dict = {}  # dict (gene,copy) -> idx
        genome = genomes[i]
        # init the sequential assignment with 1
        copy_number = {gene: 1 for gene in genome.gene_count()}
        for chrom in genome.chromosomes:
            for gene_i, copy_i in zip(map(abs,chrom.gene_order), chrom.copy_number):
                copy_dict[(gene_i, copy_i)] = copy_number[gene_i]
                copy_number[gene_i] += 1
        # now build the correct genome A IDX -> genome B IDX assignment:
        correct_assignment = {}
        genome = genomes[0]
        copy_number = {gene: 1 for gene in genome.gene_count()}
        for chrom in genome.chromosomes:
            for gene_i, copy_i in zip(map(abs,chrom.gene_order), chrom.copy_number):
                correct_assignment[(gene_i, copy_number[gene_i])] = copy_dict[(gene_i, copy_i)]
                copy_number[gene_i] += 1
        # print correct
        sol_file = "%s_%s_%s.lp.sol" % (param.genome_file, genomes[0].name, genomes[i].name)
        correct, wrong = parse_assignment_quality(sol_file, correct_assignment)
        print "I:%d TP:%d TN:%d" % (i, correct, wrong)
