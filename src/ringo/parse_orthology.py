#!/usr/bin/env python2
import collections
import os

import pyximport;

pyximport.install()
import file_ops
import argparse
import re

matching_regexp = re.compile("<variable name=\"x_A(\d+)_(\d+)h,B(\d+)_(\d+)h.*value=\"(.+)\"")
#   <variable name="x_A5_1t,B5_2t" index="150" value="1"/>


def parse_assignment_quality(sol_file, genome_a, genome_b):
    all_genes = genome_a.gene_set().union(genome_b.gene_set())
    gene_copies_a = genome_a.gene_copies()
    gene_copies_b = genome_b.gene_copies()

    correct_matching = {}
    for gene in all_genes:
        if len(gene_copies_a[gene]) > 1 or len(gene_copies_b[gene]) > 1:
            correct_matching[gene] = gene_copies_a[gene].intersection(gene_copies_b[gene])

    matching = collections.defaultdict(list)
    # open genomes to get the correct matching:
    correct = set()
    wrong = set()
    with open(sol_file) as f:
        for l in f:
            m = matching_regexp.match(l.strip())
            if m is not None:
                gene_a, copy_a, gene_b, copy_b, val = m.groups()
                if float(val) >= 0.9 and int(gene_a) > 0:
                    matching[int(gene_a)].append((int(copy_a),int(copy_b)))

    for gene, pair_list in matching.iteritems():
        if gene not in correct_matching: # if it is not, it a single copy gene, skip;
            continue
        for (c_a, c_b) in pair_list:
            # if both copies are not there, it is a balancing match;
            if c_a not in correct_matching[gene] and c_b not in correct_matching[gene]:
                continue
            if c_a == c_b and c_a in correct_matching[gene]:
                correct.add((gene, c_a))
                correct_matching[gene].remove(c_a)
            else:
                print gene, c_a, c_b
                wrong.add((gene, c_a))

    return correct, wrong, correct_matching


def build_gene_copy_assignment(genome_a, genome_b):
    # find the copy_to_ILP-idx assignment in the 2nd genome; the assignment for ILP is sequential
    copy_dict = {}  # dict (gene,copy) -> idx
    genome = genome_b
    # init the sequential assignment with 1
    copy_number = {gene: 1 for gene in genome.gene_count()}
    for chrom in genome.chromosomes:
        for gene_i, copy_i in zip(map(abs, chrom.gene_order), chrom.copy_number):
            copy_dict[(gene_i, copy_i)] = copy_number[gene_i]
            copy_number[gene_i] += 1
    # now build the correct genome A IDX -> genome B IDX assignment:
    correct_assignment = {}
    genome = genome_a
    copy_number = {gene: 1 for gene in genome.gene_count()}
    for chrom in genome.chromosomes:
        for gene_i, copy_i in zip(map(abs, chrom.gene_order), chrom.copy_number):
            if (gene_i, copy_i) in copy_dict:
                correct_assignment[(gene_i, copy_number[gene_i])] = copy_dict[(gene_i, copy_i)]
            copy_number[gene_i] += 1
    return correct_assignment


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a orthology detetion quality summary file from a list of CPLEX solution files to the DCJdupindel ILP for random walk.")
    parser.add_argument("genome_file", type=str, help="Genome file")
    param = parser.parse_args()
    # open genome file:
    genomes = file_ops.open_genomes_with_copy_number(param.genome_file)

    for i in range(1, len(genomes)):

        sol_file = "%s_%s_%s.lp.sol" % (param.genome_file, genomes[0].name, genomes[i].name)
        correct_assignment = build_gene_copy_assignment(genomes[0], genomes[i])
        correct, wrong = parse_assignment_quality(sol_file, correct_assignment)

        print "I:%d TP:%d TN:%d" % (i, len(correct), len(wrong))
        print wrong
