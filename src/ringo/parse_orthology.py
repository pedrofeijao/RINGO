#!/usr/bin/env python2
import collections
import copy
import os

import pyximport;

pyximport.install()
import file_ops
import argparse
import re

matching_regexp = re.compile("<variable name=\"x_A(\d+)_(\d+)h,B(\d+)_(\d+)h.*value=\"(.+)\"")


#   <variable name="x_A5_1t,B5_2t" index="150" value="1"/>


def build_correct_matching(genome_a, genome_b):
    all_genes = genome_a.gene_set().union(genome_b.gene_set())
    gene_copies_a = genome_a.gene_copies()
    gene_copies_b = genome_b.gene_copies()

    correct_matching = {gene: gene_copies_a[gene].intersection(gene_copies_b[gene]) for gene in all_genes if
                        len(gene_copies_a[gene]) > 1 or len(gene_copies_b[gene]) > 1}
    return correct_matching


def parse_orthology_quality(solution_matching, correct_matching):
    tp = set()
    fp = set()
    not_matched = copy.deepcopy(correct_matching)
    for gene, pair_list in solution_matching.iteritems():
        if gene not in correct_matching:  # if it is not, it a single copy gene, skip;
            continue
        for (c_a, c_b) in pair_list:
            # if both copies are not there, it is a balancing match;
            if c_a not in correct_matching[gene] and c_b not in correct_matching[gene]:
                continue
            if c_a == c_b and c_a in correct_matching[gene]:
                tp.add((gene, c_a, c_b))
                not_matched[gene].remove(c_a)
            else:
                fp.add((gene, c_a, c_b))

    fn = [(g, c_a, c_a) for g, copy_set in not_matched.iteritems() for c_a in copy_set]
    return tp, fp, fn


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a orthology detetion quality summary file from a list of CPLEX solution files to the DCJdupindel ILP for random walk.")
    parser.add_argument("genome_file", type=str, help="Genome file")
    param = parser.parse_args()
    # open genome file:
    genomes = file_ops.open_genomes_with_copy_number(param.genome_file)

    #
    # for i in range(1, len(genomes)):
    #
    #     sol_file = "%s_%s_%s.lp.sol" % (param.genome_file, genomes[0].name, genomes[i].name)
    #     correct_assignment = build_gene_copy_assignment(genomes[0], genomes[i])
    #     correct, wrong = parse_assignment_quality(sol_file, correct_assignment)
    #
    #     print "I:%d TP:%d TN:%d" % (i, len(correct), len(wrong))
    #     print wrong
