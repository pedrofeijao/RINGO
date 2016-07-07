#!/usr/bin/env python2
import pyximport; pyximport.install()
import argparse
import os
import scj
import file_ops
from model import Genome
from ringo_config import RingoConfig
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SCJ algorithms.")
    parser.add_argument("-i", "--input_genomes", required=True, type=str, help="Leaf genomes file.")
    parser.add_argument("-t", "--tree", required=True, type=str, help="Newick Tree file.")
    parser.add_argument("-o", "--output", type=str, help="Output folder. If not given, output is written to the same location of the genomes file.")
    param = parser.parse_args()

    # open files:
    extant_genomes = file_ops.open_genome_file(param.input_genomes)
    tree = file_ops.open_newick_tree(param.tree, label_internal_nodes=True)

    # run SCK
    adj_set = scj.scj_small_phylogeny_adjacencies(tree, extant_genomes)
    genomes = {label:Genome.from_adjacency_list(label, adj) for label, adj in adj_set.iteritems()}

    # write output:
    cfg = RingoConfig()
    folder = param.output if param.output is not None else os.path.dirname(param.input_genomes)
    file_ops.write_genomes_to_file(genomes, os.path.join(folder, cfg.scj_genomes()))
