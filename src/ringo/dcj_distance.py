#!/usr/bin/env python2
import argparse
import pyximport;
pyximport.install()
import file_ops
import dcj

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Returns the DCJ distance between 2 genomes.")
    parser.add_argument("g", type=str, nargs=3, help="GRIMM genomes file, idx 1 and 2 of genomes (0-indexed).")
    param = parser.parse_args()
    filename, n1, n2 = param.g
    genomes = file_ops.open_genome_file(filename, as_list=True)
    g1 = genomes[int(n1)]
    g2 = genomes[int(n2)]
    print dcj.dcj_distance(g1, g2)

