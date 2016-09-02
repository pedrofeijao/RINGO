#!/usr/bin/env python2
import pyximport;pyximport.install()
import argparse
import re
import networkx as nx
from networkx.algorithms import connected_components

from model import Ext

time_regexp = re.compile('Solution time =\s+([\d\.]+)\s*sec.\s+Iterations')
gap_regexp = re.compile('Current MIP best bound.+gap.+\s([\d.]+)%')
obj_regexp = re.compile("objectiveValue=\"(.+)\"")
balancing_regexp = re.compile("<variable name=\"w_([AB])(\d+)_(\d+)([ht]),([AB])(\d+)_(\d+)([ht]).*value=\"(.+)\"")
# simple_balancing_regexp = re.compile("<variable name=\"w_(.+?),(.+?)\".*value=\"(.+)\"")
#   <variable name="w_A121_2h,A121_2t" index="7191" value="1"/>


def parse_ilp_sol(filename):
    G = nx.Graph()
    genes = {}
    obj = 0
    with open(filename) as f:
        for l in f:
            m = obj_regexp.match(l.strip())
            if m is not None:
                obj = int(round(float(m.group(1)), 0))
            m = balancing_regexp.match(l.strip())
            if m is not None:
                # print m.groups()
                # print l
                genome_1, gene_1, copy_1, ext_1, genome_2, gene_2, copy_2, ext_2, val = m.groups()
                if float(val) >= 0.9:  # sometimes variables with values like 1e-12, 0.99999 or 1.000000001 appear.
                    # save the vertices, to connect head and tail:
                    genes[(genome_1, gene_1, copy_1)] = 1
                    genes[(genome_2, gene_2, copy_2)] = 1
                    G.add_edge((genome_1, gene_1, copy_1, ext_1), (genome_2, gene_2, copy_2, ext_2))

    # connect genes to find connected components corresponding to InDels:
    for genome, gene, copy in genes.iterkeys():
        G.add_edge((genome, gene, copy, Ext.HEAD), (genome, gene, copy, Ext.TAIL))

    indels = len(list(connected_components(G)))
    # log file:

    logfile = filename.replace("sol", "log")
    with open(logfile) as f:
        time = 0
        gap = 0
        for line in f:
            m = time_regexp.match(line.strip())
            if m is not None:
                time = float(m.group(1))

            m = gap_regexp.match(line.strip())
            if m is not None:
                gap = float(m.group(1))

        # print "GAP:", gap
        # print "TIME:", time
    return obj, obj-indels, indels, time, gap

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of CPLEX solution files to the DCJdupindel ILP.")
    parser.add_argument("files", type=str, nargs="+", help="Solution file(s).")
    param = parser.parse_args()
    print "\t".join(["File", "total", "rearrangement", "indel", "time", "gap"])
    for sol_file in param.files:
        r = parse_ilp_sol(sol_file)
        title = sol_file.replace(".lp.sol", "").replace("genomes.txt_", "")
        print "\t".join([title] + [str(x) for x in r])



