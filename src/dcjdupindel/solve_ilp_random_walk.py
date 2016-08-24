#!/usr/bin/env python2
import pyximport; pyximport.install()
import os
import sys
from ilp import dcj_dupindel_ilp

os.sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),"../ringo")))
import file_ops

if __name__ == '__main__':
    distances = []
    genomes = file_ops.open_genome_file(sys.argv[1], as_list=True)
    for i in range(len(genomes)-1):
        print >> sys.stderr, "Solving %d ..." % i
        model = dcj_dupindel_ilp(genomes[0], genomes[i])
        distances.append((genomes[i].name, model.objVal))
    print >> sys.stderr, "\n".join(["%s\t%d" % (name, o) for name,o in distances])

