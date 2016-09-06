#!/usr/bin/env python2
import collections
import os

import pyximport;

pyximport.install()
from model import Ext
import file_ops
import argparse
import re
import networkx as nx
from networkx.algorithms import connected_components
from simulation import Simulation, EventType
import ringo_config

from parse_orthology import build_gene_copy_assignment, parse_assignment_quality
from coser_to_summary import parse_coser_sol

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

    # indels = len(list(connected_components(G)))
    indels = {"A": 0, "B": 0}
    for comp in connected_components(G):
        genome, gene, copy, ext = list(comp)[0]
        indels[genome] += 1

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

    rearr = obj - indels["A"] - indels["B"]
    return {"dcj_distance": obj, "rearrangements": rearr, "indels_a": indels["A"], "indels_b": indels["B"],
            "time": time, "gap": gap}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of folders for DCJdupindel simulations.")
    parser.add_argument("folders", type=str, nargs="+", help="Sim folders.")
    parser.add_argument("-out", type=str, default="dcj_summary", help="prefix filename for output.")
    parser.add_argument("-coser", action="store_true", default=False, help="Outputs also a COSER summary.")
    param = parser.parse_args()
    # Ringo_config:
    cfg = ringo_config.RingoConfig()
    sim_cols = ["dup_length", "dup_prob", "real_distance"]
    tree_events = ["%s_%s" % (g, event) for g in ["T1", "T2"] for event in EventType.all]
    fields = sim_cols + tree_events + \
             ["dcj_distance", "rearrangements", "indels_a", "indels_b", "time", "gap", "ortho_TP", "ortho_FP"]

    coser_fields = sim_cols + tree_events + \
                   ["dcj_distance", "duplications_a", "duplications_b", "time", "gap", "ortho_TP", "ortho_FP"]

    results = collections.defaultdict(list)
    coser_results = collections.defaultdict(list)
    for folder in param.folders:
        # sim data:
        sim = Simulation.open_folder(folder)
        result = {"dup_length": "%d" % sim.sim_parameters.duplication_length,
                  "dup_prob": "%.2f" % sim.sim_parameters.duplication_p, "folder": folder,
                  "del_prob": "%.2f" % sim.sim_parameters.deletion_p,
                  "real_distance": int(sim.sim_parameters.scale * sim.sim_parameters.num_genes)
                  }
        key = "L%(dup_length)s_dup%(dup_prob)s_del%(del_prob)s" % result
        # tree events:
        result.update(
            {"%s_%s" % (g, event): sim.sim_tree.find_node_with_label(g).edge.events[event] for g in ["T1", "T2"] for
             event in EventType.all})

        # coser is the same until here:
        coser_result = dict(result)

        # distance/rearrangement results:
        sol_file = os.path.join(folder, "extant_genomes.txt_T2_T1.lp.sol")
        # sol_file = os.path.join(folder, "extant_genomes.txt_T2_T1_nobal.lp")
        if not os.path.exists(sol_file):
            continue

        r = parse_ilp_sol(sol_file)
        result.update(r)
        #result["file"] = os.path.basename(sol_file)  # .replace(".lp.sol", "")#.replace("extant_genomes.txt_", "")
        # Orthology:
        genomes = file_ops.open_genomes_with_copy_number(os.path.join(folder, cfg.sim_extant_genomes()))
        correct_assignment = build_gene_copy_assignment(genomes[0], genomes[1])
        correct, wrong = parse_assignment_quality(sol_file, correct_assignment)
        result.update({"ortho_TP": len(correct), "ortho_FP": len(wrong)})

        results[key].append(result)
        # COSER:
        if param.coser:
            coser_result.update(parse_coser_sol(os.path.join(folder, "coser.out")))
            coser_results[key].append(coser_result)


    # output:
    # DCJDUP:
    for key, result in results.iteritems():
        with open("%s_%s.txt" % (param.out, key), "w") as f:
            print >> f, "\t".join(fields)
            for line in sorted(result, key=lambda r: (r['dup_length'], r["dup_prob"], r["real_distance"])):
                print >> f, "\t".join([str(line[field]) for field in fields])
    # COSER:
    if param.coser:
        for key, result in coser_results.iteritems():
            with open("coser_%s_%s.txt" % (param.out,key), "w") as f:
                print >> f, "\t".join(fields)
                for line in sorted(result, key=lambda r: (r['dup_length'], r["dup_prob"], r["real_distance"])):
                    print >> f, "\t".join([str(line[field]) for field in coser_fields])

