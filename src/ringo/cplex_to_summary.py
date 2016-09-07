#!/usr/bin/env python2
import collections
import os

import pyximport;

pyximport.install()
from model import Ext, Chromosome
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
matching_regexp = re.compile("<variable name=\"x_A(\d+)_(\d+)h,B(\d+)_(\d+)h.*value=\"(.+)\"")


# simple_balancing_regexp = re.compile("<variable name=\"w_(.+?),(.+?)\".*value=\"(.+)\"")
#   <variable name="w_A121_2h,A121_2t" index="7191" value="1"/>


def parse_log_file(logfile):
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
    return gap, time


def parse_ilp_sol(filename):

    G = nx.Graph()
    obj = 0
    edges = []
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
                    # connect head and tail:
                    G.add_edge((genome_1, gene_1, copy_1, Ext.HEAD), (genome_1, gene_1, copy_1, Ext.TAIL))
                    G.add_edge((genome_2, gene_2, copy_2, Ext.HEAD), (genome_2, gene_2, copy_2, Ext.TAIL))
                    # balancing edge:
                    G.add_edge((genome_1, gene_1, copy_1, ext_1), (genome_2, gene_2, copy_2, ext_2))


    # Find indels as connected components in A and B:
    indels = {"A": 0, "B": 0}
    for comp in connected_components(G):
        genome, gene, copy, ext = list(comp)[0]
        indels[genome] += 1


    # log file:
    gap, time = parse_log_file(filename.replace("sol", "log"))

    rearr = obj - indels["A"] - indels["B"]
    return {"dcj_distance": obj, "rearrangements": rearr, "indels_a": indels["A"], "indels_b": indels["B"],
            "time": time, "gap": gap}


def parse_nobal_ilp_sol(filename, genome_a, genome_b):
    # build the master graph with matching and adjacency edges:

    # 1st add capping genes:
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    for genome in [genome_a, genome_b]:
        for c in genome.chromosomes:
            if not c.circular:
                c.gene_order.append(0)
                c.circular = True
        for i in range(genome.n_chromosomes(), max_chromosomes):
            genome.add_chromosome(Chromosome([0], circular=True))

    master_graph = nx.Graph()
    obj = 0
    with open(filename) as f:
        for l in f:
            m = obj_regexp.match(l.strip())
            if m is not None:
                obj = int(round(float(m.group(1)), 0))
            m = matching_regexp.match(l.strip())
            if m is not None:
                if float(m.group(5)) > 0.9:
                    # add matching edge:
                    g_a, c_a, g_b, c_b = map(int, m.groups()[:4])
                    master_graph.add_edge(("A", g_a, c_a, Ext.HEAD), ("B", g_b, c_b, Ext.HEAD))
                    master_graph.add_edge(("A", g_a, c_a, Ext.TAIL), ("B", g_b, c_b, Ext.TAIL))
    # adjacency edges:
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        for (g_i, copy_a, e_i), (g_j, copy_b, e_j) in genome.adjacency_iter_with_copies():
            master_graph.add_edge((genome_name, g_i, copy_a, e_i), (genome_name, g_j, copy_b, e_j))

    # now, graph should be a collection of cycles and paths, where path extremities are balancing extremities.
    # We need to check those extremities to classify paths in AB-, AA- and BB-
    # loop for each component:
    G = nx.Graph()
    ab_components = []
    for comp in connected_components(master_graph):
        degree_one = sorted([v for v in comp if master_graph.degree(v) == 1])
        if len(degree_one) > 1:
            assert len(degree_one) == 2
            # if AA- or BB, add the balancing edge:
            if degree_one[0][0] == degree_one[1][0]:
                # edges.append(degree_one)
                G.add_edge(*degree_one)
                # add also the gene edges, from head to tail:
                # edges.extend([[[genome, gene, copy, Ext.HEAD], (genome, gene, copy, Ext.TAIL)] for genome, gene, copy, ext in degree_one])
                for genome, gene, copy, ext in degree_one:
                    G.add_edge((genome, gene, copy, Ext.HEAD), (genome, gene, copy, Ext.TAIL))
            else:
                # AB- component, save for later;
                ab_components.append(degree_one)
    # deal with AB pairing 2-by-2 arbitrarly:
    for ab_i, ab_j in zip(*[iter(ab_components)] * 2):
        G.add_edge(ab_i[0], ab_j[0])
        G.add_edge(ab_i[1], ab_j[1])
        # edges.append([ab_i[0], ab_j[0]])
        # edges.append([ab_i[1], ab_j[1]])
        # add also the gene edges, from head to tail:
        # edges.extend([[[genome, gene, copy, Ext.HEAD], (genome, gene, copy, Ext.TAIL)] for genome, gene, copy, ext in
        #               [ab_i[0], ab_i[1], ab_j[0], ab_j[1]]])
        for genome, gene, copy, ext in [ab_i[0], ab_i[1], ab_j[0], ab_j[1]]:
            G.add_edge((genome, gene, copy, Ext.HEAD), (genome, gene, copy, Ext.TAIL))

    # Find indels as connected components in A and B:
    indels = {"A": 0, "B": 0}
    for comp in connected_components(G):
        genome, gene, copy, ext = list(comp)[0]
        indels[genome] += 1

    # log file:
    gap, time = parse_log_file(filename.replace("sol", "log"))

    # correct objective function, since each AB-component is counted as a full cycle in the ILP, but
    # should only by half-cycle.
    obj += len(ab_components)/2

    rearr = obj - indels["A"] - indels["B"]
    return {"dcj_distance": obj, "rearrangements": rearr, "indels_a": indels["A"], "indels_b": indels["B"],
            "time": time, "gap": gap}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of folders for DCJdupindel simulations.")
    parser.add_argument("folders", type=str, nargs="+", help="Sim folders.")
    parser.add_argument("-out", type=str, default="dcj_summary", help="prefix filename for output.")
    parser.add_argument("-coser", action="store_true", default=False, help="Outputs also a COSER summary.")
    parser.add_argument("-sb", "--skip_balancing", action="store_true", default=False,
                        help="Uses the no balancing ILP.")
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
        genomes = file_ops.open_genomes_with_copy_number(os.path.join(folder, cfg.sim_extant_genomes()))

        # if it has balancing edges, we can find the indels by adding the bal edges and the "gene" edges
        # and each connected component is a circular chromosome in the completion.
        if not param.skip_balancing:
            sol_file = os.path.join(folder, "extant_genomes.txt_T2_T1.lp.sol")
            if not os.path.exists(sol_file):
                continue
            r = parse_ilp_sol(sol_file)
        else:
            # without balancing edges, we have to create the master graph and check which type of
            # open components we get (AA-, BB- or AB-), and then merge them into cycles; AA and BB
            # are unique, but AB- are paired arbitrarly.
            sol_file = os.path.join(folder, "extant_genomes.txt_T2_T1_nobal.lp.sol")
            if not os.path.exists(sol_file):
                continue
            r = parse_nobal_ilp_sol(sol_file, genomes[0], genomes[1])

        result.update(r)
        # result["file"] = os.path.basename(sol_file)  # .replace(".lp.sol", "")#.replace("extant_genomes.txt_", "")
        # Orthology:
        correct_assignment = build_gene_copy_assignment(genomes[0], genomes[1])
        correct, wrong = parse_assignment_quality(sol_file, correct_assignment)
        result.update({"ortho_TP": len(correct), "ortho_FP": len(wrong)})

        results[key].append(result)
        # COSER:
        if param.coser:
            if os.path.exists(os.path.join(folder, "mapping")):
                coser_result.update(parse_coser_sol(folder))
                coser_results[key].append(coser_result)

    # output:
    # DCJDUP:
    dcj_out = param.out
    if param.skip_balancing:
        dcj_out += "_nobal"
    for key, result in results.iteritems():
        with open("%s_%s.txt" % (dcj_out, key), "w") as f:
            print >> f, "\t".join(fields)
            for line in sorted(result, key=lambda r: (r['dup_length'], r["dup_prob"], r["real_distance"])):
                print >> f, "\t".join([str(line[field]) for field in fields])
    # COSER:
    if param.coser:
        for key, result in coser_results.iteritems():
            with open("coser_%s_%s.txt" % (param.out, key), "w") as f:
                print >> f, "\t".join(fields)
                for line in sorted(result, key=lambda r: (r['dup_length'], r["dup_prob"], r["real_distance"])):
                    print >> f, "\t".join([str(line[field]) for field in coser_fields])
