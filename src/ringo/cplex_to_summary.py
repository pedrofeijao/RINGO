#!/usr/bin/env python2
import collections
import os

import pyximport;

from dcj_dupindel import add_capping_genes

pyximport.install()
from model import Ext, Chromosome
import file_ops
import argparse
import re
import networkx as nx
from networkx.algorithms import connected_components
from simulation import Simulation, EventType
import ringo_config

from parse_orthology import build_correct_matching, parse_orthology_quality
from coser_to_summary import parse_coser_sol
from dcj_dupindel import build_gene_copies_dict, CopyType

time_regexp = re.compile('Solution time =\s+([\d\.]+)\s*sec.\s+Iterations')
gap_regexp = re.compile('Current MIP best bound.+gap.+\s([\d\.]+)%')
obj_regexp = re.compile("objectiveValue=\"(.+)\"")
balancing_regexp = re.compile("<variable name=\"w_([AB])(\d+)_(\d+)([ht]),([AB])(\d+)_(\d+)([ht]).*value=\"(.+)\"")
matching_regexp = re.compile("<variable name=\"x_A(\d+)_(\d+)h,B(\d+)_(\d+)h.*value=\"(.+)\"")


# simple_balancing_regexp = re.compile("<variable name=\"w_(.+?),(.+?)\".*value=\"(.+)\"")
#   <variable name="w_A121_2h,A121_2t" index="7191" value="1"/>


def parse_log_file(logfile):
    time = 0
    gap = 0
    try:
        with open(logfile) as f:
            for line in f:
                m = time_regexp.match(line.strip())
                if m is not None:
                    time = float(m.group(1))

                m = gap_regexp.match(line.strip())
                if m is not None:
                    gap = float(m.group(1))
    except Exception as e:
        print "Did not read log file!"
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

    return {"dcj_distance": obj, "dup_a": indels["A"], "dup_b": indels["B"],
            "time": time, "gap": gap}


def obj_and_matching_from_sol(filename):
    edges = []
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
                    edges.append((g_a, c_a, c_b))
    return obj, edges


def parse_nobal_ilp_sol(filename, genome_a, genome_b):
    # build the master graph with matching and adjacency edges:

    # 1st add capping genes:
    add_capping_genes(genome_a, genome_b)

    master_graph = nx.Graph()
    # matching edges:
    obj, edges = obj_and_matching_from_sol(filename)
    for gene, c_a, c_b in edges:
        master_graph.add_edge(("A", gene, c_a, Ext.HEAD), ("B", gene, c_b, Ext.HEAD))
        master_graph.add_edge(("A", gene, c_a, Ext.TAIL), ("B", gene, c_b, Ext.TAIL))

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
                G.add_edge(*degree_one)
                # add also the gene edges, from head to tail:
                for genome, gene, copy, ext in degree_one:
                    G.add_edge((genome, gene, copy, Ext.HEAD), (genome, gene, copy, Ext.TAIL))
            else:
                # AB- component, save for later;
                ab_components.append(degree_one)
    # deal with AB pairing 2-by-2 arbitrarly:
    for ab_i, ab_j in zip(*[iter(ab_components)] * 2):
        G.add_edge(ab_i[0], ab_j[0])
        G.add_edge(ab_i[1], ab_j[1])
        # add also the gene edges, from head to tail:
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
    obj += len(ab_components) / 2

    return {"dcj_distance": obj, "dup_a": indels["A"], "dup_b": indels["B"],
            "time": time, "gap": gap}


def solution_matching_ilp(sol_file):
    # get solution matching:
    sol_matching = collections.defaultdict(list)
    with open(sol_file) as f:
        for l in f:
            m = matching_regexp.match(l.strip())
            if m is not None:
                gene_a, copy_a, gene_b, copy_b, val = m.groups()
                if float(val) >= 0.9 and int(gene_a) > 0:
                    sol_matching[int(gene_a)].append((int(copy_a), int(copy_b)))
    return sol_matching


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of folders for DCJdupindel simulations.")
    parser.add_argument("folders", type=str, nargs="+", help="Sim folders.")
    parser.add_argument("-out", type=str, default="summary.txt", help="Filename for output.")
    parser.add_argument("-coser", action="store_true", default=False, help="Outputs also a COSER summary.")
    parser.add_argument("-sb", "--skip_balancing", action="store_true", default=False,
                        help="Uses the no balancing ILP.")
    param = parser.parse_args()
    # Ringo_config:
    cfg = ringo_config.RingoConfig()
    sim_cols = ["caller", "real_distance", "dup_length", "dup_prob", "del_prob"]
    tree_events = ["%s_%s" % (g, event) for g in ["T1", "T2"] for event in EventType.all]
    time_ortho = ["time", "gap", "ortho_TP", "ortho_FP", "ortho_FN"]

    fields = sim_cols + tree_events + ["dcj_distance", "distance_diff", "dup_a", "dup_b"] + time_ortho + ["bal_A", "bal_B"]

    results = collections.defaultdict(list)
    results = []
    for folder in param.folders:
        # sim data:
        sim = Simulation.open_folder(folder)
        p = sim.sim_parameters
        result = {
              "dup_length": sim.sim_parameters.duplication_length,
              "dup_prob": sim.sim_parameters.duplication_p,
              "del_prob": sim.sim_parameters.deletion_p
        }
        # key = "L%(duplication_length)s_dup%(duplication_p)s_del%(deletion_p)s" % p.__dict__
        # tree events:
        tree_events = {"%s_%s" % (g, event): sim.sim_tree.find_node_with_label(g).edge.events[event] for g in
                       ["T1", "T2"] for
                       event in EventType.all}
        result.update(tree_events)
        result["real_distance"] = sum(map(int, tree_events.values()))

        # coser and nobal is the same until here:
        coser_result = dict(result)
        no_bal_result = dict(result)

        # caller:
        result["caller"] = "ILP"
        coser_result["caller"] = "SegDCJ"
        no_bal_result["caller"] = "ILP_NB"

        # distance/rearrangement results:
        genomes = file_ops.open_genomes_with_copy_number(os.path.join(folder, cfg.sim_extant_genomes()))
        root = file_ops.open_genomes_with_copy_number(os.path.join(folder, cfg.sim_ancestral_genomes()))[0]

        # if it has balancing edges, we can find the indels by adding the bal edges and the "gene" edges
        # and each connected component is a circular chromosome in the completion.
        sol_file = os.path.join(folder, "extant_genomes.txt_T2_T1.lp.sol")
        if not os.path.exists(sol_file):
            print ("%s does not exist, skipping ..." % sol_file)
            continue
        r = parse_ilp_sol(sol_file)
        result.update(r)
        result["distance_diff"] = result["real_distance"] - result["dcj_distance"]
        if param.skip_balancing:
            # without balancing edges, we have to create the master graph and check which type of
            # open components we get (AA-, BB- or AB-), and then merge them into cycles; AA and BB
            # are unique, but AB- are paired arbitrarly.
            sol_file_no_bal = os.path.join(folder, "extant_genomes.txt_T2_T1_nobal.lp.sol")
            if not os.path.exists(sol_file_no_bal):
                print ("%s does not exist, skipping ..." % sol_file_no_bal)
                continue
            r = parse_nobal_ilp_sol(sol_file_no_bal, genomes[0], genomes[1])
            no_bal_result.update(r)
            no_bal_result["distance_diff"] = no_bal_result["real_distance"] - no_bal_result["dcj_distance"]
        # Orthology:
        # open genomes to get the correct matching:
        # correct_matching = build_correct_matching(genomes[0], genomes[1])
        # UPDATE: get parent node genes, these are the orthologs:

        correct_matching = collections.defaultdict(set)
        for chrom in root.chromosomes:
            for gene, cn in zip(chrom.gene_order, chrom.copy_number):
                correct_matching[abs(gene)].add(cn)
        correct_matching = {gene: copies for gene, copies in correct_matching.iteritems() if len(copies) > 1}
        #
        # get solution matching:
        solution_matching = solution_matching_ilp(sol_file)
        if param.skip_balancing:
            no_bal_solution_matching = solution_matching_ilp(sol_file_no_bal)

        # compare:
        tp, fp, fn = parse_orthology_quality(solution_matching, correct_matching)
        n_assignments = sum([len(x) for x in correct_matching.itervalues()])
        result.update({"ortho_TOTAL": n_assignments, "ortho_TP": len(tp), "ortho_FP": len(fp), "ortho_FN": len(fn)})

        if param.skip_balancing:
            tp, fp, fn = parse_orthology_quality(solution_matching, no_bal_solution_matching)
            n_assignments = sum([len(x) for x in correct_matching.itervalues()])
            no_bal_result.update({"ortho_TOTAL": n_assignments, "ortho_TP": len(tp), "ortho_FP": len(fp), "ortho_FN": len(fn)})


        # balancing genes:
        # since the gene set might be different for each genome, find all genes:
        all_genes = genomes[0].gene_set().union(genomes[1].gene_set())
        # find all gene copies
        gene_copies = build_gene_copies_dict(all_genes, genomes[0], genomes[1])
        # count balancing genes:
        bal = {
            g: sum(
                [len([c for c in gene_copies[g][gene].itervalues() if c == CopyType.BALANCING]) for gene in all_genes])
            for g in ["A", "B"]}

        result["bal_A"] = bal["A"]
        result["bal_B"] = bal["B"]
        if param.skip_balancing:
            no_bal_result["bal_A"] = bal["A"]
            no_bal_result["bal_B"] = bal["B"]
        if param.coser:
            coser_result["bal_A"] = "NA"
            coser_result["bal_B"] = "NA"

        # save keep result:
        # results[key].append(result)
        results.append(result)
        if param.skip_balancing:
            results.append(no_bal_result)

        # COSER:
        if param.coser:
            if os.path.exists(os.path.join(folder, "mapping")):
                coser_result.update(parse_coser_sol(folder, correct_matching))
                coser_result["ortho_TOTAL"] = n_assignments
                coser_result["distance_diff"] = coser_result["real_distance"] - coser_result["dcj_distance"]
                results.append(coser_result)

    # output:
    with open(param.out, "w") as f:
        print >> f, "\t".join(fields)
        # ILP
        for result in sorted(results, key=lambda r: (r["dup_length"], r["dup_prob"])):
            print >> f, "\t".join([str(result[field]) for field in fields])
