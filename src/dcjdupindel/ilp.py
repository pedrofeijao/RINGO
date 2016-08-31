#!/usr/bin/env python2
import argparse
import copy

import pyximport;
from math import factorial

pyximport.install()
import os
import sys
import itertools
import operator
from model import Genome, Chromosome, Ext
import networkx as nx
from networkx.algorithms import connected_components

# TODO: this is a hack to import from other directory; should use packages
ringo_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../ringo"))
os.sys.path.insert(0, ringo_path)
import file_ops


# HELPER functions:
def vertex_name(genome, gene, copy, ext):
    return "%s%s_%s%s" % (genome, gene, copy, ext)


def matching_edge_name(gene, copyA, copyB, ext):
    return "x_%s,%s" % (vertex_name("A", gene, copyA, ext), vertex_name("B", gene, copyB, ext))


def balancing_edge_name(genome, gene1, copy1, ext1, gene2, copy2, ext2):
    if (gene1, copy1, ext1) > (gene2, copy2, ext2):
        gene1, copy1, ext1, gene2, copy2, ext2 = gene2, copy2, ext2, gene1, copy1, ext1
    return "w_%s,%s" % (vertex_name(genome, gene1, copy1, ext1), vertex_name(genome, gene2, copy2, ext2))


def self_edge_name(genome, gene, copy):
    # return "x_%s,%s" % (vertex_name(genome, gene, copy, Ext.HEAD), vertex_name(genome, gene, copy, Ext.TAIL))
    return "x_%s" % (vertex_name(genome, gene, copy, Ext.HEAD) + Ext.TAIL)


def balancing_extremities(balancing, exclude=None):
    if exclude is None:
        exclude = set()
    for gene, copies in balancing.iteritems():
        for copy in copies:
            if (gene, copy, Ext.HEAD) not in exclude:
                yield gene, copy, Ext.HEAD
            if (gene, copy, Ext.TAIL) not in exclude:
                yield gene, copy, Ext.TAIL


def build_extremity_order_lists(genome, gene_count, max_chromosomes):
    copy_number = {gene: 1 for gene in gene_count.iterkeys()}
    copy_number[0] = 1
    chromosomes = genome.chromosomes
    ext_order_list = []
    # if genomes have diff number of chromosomes, complete with empty to balance:
    if genome.n_chromosomes() < max_chromosomes:
        chromosomes += [Chromosome([])] * (max_chromosomes - genome.n_chromosomes())
    for idx, chrom in enumerate(chromosomes):
        ext_order_list.append(build_chromosome_ext_order(copy_number, chrom))
    return ext_order_list


def build_chromosome_ext_order(copy_number, chrom):
    # returns a list of tuplets (gene, copy, extremity) for a given chromosome "chrom",
    # including telomeres (0-labelled genes).

    ext_list = [(0, copy_number[0], Ext.TAIL)]
    copy_number[0] += 1
    for gene in chrom.gene_order:
        if gene > 0:
            orientation = [Ext.TAIL, Ext.HEAD]
        else:
            orientation = [Ext.HEAD, Ext.TAIL]
        ext_list.extend([(abs(gene), copy_number[abs(gene)], ext) for ext in orientation])
        copy_number[abs(gene)] += 1
    ext_list.append((0, copy_number[0], Ext.TAIL))
    copy_number[0] += 1
    return ext_list


def adjacency_list(genome, gene_count, max_chromosomes):
    # using the tuplet list generated from "build_extremity_order_lists", outputs a
    # list of pairs of tuplets, represeting the adjacencies.
    ext_order_list = build_extremity_order_lists(genome, gene_count, max_chromosomes)
    for ext_order in ext_order_list:
        a = iter(ext_order)
        return itertools.izip(a, a)


def fix_conserved_adjacencies(edges, genome_a, genome_b, total_gene_count, gene_count, max_chromosomes):
    # Matching Edge FIXING:
    # conserved adjacencies (a,b) where one of the families is fixed (either singleton or fixed in any other way)
    adj_dict = {"A": {}, "B": {}}
    for genome, name in [(genome_a, "A"), (genome_b, "B")]:
        # TODO: maybe adj_dict should point to a LIST of all matching adjacencies, otherwise it is storing
        # just the last time it sees it.
        for (gene_i, copy_i, ext_i), (gene_j, copy_j, ext_j) in adjacency_list(genome, total_gene_count,
                                                                               max_chromosomes):
            adj_i = gene_i if ext_i == Ext.HEAD else -gene_i
            adj_j = gene_j if ext_j == Ext.TAIL else -gene_j
            if adj_i < 0:
                adj_i, adj_j = -adj_j, -adj_i
            adj_dict[name][(adj_i, adj_j)] = (gene_i, copy_i, ext_i), (gene_j, copy_j, ext_j)

    # If some variable is fixed, a new round is called so it might get something new;
    new_round = True
    while new_round:
        new_round = False
        for adj in sorted(adj_dict["A"].iterkeys()):
            (gene_a_i, gene_a_j) = adj_dict["A"][adj]
            # if this adjacency is conserved in B:
            # if gene_a_i[0] == 638: #gene_a_i == (638, 2, Ext.TAIL):
            #     import ipdb;
            #     ipdb.set_trace()
            # if gene_a_j[0] == 638: #    (638, 2, Ext.TAIL):
            #     import ipdb;
            #     ipdb.set_trace()

            if adj in adj_dict["B"]:
                gene_b_i, gene_b_j = adj_dict["B"][adj]
                # reverse adjacency to match genes if necessary:
                if (gene_a_i[0], gene_a_i[2]) != (gene_b_i[0], gene_b_i[2]):
                    gene_b_i, gene_b_j = gene_b_j, gene_b_i

                # DEBUG:

                gene_i, copy_a_i = gene_a_i[0], gene_a_i[1]
                gene_j, copy_a_j = gene_a_j[0], gene_a_j[1]
                copy_b_i = gene_b_i[1]
                copy_b_j = gene_b_j[1]

                pair_a_i = (gene_i, copy_a_i)  # gene, copy (without tail/head)
                pair_a_j = (gene_j, copy_a_j)  # gene, copy (without tail/head)
                # if (a_i,b_i) are fixed, we can try to fix the a_j,b_j:
                if edges[pair_a_i] == {copy_b_i}:

                    # if not fixed yet, fix:
                    if copy_b_j in edges[pair_a_j] and len(edges[pair_a_j]) > 1:
                        edges[pair_a_j] = {copy_b_j}
                        # # remove this edge from other copies:
                        for idx in xrange(1, total_gene_count[gene_a_j[0]] + 1):
                            if idx == copy_a_j:
                                continue
                            try:
                                # if not there already, exception is thrown, that' ok
                                edges[(gene_j, idx)].remove(copy_b_j)
                            except KeyError:
                                pass

                        new_round = True
                elif edges[pair_a_j] == {copy_b_j}:
                    if copy_b_i in edges[pair_a_i] and len(edges[pair_a_i]) > 1:
                        edges[pair_a_i] = {copy_b_i}
                        # # remove this edge from other copies:
                        for idx in xrange(1, total_gene_count[gene_a_i[0]] + 1):
                            if idx == copy_a_i:
                                continue
                            try:
                                # if not there already, exception is thrown, that' ok
                                edges[(gene_i, idx)].remove(copy_b_i)
                            except KeyError:
                                pass
                        new_round = True


def define_y_label(gene_count, n_telomeres):
    y_label = {}
    idx = 1
    for genome in ["A", "B"]:
        # telomeres: (represented by TAIL always)
        for i in range(1, n_telomeres + 1):
            y_label[vertex_name(genome, 0, i, Ext.TAIL)] = idx
            idx += 1
        for gene, copies in gene_count.iteritems():
            for i in xrange(1, copies + 1):
                y_label[vertex_name(genome, gene, i, Ext.HEAD)] = idx
                idx += 1
                y_label[vertex_name(genome, gene, i, Ext.TAIL)] = idx
                idx += 1
    return y_label


####################################################################
# ILP
####################################################################
# TODO: I'm assuming the genomes are linear, and then treating telomeres;
# should just circularize with 0 genes;


def dcj_dupindel_ilp(genome_a, genome_b, output):
    # some useful vars:
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    n_telomeres = 2 * max_chromosomes

    # count of each gene on each genome
    gene_count = {"A": genome_a.gene_count(), "B": genome_b.gene_count()}
    gene_count["A"][0] = n_telomeres
    gene_count["B"][0] = n_telomeres
    # for all genes ,the "balanced" count:
    total_gene_count = {g: max(gene_count["A"][g], gene_count["B"][g]) for g in
                        set(gene_count["A"].keys()).union(set(gene_count["B"].keys()))}

    # define the y labels -> integer 1..n
    y_label = define_y_label(total_gene_count, n_telomeres)

    # list of possible edges for each vertex:
    edges = {}
    for gene, copies in total_gene_count.iteritems():
        for i in xrange(1, copies + 1):
            edges[(gene, i)] = set(range(1, copies + 1))

    # telomeres:
    for i in xrange(1, n_telomeres + 1):
        edges[(0, i)] = set(range(1, n_telomeres + 1))

    # Fix by conserved adjacencies:
    # TODO: check the edge removal; I might be creating useless matching variables that I know are 0.
    # not on the unique matching constraints, but maybe on the y_i.
    fix_conserved_adjacencies(edges, genome_a, genome_b, total_gene_count, gene_count, max_chromosomes)

    # Build the BP graph of fixed elements to try to find more variables to fix:
    master_graph = nx.Graph()
    y_fix = {}
    for (gene, copy_i), set_y in edges.iteritems():
        if len(set_y) == 1:
            copy_j = list(set_y)[0]
            for ext in [Ext.HEAD, Ext.TAIL]:
                master_graph.add_edge(("A", gene, copy_i, ext), ("B", gene, copy_j, ext))
    #
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        # for i, j in adjacency_list(genome, gene_count, max_chromosomes):
        for (g_i, c_i, e_i), (g_j, c_j, e_j) in adjacency_list(genome, total_gene_count, max_chromosomes):
            master_graph.add_edge((genome_name, g_i, c_i, e_i), (genome_name, g_j, c_j, e_j))

    AB_components = []
    balancing_fix = {"A": {}, "B": {}}
    for comp in connected_components(master_graph):
        # cycles have all degree 2
        # print " ".join(comp)
        degree_one = [v for v in comp if master_graph.degree(v) == 1]

        # if two degree one vertices, it is a path;
        if len(degree_one) == 2:
            genome_i, g_i, c_i, e_i = degree_one[0]
            genome_j, g_j, c_j, e_j = degree_one[1]
            # 1 - check if nodes are balancing, to find AA-, BB- and AB- paths that can be fixed.
            i_is_balancing = g_i != 0 and c_i > gene_count[genome_i][g_i]
            j_is_balancing = g_j != 0 and c_j > gene_count[genome_j][g_j]
            if i_is_balancing and j_is_balancing:
                if genome_i == genome_j:  # AA- or BB-path, close it
                    balancing_fix[genome_i][degree_one[0][1:]] = degree_one[1][1:]
                    balancing_fix[genome_i][degree_one[1][1:]] = degree_one[0][1:]
                    degree_one = []
                else:
                    # TODO: deal with AB-components;
                    pass

            # TODO: investigate components where both degree-1 vertices have the same gene on diff genomes;
            # could it be that I can close them? I tried quickly and I got a larger ILP and run time...
            # UPDATE:If I just close the cycle, I get a higher value; so apparently not optimal. If I
            # also fix the y_i on this cycle, I get an infeasible model, which is very strange. I have
            # to check this better. I think it would even be possible to do even the adjacency fix with
            # the graph, by closing components, making it all simpler in the end.
            # UPDATE2: aftex improving the edge fixing, this is working now!
            # Which makes me think again that it can be from the start.

            elif genome_i != genome_j and g_i == g_j: # and len(comp) == 4:
                if genome_j == "A":
                    genome_i, g_i, c_i, e_i, genome_j, g_j, c_j, e_j = g_j, c_j, e_j, genome_i, g_i, c_i, e_i, genome_j

                    degree_one = []
                    # print "SAME deg1:",
                    # print [(v, master_graph.degree(v)) for v in sorted(comp)]
                    # for g, gen, c, e in comp:
                        # if g == "A":
                        #     print "EDGE:", gen,c, "->", edges[(gen, c)]
                    # import ipdb; ipdb.set_trace()
                    edges[(g_i, c_i)] = {c_j}
                    for idx in xrange(1, total_gene_count[g_i] + 1):
                        if idx == c_i:
                            continue
                        try:
                            # if not there already, exception is thrown, that' ok
                            edges[(g_i, idx)].remove(c_j)
                        except KeyError:
                            pass

        # if no degree one vertices, it is a cycle, I can fix the y_i:
        if len(degree_one) == 0:
            # get indexes of the y_i:
            indexes = [(v, y_label[vertex_name(*v)]) for v in comp]
            min_label = min([x[1] for x in indexes])
            for v, label in indexes:
                y_fix[label] = min_label

    # for x,y in sorted(y_fix.items(), key=operator.itemgetter(1)):
    #     print x,y

    # sys.exit()

    # build constraints:
    constraints = []

    # fix y:
    for i, label in y_fix.iteritems():
        constraints.append("y_%d = %d" % (i, label))
    # consistency:
    constraints.append("\ Consistency constraints")
    for gene, copies in total_gene_count.iteritems():
        constraints.append("\\ gene %s" % gene)
        for i, j in itertools.product(xrange(1, copies + 1), xrange(1, copies + 1)):
            if edges[(gene, i)] != {j}:
                constraints.append(
                    "%s - %s = 0" % (
                        matching_edge_name(gene, i, j, Ext.TAIL), matching_edge_name(gene, i, j, Ext.HEAD)))

    # unique matching: (1-to-1) (only tail, it should suffice, together with the consistency)
    constraints.append("\ Unique matching:")
    for gene, copies in total_gene_count.iteritems():
        # degree 1 genes are already fixed:
        if copies == 1:
            continue
        for i in xrange(1, copies + 1):
            constraints.append(
                " + ".join([matching_edge_name(gene, i, j, Ext.TAIL) for j in xrange(1, copies + 1)]) + " = 1")

    constraints.append("\ Telomere matching:")  # special case, "0" gene.
    for i in range(1, n_telomeres + 1):
        constraints.append(
            " + ".join([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)]) + " = 1")

    constraints.append("\ Balancing:")
    balancing_A = {g: range(gene_count["A"][g] + 1, gene_count["B"][g] + 1) for g in total_gene_count.iterkeys() if
                   gene_count["A"][g] < gene_count["B"][g]}
    balancing_B = {g: range(gene_count["B"][g] + 1, gene_count["A"][g] + 1) for g in total_gene_count.iterkeys() if
                   gene_count["B"][g] < gene_count["A"][g]}
    # import ipdb; ipdb.set_trace()
    for genome, balancing in [("A", balancing_A), ("B", balancing_B)]:
        constraints.append("\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing):
            if (gene_i, copy_i, ext_i) in balancing_fix[genome]:
                gene_j, copy_j, ext_j = balancing_fix[genome][(gene_i, copy_i, ext_i)]
                if (gene_i, copy_i, ext_i) < (gene_j, copy_j, ext_j):
                    constraints.append(
                        "%s = 1" % balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j))
            else:
                constraints.append(
                    " + ".join([balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j) for
                                gene_j, copy_j, ext_j in
                                balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) if
                                (gene_i, copy_i, ext_i) != (gene_j, copy_j, ext_j)]) + " = 1")

    constraints.append("\ Labelling")

    # for each adjacency, fix label:
    constraints.append("\\ Adjacency have the same label:")

    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        for i, j in adjacency_list(genome, total_gene_count, max_chromosomes):
            v_i = vertex_name(genome_name, *i)
            v_j = vertex_name(genome_name, *j)
            # if already fixed, skip
            if y_label[v_i] in y_fix and y_label[v_j] in y_fix:
                continue
            constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (y_label[v_i], y_label[v_j], v_i, v_j))

    constraints.append("\\ Matching edges with the same label:")
    for gene, copies in total_gene_count.iteritems():
        for i in xrange(1, copies + 1):
            for j in xrange(1, copies + 1):
                constraints.append("\ Edges (%s%s_%s, %s%s_%s)" % ("A", gene, i, "B", gene, j))
                for ext in [Ext.HEAD, Ext.TAIL]:
                    y_i = y_label[vertex_name("A", gene, i, ext)]
                    y_j = y_label[vertex_name("B", gene, j, ext)]
                    # if already fixed, skip
                    if y_i in y_fix and y_j in y_fix:
                        continue
                    if edges[(gene, i)] == {j}:
                        constraints.append("y_%s - y_%s = 0 " % (y_i, y_j))
                    else:
                        constraints.append(
                            "y_%s - y_%s + %s %s <= %d" % (y_i, y_j, y_i, matching_edge_name(gene, i, j, ext), y_i))
                        constraints.append(
                            "y_%s - y_%s + %s %s <= %d" % (y_j, y_i, y_j, matching_edge_name(gene, i, j, ext), y_j))

    constraints.append("\\ Balancing edges with same label:")
    for genome, balancing in [("A", balancing_A), ("B", balancing_B)]:
        constraints.append("\\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()):
            for gene_j, copy_j, ext_j in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()):
                if (gene_i, copy_i, ext_i) == (gene_j, copy_j, ext_j):
                    continue
                y_i = y_label[vertex_name(genome, gene_i, copy_i, ext_i)]
                y_j = y_label[vertex_name(genome, gene_j, copy_j, ext_j)]
                # should not have someone here if I'm excluding fixed edges:
                if y_i in y_fix and y_j in y_fix:
                    continue
                constraints.append("y_%s - y_%s + %s %s <= %d" % (
                    y_i, y_j, y_i, balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j), y_i))
                constraints.append("y_%s - y_%s + %s %s <= %d" % (
                    y_j, y_i, y_j, balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j), y_j))

    # z variables: since all cycles have to contains vertices from both genomes, we only add z variables
    # for genome A, that have smallest labels, so a genome B z variable will never be =1.
    constraints.append("\\ Z variables")
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "A":
            constraints.append("%d z_%s - y_%s <= 0" % (i, i, i))

    # number of genes, to fix distance:
    constraints.append("n = %d" % (sum(total_gene_count.itervalues()) + n_telomeres / 2))

    # bounds:
    bounds = []
    for i in sorted(y_label.itervalues()):
        bounds.append("y_%d <= %d" % (i, i))

    # variables:
    binary = []

    # matching edges

    # telomere:
    for i in range(1, n_telomeres + 1):
        binary.extend([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)])

    # matching edges, skipping fixed pairs.
    matching = [matching_edge_name(gene, i, j, ext) for gene, copies in total_gene_count.iteritems() for i in
                xrange(1, copies + 1)
                for j in xrange(1, copies + 1) if edges[(gene, i)] != {j} for ext in [Ext.HEAD, Ext.TAIL]]
    print "%d matching edges" % len(matching)
    # print "Potentially %d matching edges" % sum([2*x ** 2 for x in gene_count.itervalues()])
    binary.extend(matching)

    # balancing edges:
    balancing_edges = [balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j) for genome, balancing
                       in [("A", balancing_A), ("B", balancing_B)] for gene_i, copy_i, ext_i in
                       balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) for gene_j, copy_j, ext_j
                       in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) if
                       (gene_i, copy_i, ext_i) < (gene_j, copy_j, ext_j)]
    print "%d balancing edges" % len(balancing_edges)
    binary.extend(balancing_edges)

    # z cycles:
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "B":
            continue
        binary.append("z_%d" % i)

    # Y label are general:
    general = []
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        general.append("y_%d" % i)

    # number of genes:
    general.append("n")
    # objective function:
    objective = ["obj: n - " + " - ".join(
        ["z_%d" % i for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)) if vertex[0] != "B"])]

    # write:
    with open(output, "w") as f:
        for header, lines in [("Minimize", objective), ("Subject to", constraints),
                              ("Bounds", bounds), ("Binary", binary), ("General", general)]:
            print >> f, header
            print >> f, "\n".join(lines)


def dcj_dupindel_ilp_unitary(genome_a, genome_b, output):
    # count of each gene on each genome
    gene_count_a = genome_a.gene_count()
    gene_count_b = genome_b.gene_count()

    # find all genes
    gene_count = {g: max(gene_count_a[g], gene_count_b[g]) for g in
                  set(gene_count_a.keys()).union(set(gene_count_b.keys()))}

    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    n_telomeres = 2 * max_chromosomes
    y_label = define_y_label(gene_count, n_telomeres)
    # build constraints:
    constraints = []
    # consistency:
    constraints.append("\ Consistency constraints")
    for gene, copies in gene_count.iteritems():
        constraints.append("\\ gene %s" % gene)
        for i, j in itertools.product(xrange(1, gene_count_a[gene] + 1), xrange(1, gene_count_b[gene] + 1)):
            constraints.append(
                "%s - %s = 0" % (matching_edge_name(gene, i, j, Ext.TAIL), matching_edge_name(gene, i, j, Ext.HEAD)))

    # unique matching: (1-to-1) (only tail, it should suffice, together with the consistency)
    constraints.append("\\ Unique matching on each vertex (with self edges):")
    for genome, count_1, count_2 in [("A", gene_count_a, gene_count_b), ("B", gene_count_b, gene_count_a)]:
        constraints.append("\ Genome %s" % genome)
        for gene, copies in gene_count.iteritems():
            for i in xrange(1, count_1[gene] + 1):
                # for ext in [Ext.HEAD, Ext.TAIL]:
                for ext in [Ext.TAIL]:
                    c = " + ".join([matching_edge_name(gene, i, j, ext) for j in xrange(1, count_2[gene] + 1)])
                    # if there are more copies on this genome, add the edge:
                    if count_1[gene] > count_2[gene]:
                        c += "%s%s" % (" + " if len(c) > 0 else "", self_edge_name(genome, gene, i))
                    c += " = 1 \ %s" % vertex_name(genome, gene, i, ext)
                    # print "C:", c
                    constraints.append(c)

    ######
    constraints.append("\\ Telomere matching:")  # special case, "0" gene.

    for i in range(1, n_telomeres + 1):
        constraints.append(
            " + ".join([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)]) + " = 1")

    constraints.append("\\ Labelling")

    # for each adjacency, fix label:
    constraints.append("\\ Adjacency have the same label:")
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        copy_number = {gene: 1 for gene in gene_count.iterkeys()}
        copy_number[0] = 1
        chromosomes = genome.chromosomes
        # if genomes have diff number of chromosomes, complete with empty to balance:
        if genome.n_chromosomes() < max_chromosomes:
            chromosomes += [Chromosome([])] * (max_chromosomes - genome.n_chromosomes())
        for idx, chrom in enumerate(chromosomes):
            adj = [(0, copy_number[0], Ext.TAIL)]
            copy_number[0] += 1
            for gene in chrom.gene_order:
                if gene > 0:
                    orientation = [Ext.TAIL, Ext.HEAD]
                else:
                    orientation = [Ext.HEAD, Ext.TAIL]
                adj.extend([(abs(gene), copy_number[abs(gene)], ext) for ext in orientation])
                copy_number[abs(gene)] += 1
            adj.append((0, copy_number[0], Ext.TAIL))
            copy_number[0] += 1

            a = iter(range(len(adj)))
            for i, j in itertools.izip(a, a):
                v_i = vertex_name(genome_name, *adj[i])
                v_j = vertex_name(genome_name, *adj[j])
                constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (y_label[v_i], y_label[v_j], v_i, v_j))

    constraints.append("\\ Matching edges with the same label:")
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, gene_count_a[gene] + 1):
            for j in xrange(1, gene_count_b[gene] + 1):
                constraints.append("\\ Edges (%s%s_%s, %s%s_%s)" % ("A", gene, i, "B", gene, j))
                for ext in [Ext.HEAD, Ext.TAIL]:
                    y_i = y_label[vertex_name("A", gene, i, ext)]
                    y_j = y_label[vertex_name("B", gene, j, ext)]
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_i, y_j, y_i, matching_edge_name(gene, i, j, ext), y_i))
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_j, y_i, y_j, matching_edge_name(gene, i, j, ext), y_j))

    # self edges with same label:
    constraints.append("\ Self edges, same label")

    for genome, count_1, count_2 in [("A", gene_count_a, gene_count_b), ("B", gene_count_b, gene_count_a)]:
        constraints.append("\ Genome %s" % genome)
        for gene, copies in gene_count.iteritems():
            if count_1[gene] > count_2[gene]:
                for i in range(1, count_1[gene] + 1):
                    y_i = y_label[vertex_name(genome, gene, i, Ext.HEAD)]
                    y_j = y_label[vertex_name(genome, gene, i, Ext.TAIL)]
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_i, y_j, y_i, self_edge_name(genome, gene, i), y_i))
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_j, y_i, y_j, self_edge_name(genome, gene, i), y_j))

    # z variables: since all cycles have to contains vertices from both genomes, we only add z variables
    # for genome A, that have smallest labels, so a genome B z variable will never be =1.
    constraints.append("\\ Z variables")
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "A":
            constraints.append("%d z_%s - y_%s <= 0" % (i, i, i))

    # number of genes, to fix distance:
    constraints.append("n = %d" % (sum(gene_count.itervalues()) + n_telomeres / 2))

    # bounds:
    bounds = []
    for i in sorted(y_label.itervalues()):
        bounds.append("y_%d <= %d" % (i, i))

    # variables:
    binary = []

    # matching variables:
    # telomere matching:
    for i in range(1, n_telomeres + 1):
        binary.extend([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)])

    # matching genes:
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, gene_count_a[gene] + 1):
            for j in xrange(1, gene_count_b[gene] + 1):
                for ext in [Ext.HEAD, Ext.TAIL]:
                    binary.append(matching_edge_name(gene, i, j, ext))
    # self edges:
    self_edges = []
    for gene, copies in gene_count.iteritems():
        for genome, count_1, count_2 in [("A", gene_count_a, gene_count_b), ("B", gene_count_b, gene_count_a)]:
            if count_1[gene] > count_2[gene]:
                for i in range(1, count_1[gene] + 1):
                    binary.append(self_edge_name(genome, gene, i))
                    self_edges.append(self_edge_name(genome, gene, i))

    # z cycles:
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "B":
            continue
        binary.append("z_%d" % i)

    # Y label are general:
    general = []
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        general.append("y_%d" % i)

    # number of genes:
    general.append("n")
    # objective function:
    objective = ["obj: n - " + " - ".join(
        ["z_%d" % i for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)) if vertex[0] != "B"])]

    # write:
    with open(output, "w") as f:
        for header, lines in [("Minimize", objective), ("Subject to", constraints), ("Bounds", bounds),
                              ("Binary", binary), ("General", general)]:
            print >> f, header
            print >> f, "\n".join(lines)


def solve_ilp(filename, timelimit=60):
    # import here, so only if actually solving we will need gurobi.
    from gurobipy import read, GRB
    # pycharm complains of gurobi commands, cannot see them from the import
    model = read(filename)

    # set some options:
    # time limit in seconds:
    model.params.timeLimit = timelimit

    model.optimize()

    if model.status != GRB.Status.INFEASIBLE:
        print('FINISHED: Best objective: %g' % model.objVal)
        print('Optimization ended with status %d' % model.status)
        model.write(filename + '.sol')

    return model


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generates and optionally solve an ILP for the DCJ duplication and indel distance.")
    parser.add_argument("-s", "--solve", action="store_true", default=False, help="Solve the model with Gurobi.")
    parser.add_argument("-u", "--unitary", action="store_true", default=False, help="Use the unitary model.")
    parser.add_argument("file", type=str, help="Genomes file.")
    parser.add_argument("g1", type=int, help="Index of genome 1 in the file. (0-indexed).")
    parser.add_argument("g2", type=int, help="Index of genome 1 in the file. (0-indexed).")
    param = parser.parse_args()

    genomes = file_ops.open_genome_file(param.file, as_list=True)
    filename = "%s_%d_%d%s.lp" % (os.path.basename(param.file), param.g1, param.g2, "_unitary" if param.unitary else "")
    if param.unitary:
        dcj_dupindel_ilp_unitary(genomes[param.g1], genomes[param.g2], filename)
    else:
        dcj_dupindel_ilp(genomes[param.g1], genomes[param.g2], filename)

    if param.solve:
        model = solve_ilp(filename)
