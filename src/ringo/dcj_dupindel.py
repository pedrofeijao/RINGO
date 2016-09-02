#!/usr/bin/env python2
import pyximport;pyximport.install()
from model import Ext, Chromosome
import argparse
import copy
import os
import itertools
import operator
import networkx as nx
from networkx.algorithms import connected_components
import matplotlib.pyplot as plt
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


def define_y_label(gene_count):
    y_label = {}
    idx = 1
    for genome in ["A", "B"]:
        for gene, copies in gene_count.iteritems():
            for i in xrange(1, copies + 1):
                y_label[vertex_name(genome, gene, i, Ext.HEAD)] = idx
                idx += 1
                y_label[vertex_name(genome, gene, i, Ext.TAIL)] = idx
                idx += 1
    return y_label


def build_extremity_order_lists(genome, gene_count):
    copy_number = {gene: 1 for gene in gene_count.iterkeys()}
    ext_order_list = []
    for idx, chrom in enumerate(genome.chromosomes):
        ext_order_list.append(build_chromosome_ext_order(copy_number, chrom))
    return ext_order_list


def build_chromosome_ext_order(copy_number, chromosome):
    # returns a list of tuplets (gene, copy, extremity) for the extremities of a given chromosome "chrom"
    ext_list = []
    for gene in chromosome.gene_order:
        if gene >= 0:
            orientation = [Ext.TAIL, Ext.HEAD]
        else:
            orientation = [Ext.HEAD, Ext.TAIL]
        ext_list.extend([(abs(gene), copy_number[abs(gene)], ext) for ext in orientation])
        copy_number[abs(gene)] += 1
    return ext_list


def adjacency_list(genome, gene_count):
    # using the tuplet list generated from "build_extremity_order_lists", outputs a
    # list of pairs of tuplets, represeting the adjacencies.
    ext_order_list = build_extremity_order_lists(genome, gene_count)
    for ext_order in ext_order_list:
        # rotate 1 to make the adjacencies:
        a = iter(ext_order[1:] + ext_order[:1])
        # yield
        for i, j in itertools.izip(a, a):
            yield i, j


def balancing_extremities(balancing, exclude=None):
    if exclude is None:
        exclude = set()
    for gene, copies in balancing.iteritems():
        for copy in copies:
            if (gene, copy, Ext.HEAD) not in exclude:
                yield gene, copy, Ext.HEAD
            if (gene, copy, Ext.TAIL) not in exclude:
                yield gene, copy, Ext.TAIL


######################################################################
# MAIN function
######################################################################

def dcj_dupindel_ilp(genome_a, genome_b, output):
    # copy genomes to possibly make some changes:
    genome_a = copy.deepcopy(genome_a)
    genome_b = copy.deepcopy(genome_b)
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())

    # add capping genes:
    for genome in [genome_a, genome_b]:
        for c in genome.chromosomes:
            if not c.circular:
                c.gene_order.append(0)
                c.circular = True
        for i in range(genome.n_chromosomes(), max_chromosomes):
            genome.add_chromosome(Chromosome([0], circular=True))

    # count of each gene on each genome
    gene_count = {"A": genome_a.gene_count(), "B": genome_b.gene_count()}

    # for all genes ,the total "balanced" count:
    total_gene_count = {g: max(gene_count["A"][g], gene_count["B"][g]) for g in
                        set(gene_count["A"].keys()).union(set(gene_count["B"].keys()))}

    # define the y labels -> integer 1..n
    y_label = define_y_label(total_gene_count)

    # list of possible edges for each vertex:
    edges = {}
    for gene, copies in total_gene_count.iteritems():
        for i in xrange(1, copies + 1):
            edges[(gene, i)] = set(range(1, copies + 1))

    # try to fix variables:
    # Build the BP graph of fixed elements to try to find more variables to fix:
    master_graph = nx.Graph()
    # fixed vars:
    y_fix = {}
    z_fix = {}
    balancing_fix = {"A": {}, "B": {}}

    # add matching edges of genes with single copy:
    for (gene, copy_a), set_y in edges.iteritems():
        if len(set_y) == 1:
            copy_b = list(set_y)[0]
            for ext in [Ext.HEAD, Ext.TAIL]:
                master_graph.add_edge(("A", gene, copy_a, ext), ("B", gene, copy_b, ext))

    # add adjacency edges:
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        for (g_i, copy_a, e_i), (g_j, copy_b, e_j) in adjacency_list(genome, total_gene_count):
            master_graph.add_edge((genome_name, g_i, copy_a, e_i), (genome_name, g_j, copy_b, e_j))

    # Search components to fix:
    rescan = True
    edges_to_add = []
    vertices_to_remove = []
    while rescan:
        rescan = False
        master_graph.add_edges_from(edges_to_add)
        master_graph.remove_nodes_from(vertices_to_remove)
        edges_to_add = []
        vertices_to_remove = []

        # check each connected component:
        for comp in connected_components(master_graph):
            # get degree-1 vertices:
            degree_one = [v for v in comp if master_graph.degree(v) == 1]

            # if two degree one vertices, it is a path;
            if len(degree_one) == 2:
                genome_i, g_i, copy_a, e_i = degree_one[0]
                genome_j, g_j, copy_b, e_j = degree_one[1]
                # 1 - check if nodes are balancing, to find AA-, BB- and AB- paths that can be fixed.
                i_is_balancing = g_i != 0 and copy_a > gene_count[genome_i][g_i]
                j_is_balancing = g_j != 0 and copy_b > gene_count[genome_j][g_j]
                if i_is_balancing and j_is_balancing:
                    if genome_i == genome_j:  # AA- or BB-path, close it
                        balancing_fix[genome_i][degree_one[0][1:]] = degree_one[1][1:]
                        balancing_fix[genome_i][degree_one[1][1:]] = degree_one[0][1:]
                    else:
                        # TODO: deal with AB-components;
                        pass

                # if the path has homologous genes at the ends, I can join:
                elif genome_i != genome_j and g_i == g_j:
                    # invert to put genome A always in variables _i :
                    if genome_j == "A":
                        genome_i, g_i, copy_a, e_i, genome_j, g_j, copy_b, e_j = genome_j, g_j, copy_b, e_j, genome_i, g_i, copy_a, e_i
                    # check conflict, only add edge if ok:
                    if copy_b in edges[(g_i, copy_a)]:
                        edges[(g_i, copy_a)] = {copy_b}
                        # save edges to add to graph:
                        for ext in [Ext.HEAD, Ext.TAIL]:
                            edges_to_add.append((("A", g_i, copy_a, ext), ("B", g_i, copy_b, ext)))
                        # new edges, re-scan:
                        rescan = True

                        # remove possible edges from other copies:
                        for idx in xrange(1, total_gene_count[g_i] + 1):
                            if idx == copy_a:
                                continue
                            try:
                                # if not there already, exception is thrown, that' ok
                                edges[(g_i, idx)].remove(copy_b)
                                # Add new edges to graph, if the removal created degree 1 vertices:
                                if len(edges[(g_i, idx)]) == 1:
                                    idx_c = list(edges[(g_i, idx)])[0]
                                    for ext in [Ext.HEAD, Ext.TAIL]:
                                        edges_to_add.append((("A", g_i, idx, ext), ("B", g_i, idx_c, ext)))
                            except KeyError:
                                pass
            # if no degree one vertices, it is a cycle, I can fix the y_i:
            elif len(degree_one) == 0:
                # get indexes of the y_i:
                indexes = [(v, y_label[vertex_name(*v)]) for v in comp]
                min_label = min([x[1] for x in indexes])
                for v, label in indexes:
                    y_fix[label] = min_label
                    z_fix[label] = 0
                z_fix[min_label] = 1
                vertices_to_remove.extend(comp)

    # DRAW?
    # nx.draw_circular(master_graph, font_size=8, width=0.5, node_shape="8", node_size=1, with_labels=True)
    # nx.draw_spring(master_graph, font_size=8, width=0.5, node_shape="8", node_size=20, with_labels=True)
    # nx.draw_spectral(master_graph, font_size=8, width=0.5, node_shape="8", node_size=20, with_labels=True)
    # nx.draw_graphviz(master_graph, font_size=8, width=0.5, node_shape="8", node_size=20, with_labels=True)
    # plt.savefig('graph.pdf', bbox_inches='tight')

    # all fixed, generate ILP:
    constraints = []

    # consistency and matching 1-to-1
    constraints.append("\ Matching and consistency constraints")
    # sorting just to make it nicer looking:
    for (gene, copy_a) in sorted(edges):
        copy_set_b = edges[(gene, copy_a)]
        if len(copy_set_b) == 1:
            constraints.append("%s = 1" % matching_edge_name(gene, copy_a, copy_b, Ext.HEAD))
        else:
            for copy_b in copy_set_b:
                constraints.append("%s - %s = 0" % (
                    matching_edge_name(gene, copy_a, copy_b, Ext.TAIL),
                    matching_edge_name(gene, copy_a, copy_b, Ext.HEAD)))
            constraints.append(
                " + ".join([matching_edge_name(gene, copy_a, copy_b, Ext.TAIL) for copy_b in copy_set_b]) + " = 1")

    constraints.append("\ Balancing:")
    balancing_genes_A = {g: range(gene_count["A"][g] + 1, gene_count["B"][g] + 1) for g in total_gene_count.iterkeys()
                         if gene_count["A"][g] < gene_count["B"][g]}
    balancing_genes_B = {g: range(gene_count["B"][g] + 1, gene_count["A"][g] + 1) for g in total_gene_count.iterkeys()
                         if gene_count["B"][g] < gene_count["A"][g]}

    for genome, balancing in [("A", balancing_genes_A), ("B", balancing_genes_B)]:
        constraints.append("\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing):
            # check if fixed:
            if (gene_i, copy_i, ext_i) in balancing_fix[genome]:
                gene_j, copy_j, ext_j = balancing_fix[genome][(gene_i, copy_i, ext_i)]
                if (gene_i, copy_i, ext_i) < (gene_j, copy_j, ext_j):
                    constraints.append(
                        "%s = 1" % balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j))
            # if not, matching 1-to-1:
            else:
                constraints.append(
                    " + ".join([balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j) for
                                gene_j, copy_j, ext_j in
                                balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) if
                                (gene_i, copy_i, ext_i) != (gene_j, copy_j, ext_j)]) + " = 1")
    constraints.append("\ Labelling")

    #
    # for each adjacency, fix label:
    constraints.append("\\ Adjacency have the same label:")
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        for i, j in adjacency_list(genome, total_gene_count):
            v_i = vertex_name(genome_name, *i)
            v_j = vertex_name(genome_name, *j)
            # if already fixed, skip
            if y_label[v_i] in y_fix and y_label[v_j] in y_fix:
                continue
            # if the edge is 0 for sure, also skip:
            constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (y_label[v_i], y_label[v_j], v_i, v_j))
    #
    constraints.append("\\ Matching edges with the same label:")
    for (gene, copy_a) in sorted(edges):
        copy_set_b = edges[(gene, copy_a)]
        for ext in [Ext.HEAD, Ext.TAIL]:
            y_i = y_label[vertex_name("A", gene, copy_a, ext)]
            # if edge is set, just make the y_i's equal;
            if len(copy_set_b) == 1:
                y_j = y_label[vertex_name("B", gene, list(copy_set_b)[0], ext)]
                # skip if this y_i's are already fixed
                if y_i in y_fix and y_j in y_fix:
                    continue
                constraints.append("y_%s - y_%s = 0 " % (y_i, y_j))
            else:
                # if edge not set, add both ineqs.
                for copy_b in copy_set_b:
                    y_j = y_label[vertex_name("B", gene, copy_b, ext)]
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (
                        y_i, y_j, y_i, matching_edge_name(gene, copy_a, copy_b, ext), y_i))
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (
                        y_j, y_i, y_j, matching_edge_name(gene, copy_a, copy_b, ext), y_j))

    constraints.append("\\ Balancing edges with same label:")
    for genome, balancing in [("A", balancing_genes_A), ("B", balancing_genes_B)]:
        constraints.append("\\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()):
            for gene_j, copy_j, ext_j in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()):
                if (gene_i, copy_i, ext_i) >= (gene_j, copy_j, ext_j):
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
            # if i in z_fix and z_fix[i] == 0:
            #     continue
            # if i in z_fix and z_fix[i] == 1:
            #     constraints.append("z_%s = 1" % i)
            if i not in z_fix:
                constraints.append("%d z_%s - y_%s <= 0" % (i, i, i))
    #
    # # number of genes, to fix distance:
    constraints.append("n = %d" % (sum(total_gene_count.itervalues())))
    # # number of fixed cycles
    constraints.append("c = %d" % (sum(z_fix.itervalues())))
    # for g in sorted(total_gene_count):
    #     print g,total_gene_count[g]

    #
    # # bounds:
    bounds = []
    for i in sorted(y_label.itervalues()):
        if i not in y_fix:
            bounds.append("y_%d <= %d" % (i, i))
    #
    # # variables:
    binary = []
    #
    # # matching edges
    # matching edges, skipping fixed pairs.
    matching = ["\ match"]
    # for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
    for (gene, copy_a), copy_set_b in sorted(edges.items(), key=operator.itemgetter(0)):
        # fixed vars, just the head, to know the fixed value when parsing;
        if len(copy_set_b) == 1:
            matching.append(matching_edge_name(gene, copy_a, copy_b, Ext.HEAD))
        # non fixed, both head and tail;
        else:
            for copy_b in copy_set_b:
                for ext in [Ext.HEAD, Ext.TAIL]:
                    matching.append(matching_edge_name(gene, copy_a, copy_b, ext))

    print "%d matching edges" % len(matching)
    # print "Potentially %d matching edges" % sum([2*x ** 2 for x in gene_count.itervalues()])
    binary.extend(matching)
    #
    # balancing edges:
    balancing_edges = [balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j) for genome, balancing
                       in [("A", balancing_genes_A), ("B", balancing_genes_B)] for gene_i, copy_i, ext_i in
                       balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) for gene_j, copy_j, ext_j
                       in balancing_extremities(balancing, exclude=balancing_fix[genome].keys()) if
                       (gene_i, copy_i, ext_i) < (gene_j, copy_j, ext_j)]
    print "%d balancing edges" % len(balancing_edges)
    binary.extend(balancing_edges)
    #
    # z cycles:
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if i in z_fix:  # and z_fix[i] == 0:
            continue
        if vertex[0] == "B":
            continue
        binary.append("z_%d" % i)
    #
    # # Y label are general:
    # TODO: remove unused y' and z's from model. If y=1, it can be removed, just set z=1.
    general = []
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if i not in y_fix:
            general.append("y_%d" % i)
    #
    # # number of genes and fixed cycles:
    general.append("n")
    general.append("c")
    # # objective function:
    objective = ["obj: n - c - " + " - ".join(
        ["z_%d" % i for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)) if
         vertex[0] == "A" and i not in z_fix])]

    # write:
    with open(output, "w") as f:
        for header, lines in [("Minimize", objective), ("Subject to", constraints),
                              ("Bounds", bounds), ("Binary", binary), ("General", general)]:
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

    # not verbose:
    # model.setParam('OutputFlag', False)
    # MIP focus, from 0 to 3:
    model.params.MIPFocus = 1  # best solutions, less focus on bounds.
    model.optimize()

    if model.status != GRB.Status.INFEASIBLE:
        print('FINISHED: Best objective: %g' % model.objVal)
        print('Optimization ended with status %d' % model.status)
        model.write(filename + '.sol')

    if model.status == GRB.INFEASIBLE:
        model.computeIIS()
        model.write("unfeasible.lp")
        print('\nThe following constraint(s) cannot be satisfied:')
        for c in model.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
    else:
        z = 0
        n = 0
        for v in model.getVars():
            if v.varName == "n":
                n = v.x
            if v.varName.startswith("z") and v.x == 1:
                z += 1
        print "N: %d  cycles:%d" % (n, z)

    return model


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generates and optionally solve an ILP for the DCJ duplication and indel distance.")
    parser.add_argument("-s", "--solve", action="store_true", default=False, help="Solve the model with Gurobi.")
    parser.add_argument("file", type=str, help="Genomes file.")
    parser.add_argument("g1", type=int, help="Index of genome 1 in the file. (0-indexed).")
    parser.add_argument("g2", type=int, help="Index of genome 1 in the file. (0-indexed).")
    param = parser.parse_args()

    genomes = file_ops.open_genome_file(param.file, as_list=True)
    g1 = genomes[param.g1]
    g2 = genomes[param.g2]

    filename = "%s_%s_%s.lp" % (os.path.basename(param.file), g1.name,  g2.name)
    dcj_dupindel_ilp(g1, g2, filename)

    if param.solve:
        model = solve_ilp(filename)
