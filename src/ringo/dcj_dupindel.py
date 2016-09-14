#!/usr/bin/env python2
import collections

import pyximport;
import re

pyximport.install()
from model import Ext, Chromosome
import argparse
import copy
import operator
import networkx as nx
from networkx.algorithms import connected_components, dfs_successors
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
        for gene, copy_dict in gene_count[genome].iteritems():
            for copy_i, type_i in copy_dict.iteritems():
                y_label[(genome, gene, copy_i, Ext.HEAD)] = idx
                idx += 1
                y_label[(genome, gene, copy_i, Ext.TAIL)] = idx
                idx += 1
    return y_label


def balancing_extremities(gene_copies, exclude=None):
    if exclude is None:
        exclude = set()
    for gene, copy_dict in gene_copies.iteritems():
        for copy_i, type_i in copy_dict.iteritems():
            if type_i == CopyType.BALANCING:
                if (gene, copy_i, Ext.HEAD) not in exclude:
                    yield (gene, copy_i, Ext.HEAD)
                if (gene, copy_i, Ext.TAIL) not in exclude:
                    yield (gene, copy_i, Ext.TAIL)


class CopyType:
    REAL, BALANCING = ['r', 'b']


# just to better print cycles and paths in order:
def sort_component(G, comp, fmt=True):
    # get degree-1 vertices:
    degree_one = [v for v in comp if G.degree(v) == 1]
    if len(degree_one) > 0:
        v = degree_one[0]
    else:
        v = list(comp)[0]
    # import ipdb; ipdb.set_trace()
    succ = dfs_successors(G, v)
    initial = v
    sort = []
    while True:
        if fmt:
            sort.append("%s%s_%s%s" % v)
        else:
            sort.append(v)
        if v not in succ:
            break
        v = succ[v][0]
        if v == initial:
            break
    return sort


def add_capping_genes(genome_a, genome_b):
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    for genome in [genome_a, genome_b]:
        copy_idx = 1
        for c in genome.chromosomes:
            if not c.circular:
                c.gene_order.append(0)
                c.copy_number.append(copy_idx)
                copy_idx += 1
                c.circular = True
        for i in range(genome.n_chromosomes(), max_chromosomes):
            genome.add_chromosome(Chromosome([0], copy_number=[copy_idx], circular=True))
            copy_idx += 1


def fix_cycle_y_z(comp, y_label, y_fix, z_fix, vertices_to_remove):
    # get indexes of the y_i:
    indexes = [(v, y_label[v]) for v in comp]
    min_label = min([x[1] for x in indexes])
    for v, label in indexes:
        y_fix[label] = min_label
        z_fix[label] = 0
    z_fix[min_label] = 1
    vertices_to_remove.extend(comp)


def build_gene_copies_dict(all_genes, genome_a, genome_b):
    # store the copy number for each 'real' gene in each genome:
    gene_copies_a = genome_a.gene_copies()
    gene_copies_b = genome_b.gene_copies()
    gene_copies = {"A": {}, "B": {}}
    for gene in all_genes:
        gene_copies["A"][gene] = {cn: CopyType.REAL for cn in gene_copies_a[gene]}
        gene_copies["B"][gene] = {cn: CopyType.REAL for cn in gene_copies_b[gene]}

    # now complete the list with 'balancing' genes
    for gene in all_genes:
        copy_a = len(gene_copies["A"][gene])
        copy_b = len(gene_copies["B"][gene])
        # get the max copy number, or 0 if set is empty, to know which label is free for the newly added
        # balancing copies;
        max_copy_a = max([cn for cn, c_type in gene_copies["A"][gene].iteritems()] + [0])
        max_copy_b = max([cn for cn, c_type in gene_copies["B"][gene].iteritems()] + [0])

        if copy_a < copy_b:
            gene_copies["A"][gene].update({bal_copy: CopyType.BALANCING for bal_copy in
                                           xrange(max_copy_a + 1, max_copy_a + 1 + copy_b - copy_a)})
        if copy_b < copy_a:
            gene_copies["B"][gene].update({bal_copy: CopyType.BALANCING for bal_copy in
                                           xrange(max_copy_b + 1, max_copy_b + 1 + copy_a - copy_b)})

    return gene_copies


def plot_bp(filename, master_graph, gene_copies, fixed_matching):
    # add isolated vertices (balancing extremities that are not fixed already)
    for genome_i in ["A", "B"]:
        for gene_i, copy_i, ext_i in balancing_extremities(gene_copies[genome_i]):
            if (gene_i, copy_i) not in fixed_matching:
                master_graph.add_node((genome_i, gene_i, copy_i, ext_i))
    # Relabel nodes to make it easier to read:
    mapping = {}
    normal = []
    balancing = []
    be = {genome_i: list(balancing_extremities(gene_copies[genome_i])) for genome_i in ["A", "B"]}
    for v in master_graph.nodes():
        genome_i, gene_i, copy_i, ext_i = v
        mapping[v] = "$%s%s_{(%s)}^%s$" % v
        if (gene_i, copy_i, ext_i) in be[genome_i]:
            balancing.append(mapping[v])
        else:
            normal.append(mapping[v])
    master_graph = nx.relabel_nodes(master_graph, mapping)

    # Graphviz position:
    # pos = nx.nx_agraph.graphviz_layout(master_graph, prog="fdp")

    # custom position:
    x_pos = 0
    y_pos = {"A": 1, "B": 0}
    pos = {}
    for comp in sorted(connected_components(master_graph), key=lambda c: (-len(c), min(c))):
        last_v = None
        for v in sort_component(master_graph, comp, fmt=False):
            if last_v == v[1]:
                x_pos += 1
            last_v = v[1]
            pos[v] = (x_pos, y_pos[v[1]])
        x_pos += 1
        if x_pos > 7:
            x_pos = 0
            y_pos["A"] -= 2
            y_pos["B"] -= 2

    # draw and save:
    for nodelist, color in [(normal, "lightgray"), (balancing, "lightblue")]:
        nx.draw(master_graph, pos, font_size=5, nodelist=nodelist, node_color=color, linewidths=0.1, width=0.5,
                node_size=400,
                with_labels=True)
    plt.savefig(filename, bbox_inches='tight')


######################################################################
# MAIN function
######################################################################

def dcj_dupindel_ilp(genome_a, genome_b, output, skip_balancing=False, fix_vars=True, solve=False):
    def solve_ilp(timelimit=60):
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

            z = n = c = 0
            solution_matching = collections.defaultdict(list)
            matching_regexp = re.compile("x_A(\d+)_(\d+)h,B(\d+)_(\d+)h")
            # get basic vars and matching:
            for v in model.getVars():
                if v.varName == "n":
                    n = v.x
                elif v.varName == "c":
                    c = v.x
                elif v.varName.startswith("z") and v.x >= 0.9:
                    z += 1
                else:
                    m = matching_regexp.match(v.varName)
                    if m is not None and v.x == 1:
                        g_a, c_a, g_b, c_b = map(int, m.groups())
                        solution_matching[g_a].append((c_a, c_b))

            from parse_orthology import build_correct_matching, parse_orthology_quality
            correct_matching = build_correct_matching(genome_a, genome_b)
            tp, fp, fn = parse_orthology_quality(solution_matching, correct_matching)

            print "N: %d  cycles:%d (%d fixed, %d from opt)" % (n, z + c, c, z)
            print "Orthology. TP:%d  FP:%d  FN:%d" % (len(tp), len(fp), len(fn))
            # print match_edges
            # Now, analyse the BP graph, for the incomplete matching model, to find AA-, BB- and AB- components:
            master_graph = nx.Graph()
            # fixed vars:
            # add matching edges of genes with single copy:
            # for (gene, copy_a), copy_j in match_edges.iteritems():
            for gene, pair_list in solution_matching.iteritems():
                for copy_a, copy_b in pair_list:
                    for ext in [Ext.HEAD, Ext.TAIL]:
                        master_graph.add_edge(("A", gene, copy_a, ext), ("B", gene, copy_b, ext))

            # add adjacency edges:
            for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
                for (g_i, copy_a, e_i), (g_j, copy_j, e_j) in genome.adjacency_iter_with_copies():
                    master_graph.add_edge((genome_name, g_i, copy_a, e_i), (genome_name, g_j, copy_j, e_j))

            count = {"A": 0, "B": 0, "AB": 0}
            c = 0
            # print "C:", len([x for x in connected_components(master_graph)])
            for comp in connected_components(master_graph):
                degree_one = [v for v in comp if master_graph.degree(v) == 1]
                if len(degree_one) == 0:
                    c += 1
                else:
                    if len(degree_one) != 2:
                        import ipdb;
                        ipdb.set_trace()
                    if degree_one[0][0] == degree_one[1][0]:
                        count[degree_one[0][0]] += 1
                    else:
                        count["AB"] += 1
            print count
            if skip_balancing:
                print "Corrected distance: %d" % (model.objVal + count["AB"] / 2)

        return model

    # copy genomes to possibly make some changes:
    genome_a = copy.deepcopy(genome_a)
    genome_b = copy.deepcopy(genome_b)

    add_capping_genes(genome_a, genome_b)

    # since the gene set might be different for each genome, find all genes:
    all_genes = genome_a.gene_set().union(genome_b.gene_set())
    # find all gene copies
    gene_copies = build_gene_copies_dict(all_genes, genome_a, genome_b)
    # count balancing genes:
    bal = {
        g: sum([len([c for c in gene_copies[g][gene].itervalues() if c == CopyType.BALANCING]) for gene in all_genes])
        for g in ["A", "B"]}

    print "Balancing genes:A=%(A)d, B=%(B)d" % bal
    # define the y labels (vertex = genome,gene,copy,ext) -> integer 1..n
    y_label = define_y_label(gene_copies)

    # store all possible matchings (edges) from each family:
    fixed_matching = {}
    possible_matching = {}
    for gene in all_genes:
        # if only 1 copy, matching is fixed:
        if len(gene_copies["A"][gene]) == 1:
            # fix the matching, then remove from the available copies
            copy_a, type_a = gene_copies["A"][gene].items()[0]
            copy_j, type_b = gene_copies["B"][gene].items()[0]
            fixed_matching[(gene, copy_a)] = copy_j
        else:
            possible_matching[gene] = {"A": {copy_i for copy_i, type_i in gene_copies["A"][gene].items()},
                                       "B": {copy_i for copy_i, type_i in gene_copies["B"][gene].items()}}

    # Build the BP graph of fixed matchings to try to find more variables to fix:
    y_fix = {}
    z_fix = {}
    balancing_fix = {"A": {}, "B": {}}

    if fix_vars:
        master_graph = nx.Graph()
        # fixed vars:

        # add matching edges of genes with single copy:
        for (gene, copy_a), copy_j in fixed_matching.iteritems():
            for ext in [Ext.HEAD, Ext.TAIL]:
                master_graph.add_edge(("A", gene, copy_a, ext), ("B", gene, copy_j, ext))

        # add adjacency edges:
        for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
            for (g_i, copy_a, e_i), (g_j, copy_j, e_j) in genome.adjacency_iter_with_copies():
                master_graph.add_edge((genome_name, g_i, copy_a, e_i), (genome_name, g_j, copy_j, e_j))

        # Search components to fix:
        rescan = True
        edges_to_add = []
        vertices_to_remove = []
        ab_components = set()
        while rescan:
            rescan = False
            # Pre-scan:
            # add and remove edges detected from previous rounds:
            master_graph.add_edges_from(edges_to_add)
            master_graph.remove_nodes_from(vertices_to_remove)
            edges_to_add = []
            vertices_to_remove = []
            # fix AB-components; while I have at least 2, join pairs arbitrarily:
            while len(ab_components) > 1:
                a_i, b_i = ab_components.pop()
                a_j, b_j = ab_components.pop()
                master_graph.add_edge(a_i, a_j)
                balancing_fix["A"][a_i[1:]] = a_j[1:]
                balancing_fix["A"][a_j[1:]] = a_i[1:]
                master_graph.add_edge(b_i, b_j)
                balancing_fix["B"][b_i[1:]] = b_j[1:]
                balancing_fix["B"][b_j[1:]] = b_i[1:]

            # Now I search for vertices that have only balancing vertices as matching
            # candidates. If that is the case, I can fix them arbitrarly.
            fix_only_bal = True
            if fix_only_bal:
                for gene in sorted(possible_matching):
                    set_a = possible_matching[gene]["A"]
                    set_b = possible_matching[gene]["B"]
                    if all([gene_copies["A"][gene][copy_a] == CopyType.BALANCING for copy_a in set_a]) or all(
                            [gene_copies["B"][gene][copy_b] == CopyType.BALANCING for copy_b in set_b]):
                        for copy_a, copy_b in zip(set_a, set_b):
                            fixed_matching[(gene, copy_a)] = copy_b
                            # save edges to add to graph:
                            for ext in [Ext.HEAD, Ext.TAIL]:
                                # edges_to_add.append((("A", gene, copy_a, ext), ("B", gene, copy_b, ext)))
                                master_graph.add_edge(("A", gene, copy_a, ext), ("B", gene, copy_b, ext))
                        rescan = True
                        # remove from possible matching:
                        del possible_matching[gene]

            # now loop for each connected component, fixing cycles and trying to close paths to cycles when possible.
            for comp in connected_components(master_graph):
                # can only consider even components;
                if len(comp) % 2 != 0:
                    continue
                # get degree-1 vertices:
                degree_one = [v for v in comp if master_graph.degree(v) == 1]
                # if two degree one vertices, it is a path;
                if len(degree_one) == 2:
                    genome_i, g_i, copy_a, e_i = degree_one[0]
                    genome_j, g_j, copy_j, e_j = degree_one[1]

                    # 1 - check if both nodes are balancing, to find AA-, BB- and AB- paths that can be fixed.
                    i_is_balancing = g_i != 0 and gene_copies[genome_i][g_i][copy_a] == CopyType.BALANCING
                    j_is_balancing = g_j != 0 and gene_copies[genome_j][g_j][copy_j] == CopyType.BALANCING

                    if i_is_balancing and j_is_balancing:
                        # open-path, both ends are balancing.
                        # If AA- or BB-path, close it to a cycle:
                        if genome_i == genome_j:
                            # fix the cycle:
                            fix_cycle_y_z(comp, y_label, y_fix, z_fix, vertices_to_remove)
                            # fix the balancing variables if we have them:
                            if not skip_balancing:
                                balancing_fix[genome_i][degree_one[0][1:]] = degree_one[1][1:]
                                balancing_fix[genome_i][degree_one[1][1:]] = degree_one[0][1:]
                        else:
                            # If not, it is AB-, add to the list to try to make pairs.
                            if skip_balancing:  # if not using balancing edges, I can fix the AB directly, instead of
                                # doing the merge in pairs;
                                fix_cycle_y_z(comp, y_label, y_fix, z_fix, vertices_to_remove)
                            else:
                                # merge in pairs:
                                ab_components.add(tuple(sorted(degree_one)))
                                if len(ab_components) > 1:
                                    rescan = True

                    # Not open path; then, check if the path has homologous extremities at both ends, so I can close
                    # to a path:
                    elif genome_i != genome_j and g_i == g_j and e_i == e_j:
                        # invert to put genome A always in variables _i :
                        if genome_j == "A":
                            genome_i, g_i, copy_a, e_i, genome_j, g_j, copy_j, e_j = genome_j, g_j, copy_j, e_j, genome_i, g_i, copy_a, e_i

                        # check conflict, only add edge if it's in the allowed edges:
                        if g_i in possible_matching and copy_a in possible_matching[g_i]["A"] and copy_j in \
                                possible_matching[g_i]["B"]:
                            fixed_matching[(g_i, copy_a)] = copy_j
                            # save edges to add to graph:
                            for ext in [Ext.HEAD, Ext.TAIL]:
                                edges_to_add.append((("A", g_i, copy_a, ext), ("B", g_i, copy_j, ext)))
                            # new edges, re-scan:
                            rescan = True

                            # remove possible edges from other copies:
                            possible_matching[g_i]["A"].remove(copy_a)
                            possible_matching[g_i]["B"].remove(copy_j)
                            # if now there is just one possibility, also fix:
                            if len(possible_matching[g_i]["A"]) == 1:
                                copy_a = possible_matching[g_i]["A"].pop()
                                copy_b = possible_matching[g_i]["B"].pop()
                                fixed_matching[(g_i, copy_a)] = copy_b
                                del possible_matching[g_i]
                                # save edges to add to graph:
                                for ext in [Ext.HEAD, Ext.TAIL]:
                                    edges_to_add.append((("A", g_i, copy_a, ext), ("B", g_i, copy_b, ext)))

                # if there are no degree one vertices, it is a cycle; I can fix the y_i and z_i for this cycle:
                elif len(degree_one) == 0:
                    fix_cycle_y_z(comp, y_label, y_fix, z_fix, vertices_to_remove)
                    rescan = True

    # DRAW:
    draw_bp = False
    if draw_bp:
        plot_bp('graph.pdf')

    # all fixed, generate ILP

    # to make it easier to find the matching edges, specially when limiting edges from balancing genes,
    # I will build a gene connections graph;
    gene_connection = nx.DiGraph()  # make it directed, so the vertex of A is always 1st on the edge tuple.
    for gene in possible_matching.iterkeys():
        set_a = possible_matching[gene]["A"]
        set_b = possible_matching[gene]["B"]
        # All vs all model:
        for copy_a in set_a:
            for copy_b in set_b:
                gene_connection.add_edge(("A", gene, copy_a), ("B", gene, copy_b))

    # Start building constraints:
    constraints = []

    # consistency and matching 1-to-1

    # Fixed matching:
    # sorting just to make it nicer looking:
    constraints.append("\ Fixed matching:")
    for (gene, copy_a), copy_b in sorted(fixed_matching.items(), key=lambda pair: pair[0]):
        constraints.append("%s = 1" % matching_edge_name(gene, copy_a, copy_b, Ext.TAIL))
        constraints.append("%s = 1" % matching_edge_name(gene, copy_a, copy_b, Ext.HEAD))

    # HEAD TAIL consistency:
    constraints.append("\ Head/Tail consistency:")
    for (_, gene_a, copy_a), (_, gene_b, copy_b) in gene_connection.edges_iter():
        constraints.append("%s - %s = 0" % (
            matching_edge_name(gene_a, copy_a, copy_b, Ext.TAIL),
            matching_edge_name(gene_a, copy_a, copy_b, Ext.HEAD)))

    # 1 Matching per node :
    constraints.append("\ Degree 1 per node (Matching):")
    # for all vertices:
    for v in gene_connection.nodes_iter():
        # find the incident edges:
        if v[0] == "A":
            edges = gene_connection.out_edges_iter
        else:
            edges = gene_connection.in_edges_iter
        incident = [matching_edge_name(gene_a, copy_a, copy_b, Ext.TAIL) for
                    (_, gene_a, copy_a), (_, gene_b, copy_b) in edges(v)]
        # sum of incidents is 1:
        constraints.append("%s = 1" % (" + ".join(incident)))

    if not skip_balancing:
        constraints.append("\ Balancing:")

        for genome in ["A", "B"]:
            constraints.append("\ Genome %s" % genome)
            for gene_i, copy_a, ext_i in balancing_extremities(gene_copies[genome]):
                # check if fixed:
                if (gene_i, copy_a, ext_i) in balancing_fix[genome]:
                    gene_j, copy_j, ext_j = balancing_fix[genome][(gene_i, copy_a, ext_i)]
                    if (gene_i, copy_a, ext_i) < (gene_j, copy_j, ext_j):
                        constraints.append(
                            "%s = 1" % balancing_edge_name(genome, gene_i, copy_a, ext_i, gene_j, copy_j, ext_j))
                # if not, matching 1-to-1:
                else:
                    constraints.append(
                        " + ".join([balancing_edge_name(genome, gene_i, copy_a, ext_i, gene_j, copy_j, ext_j) for
                                    gene_j, copy_j, ext_j in
                                    balancing_extremities(gene_copies[genome], exclude=balancing_fix[genome].keys())
                                    if
                                    (gene_i, copy_a, ext_i) != (gene_j, copy_j, ext_j)]) + " = 1")


    constraints.append("\ Labelling")
    # for each adjacency, fix the label of adjacent genes:

    constraints.append("\\ Adjacent nodes have the same label:")
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        for (g_i, copy_a, ext_i), (g_j, copy_j, ext_j) in genome.adjacency_iter_with_copies():
            v_i = genome_name, g_i, copy_a, ext_i
            v_j = genome_name, g_j, copy_j, ext_j
            # if already fixed, skip
            if y_label[v_i] in y_fix and y_label[v_j] in y_fix:
                continue
            # if the edge is 0 for sure, also skip:
            constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (y_label[v_i], y_label[v_j], v_i, v_j))
    #
    constraints.append("\\ Matching extremities have the same label:")

    # if extremities are matched, but I don't know the y_i (cycle was not closed in the fixing phase),
    # then I know that the y_i's of these extremities are equal:
    constraints.append("\\ Fixed matching without fixed y_i:")
    for (gene, copy_a) in sorted(fixed_matching):
        copy_j = fixed_matching[(gene, copy_a)]
        for ext in [Ext.HEAD, Ext.TAIL]:
            y_i = y_label[("A", gene, copy_a, ext)]
            y_j = y_label[("B", gene, copy_j, ext)]
            # only add if this y_i's aren't already fixed
            if y_i not in y_fix and y_j not in y_fix:
                constraints.append("y_%s - y_%s = 0 " % (y_i, y_j))

    # for the "open" matching, for each edge I add the "y fixing" restrictions, that force the y_i's
    # to be equal whenever the edge variable is set to 1.
    constraints.append("\\ Open matching:")
    for (_, gene_a, copy_a), (_, gene_b, copy_b) in gene_connection.edges_iter():
        for ext in [Ext.HEAD, Ext.TAIL]:
            y_a = y_label[("A", gene_a, copy_a, ext)]
            y_b = y_label[("B", gene_b, copy_b, ext)]
            constraints.append(
                "y_%s - y_%s + %s %s <= %d" % (
                    y_a, y_b, y_a, matching_edge_name(gene_a, copy_a, copy_b, ext), y_a))
            constraints.append(
                "y_%s - y_%s + %s %s <= %d" % (
                    y_b, y_a, y_b, matching_edge_name(gene_a, copy_a, copy_b, ext), y_b))

    if not skip_balancing:
        constraints.append("\\ Balancing edges have same label:")
        for genome in ["A", "B"]:
            constraints.append("\\ Genome %s" % genome)
            for gene_i, copy_a, ext_i in balancing_extremities(gene_copies[genome],
                                                               exclude=balancing_fix[genome].keys()):
                for gene_j, copy_j, ext_j in balancing_extremities(gene_copies[genome],
                                                                   exclude=balancing_fix[genome].keys()):
                    if (gene_i, copy_a, ext_i) >= (gene_j, copy_j, ext_j):
                        continue
                    y_i = y_label[(genome, gene_i, copy_a, ext_i)]
                    y_j = y_label[(genome, gene_j, copy_j, ext_j)]
                    # should not have someone here if I'm excluding fixed edges:
                    if y_i in y_fix and y_j in y_fix:
                        continue
                    constraints.append("y_%s - y_%s + %s %s <= %d" % (
                        y_i, y_j, y_i, balancing_edge_name(genome, gene_i, copy_a, ext_i, gene_j, copy_j, ext_j), y_i))
                    constraints.append("y_%s - y_%s + %s %s <= %d" % (
                        y_j, y_i, y_j, balancing_edge_name(genome, gene_i, copy_a, ext_i, gene_j, copy_j, ext_j), y_j))

    # z variables: since all cycles have to contains vertices from both genomes, we only add z variables
    # for genome A, that have smallest labels, so a genome B z variable will never be =1.
    constraints.append("\\ Z variables")
    for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "A":
            if i not in z_fix:
                constraints.append("%d z_%s - y_%s <= 0" % (i, i, i))
    #
    # # number of genes, to fix distance:
    n_genes = sum([len(copies) for copies in gene_copies["A"].itervalues()])
    constraints.append("n = %d" % n_genes)
    # # number of fixed cycles
    constraints.append("c = %d" % (sum(z_fix.itervalues())))

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
    matching = ["\ Fixed matching:"]
    for (gene, copy_a), copy_b in fixed_matching.iteritems():
        matching.append(matching_edge_name(gene, copy_a, copy_b, Ext.TAIL))
        matching.append(matching_edge_name(gene, copy_a, copy_b, Ext.HEAD))

    matching.append("\ Open matching:")
    for (_, gene_a, copy_a), (_, gene_b, copy_b) in gene_connection.edges_iter():
        for ext in [Ext.HEAD, Ext.TAIL]:
            matching.append(matching_edge_name(gene_a, copy_a, copy_b, ext))

    print "%d fixed matching edges" % (len(fixed_matching) * 2)
    print "%d open matching edges" % (len(gene_connection.edges()) * 2)
    binary.extend(matching)
    if not skip_balancing:
        # balancing edges:
        balancing_edges = [balancing_edge_name(genome, gene_i, copy_a, ext_i, gene_j, copy_j, ext_j) for
                           genome
                           in ["A", "B"] for gene_i, copy_a, ext_i in
                           balancing_extremities(gene_copies[genome], exclude=balancing_fix[genome].keys()) for
                           gene_j, copy_j, ext_j
                           in balancing_extremities(gene_copies[genome], exclude=balancing_fix[genome].keys()) if
                           (gene_i, copy_a, ext_i) < (gene_j, copy_j, ext_j)]
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
    z_obj = " - ".join(["z_%d" % i for vertex, i in sorted(y_label.items(), key=operator.itemgetter(1)) if
                        vertex[0] == "A" and i not in z_fix])

    objective = ["obj: n - c %s" % ("- " + z_obj if len(z_obj) > 0 else "")]
    # write ILP:
    with open(output, "w") as f:
        for header, lines in [("Minimize", objective), ("Subject to", constraints),
                              ("Bounds", bounds), ("Binary", binary), ("General", general)]:
            print >> f, header
            print >> f, "\n".join(lines)

    if solve:
        model = solve_ilp(timelimit=60)
        return model


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generates and optionally solve an ILP for the DCJ duplication and indel distance.")
    parser.add_argument("-s", "--solve", action="store_true", default=False, help="Solve the model with Gurobi.")
    parser.add_argument("-sb", "--skip_balancing", action="store_true", default=False,
                        help="Do not add balancing edges.")
    parser.add_argument("-sf", "--skip_fixing", action="store_false", default=True,
                        help="Do not try to fix variables.")
    parser.add_argument("-t", "--timelimit", type=int, default=60,
                        help="Time limit in seconds for the solver. (default 60 secs.)")
    input_type = parser.add_mutually_exclusive_group(required=True)
    input_type.add_argument("-g", type=str, nargs=3, help="Genomes file, idx 1 and 2 of genomes (0-indexed).")
    input_type.add_argument("-c", type=str, nargs=2, help="Two coser files.")

    param = parser.parse_args()
    if param.g is not None:
        filename, n1, n2 = param.g
        genomes = file_ops.open_genomes_with_copy_number(filename)
        g1 = genomes[int(n1)]
        g2 = genomes[int(n2)]
    elif param.c is not None:
        g1 = file_ops.open_coser_genome(param.c[0])
        g2 = file_ops.open_coser_genome(param.c[1])
        filename = "ilp"
    filename = "%s_%s_%s%s.lp" % (filename, g1.name, g2.name, "_nobal" if param.skip_balancing else "")
    dcj_dupindel_ilp(g1, g2, filename, skip_balancing=param.skip_balancing, fix_vars=param.skip_fixing,
                     solve=param.solve)
