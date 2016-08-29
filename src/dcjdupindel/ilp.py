#!/usr/bin/env python2
import argparse
import pyximport; pyximport.install()
import os
import sys
import itertools
import operator
from model import Genome, Chromosome

# TODO: this is a hack to import from other directory; should use packages
ringo_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../ringo"))
print ringo_path
os.sys.path.insert(0, ringo_path)
import file_ops

class Ext:
    HEAD = 'h'
    TAIL = 't'


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


def balancing_extremities(balancing):
    for gene, copies in balancing.iteritems():
        for copy in copies:
            yield gene, copy, Ext.HEAD
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
    adj_list = [(0, copy_number[0], Ext.TAIL)]
    copy_number[0] += 1
    for gene in chrom.gene_order:
        if gene > 0:
            orientation = [Ext.TAIL, Ext.HEAD]
        else:
            orientation = [Ext.HEAD, Ext.TAIL]
        adj_list.extend([(abs(gene), copy_number[abs(gene)], ext) for ext in orientation])
        copy_number[abs(gene)] += 1
    adj_list.append((0, copy_number[0], Ext.TAIL))
    copy_number[0] += 1
    return adj_list


####################################################################
# ILP
####################################################################
def dcj_dupindel_ilp(genome_a, genome_b, output):
    # count of each gene on each genome
    gene_count_a = genome_a.gene_count()
    gene_count_b = genome_b.gene_count()

    # find all genes
    gene_count = {g: max(gene_count_a[g], gene_count_b[g]) for g in
                  set(gene_count_a.keys()).union(set(gene_count_b.keys()))}

    # build constraints:
    constraints = []
    # consistency:
    constraints.append("\ Consistency constraints")
    for gene, copies in gene_count.iteritems():
        constraints.append("\\ gene %s" % gene)
        for i, j in itertools.product(xrange(1, copies + 1), xrange(1, copies + 1)):
            constraints.append(
                "%s - %s = 0" % (matching_edge_name(gene, i, j, Ext.TAIL), matching_edge_name(gene, i, j, Ext.HEAD)))

    # unique matching: (1-to-1) (only tail, it should suffice, together with the consistency)
    constraints.append("\\ Unique matching:")
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, copies + 1):
            constraints.append(
                " + ".join([matching_edge_name(gene, i, j, Ext.TAIL) for j in xrange(1, copies + 1)]) + " = 1")
            # adding the HEAD usually is worse, but sometimes makes things a bit faster, which is strange
            # because these are redundant constraints.
            # constraints.append(
            #     " + ".join([matching_edge_name(gene, i, j, Ext.HEAD) for j in xrange(1, copies + 1)]) + " = 1")

    constraints.append("\\ Telomere matching:")  # special case, "0" gene.
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    n_telomeres = 2 * max_chromosomes
    for i in range(1, n_telomeres + 1):
        constraints.append(
            " + ".join([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)]) + " = 1")

    constraints.append("\\ Balancing:")
    balancing_A = {g: range(gene_count_a[g] + 1, gene_count_b[g] + 1) for g in gene_count.iterkeys() if
                   gene_count_a[g] < gene_count_b[g]}
    balancing_B = {g: range(gene_count_b[g] + 1, gene_count_a[g] + 1) for g in gene_count.iterkeys() if
                   gene_count_b[g] < gene_count_a[g]}
    for genome, balancing in [("A", balancing_A), ("B", balancing_B)]:
        constraints.append("\\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing):
            constraints.append(" + ".join([balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j) for
                                           gene_j, copy_j, ext_j in balancing_extremities(balancing) if
                                           (gene_i, copy_i, ext_i) != (gene_j, copy_j, ext_j)]) + " = 1")

    constraints.append("\\ Labelling")
    # define the label -> integer 1..n
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

    # for each adjacency, fix label:
    constraints.append("\\ Adjacency have the same label:")

    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        ext_order_list = build_extremity_order_lists(genome, gene_count, max_chromosomes)
        for ext_order in ext_order_list:
            a = iter(range(len(ext_order)))
            for i, j in itertools.izip(a, a):
                v_i = vertex_name(genome_name, *ext_order[i])
                v_j = vertex_name(genome_name, *ext_order[j])
                constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (y_label[v_i], y_label[v_j], v_i, v_j))

    constraints.append("\\ Matching edges with the same label:")
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, copies + 1):
            for j in xrange(1, copies + 1):
                constraints.append("\\ Edges (%s%s_%s, %s%s_%s)" % ("A", gene, i, "B", gene, j))
                for ext in [Ext.HEAD, Ext.TAIL]:
                    y_i = y_label[vertex_name("A", gene, i, ext)]
                    y_j = y_label[vertex_name("B", gene, j, ext)]
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_i, y_j, y_i, matching_edge_name(gene, i, j, ext), y_i))
                    constraints.append(
                        "y_%s - y_%s + %s %s <= %d" % (y_j, y_i, y_j, matching_edge_name(gene, i, j, ext), y_j))

    constraints.append("\\ Balancing edges with same label:")
    for genome, balancing in [("A", balancing_A), ("B", balancing_B)]:
        constraints.append("\\ Genome %s" % genome)
        for gene_i, copy_i, ext_i in balancing_extremities(balancing):
            for gene_j, copy_j, ext_j in balancing_extremities(balancing):
                if (gene_i, copy_i, ext_i) == (gene_j, copy_j, ext_j):
                    continue
                constraints.append(
                    "\\ Edge (%s%s_%s%s, %s%s_%s%s)" % (genome, gene_i, copy_i, ext_i, genome, gene_j, copy_j, ext_j))
                y_i = y_label[vertex_name(genome, gene_i, copy_i, ext_i)]
                y_j = y_label[vertex_name(genome, gene_j, copy_j, ext_j)]
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
    constraints.append("n = %d" % (sum(gene_count.itervalues()) + n_telomeres / 2))

    # bounds:
    bounds = []
    for i in sorted(y_label.itervalues()):
        bounds.append("y_%d <= %d" % (i, i))

    # EDGE FIXING:
    # 1 - conserved adjacencies (a,b) where one of the families is singleton;
    singleton_genes = {g for g in gene_count.iterkeys() if gene_count_a[g] == gene_count_b[g] == 1}
    adj_dict = {"A": {}, "B": {}}
    for genome, gene_count, name in [(genome_a, gene_count_a, "A"), (genome_b, gene_count_b, "B")]:
        ext_list = build_extremity_order_lists(genome, gene_count, max_chromosomes)
        for ext_order in ext_list:
            a = iter(range(len(ext_order)))
            for i, j in itertools.izip(a, a):
                gene_i, copy_i, ext_i = ext_order[i]
                gene_j, copy_j, ext_j = ext_order[j]
                if gene_count[gene_i] == 1 or gene_count[gene_j] == 1:
                    adj_i = gene_i if ext_i == Ext.HEAD else -gene_i
                    adj_j = gene_j if ext_j == Ext.TAIL else -gene_j
                    if adj_i < 0:
                        adj_i, adj_j = -adj_j, -adj_i
                    adj_dict[name][(adj_i,adj_j)] = (gene_i, copy_i, ext_i), (gene_j, copy_j, ext_j)

    # print adj_dict
    mates = {}
    for adj in sorted(adj_dict["A"].iterkeys()):
        (gene_a_i, gene_a_j) = adj_dict["A"][adj]
        if adj in adj_dict["B"]:
            # import ipdb; ipdb.set_trace()
            gene_b_i, gene_b_j = adj_dict["B"][adj]
            if (gene_a_i[0], gene_a_i[2]) != (gene_b_i[0], gene_b_i[2]):
                gene_b_i, gene_b_j = gene_b_j, gene_b_i
            if (gene_count_a[gene_a_i[0]] == gene_count_b[gene_b_i[0]] == 1) or \
                    (gene_count_a[gene_a_j[0]] == gene_count_b[gene_b_j[0]] == 1):
                # test conflict:
                if gene_a_i not in mates or mates[gene_a_i] == gene_b_i:
                    constraints.append("%s = 1" % matching_edge_name(gene_a_i[0], gene_a_i[1], gene_b_i[1], gene_a_i[2]))
                    # print("%s = 1" % matching_edge_name(gene_a_i[0], gene_a_i[1], gene_b_i[1], gene_a_i[2]))
                    # print gene_a_i, gene_b_i
                    rev_ext = Ext.TAIL if gene_a_i[2] == Ext.HEAD else Ext.HEAD
                    mates[(gene_a_i[0], gene_a_i[1], rev_ext)] = (gene_b_i[0], gene_b_i[1], rev_ext)
                    mates[(gene_b_i[0], gene_b_i[1], rev_ext)] = (gene_a_i[0], gene_a_i[1], rev_ext)
                    # print mates
                    # print

                if gene_a_j not in mates or mates[gene_a_j] == gene_b_j:
                    constraints.append("%s = 1" % matching_edge_name(gene_a_j[0], gene_a_j[1], gene_b_j[1], gene_a_j[2]))
                    # print("%s = 1" % matching_edge_name(gene_a_j[0], gene_a_j[1], gene_b_j[1], gene_a_j[2]))
                    # print gene_a_j, gene_b_j
                    rev_ext = Ext.TAIL if gene_a_j[2] == Ext.HEAD else Ext.HEAD
                    mates[(gene_a_j[0], gene_a_j[1], rev_ext)] = (gene_b_j[0], gene_b_j[1], rev_ext)
                    mates[(gene_b_j[0], gene_b_j[1], rev_ext)] = (gene_a_j[0], gene_a_j[1], rev_ext)
                    # print mates
                    # print


    # variables:
    binary = []

    # matching edges

    # telomere:
    for i in range(1, n_telomeres + 1):
        binary.extend([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)])

    # matching genes:
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, copies + 1):
            for j in xrange(1, copies + 1):
                for ext in [Ext.HEAD, Ext.TAIL]:
                    binary.append(matching_edge_name(gene, i, j, ext))
    # balancing edges:
    for genome, balancing in [("A", balancing_A), ("B", balancing_B)]:
        for gene_i, copy_i, ext_i in balancing_extremities(balancing):
            for gene_j, copy_j, ext_j in balancing_extremities(balancing):
                if (gene_i, copy_i, ext_i) < (gene_j, copy_j, ext_j):
                    binary.append(balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j))

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


def dcj_dupindel_ilp_unitary(genome_a, genome_b, output):
    # count of each gene on each genome
    gene_count_a = genome_a.gene_count()
    gene_count_b = genome_b.gene_count()

    # find all genes
    gene_count = {g: max(gene_count_a[g], gene_count_b[g]) for g in
                  set(gene_count_a.keys()).union(set(gene_count_b.keys()))}

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
    max_chromosomes = max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
    n_telomeres = 2 * max_chromosomes
    for i in range(1, n_telomeres + 1):
        constraints.append(
            " + ".join([matching_edge_name(0, i, j, Ext.TAIL) for j in xrange(1, n_telomeres + 1)]) + " = 1")

    constraints.append("\\ Labelling")
    # define the label -> integer 1..n
    y_label = {}
    idx = 1
    for genome, genome_count in [("A", gene_count_a), ("B", gene_count_b)]:
        # telomeres: (represented by TAIL always)
        for i in range(1, n_telomeres + 1):
            y_label[vertex_name(genome, 0, i, Ext.TAIL)] = idx
            idx += 1
        for gene, copies in gene_count.iteritems():
            for i in xrange(1, genome_count[gene] + 1):
                y_label[vertex_name(genome, gene, i, Ext.HEAD)] = idx
                idx += 1
                y_label[vertex_name(genome, gene, i, Ext.TAIL)] = idx
                idx += 1

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


def solve_ilp(filename):
    # import here, so only if actually solving we will need gurobi.
    from gurobipy import read, GRB
    # pycharm complains of gurobi commands, cannot see them from the import
    model = read(filename)

    # set some options:
    # time limit in seconds:
    model.params.timeLimit = 60

    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        print('FINISHED: Optimal objective: %g' % model.objVal)
    else:
        print('Optimization ended with status %d' % model.status)

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
    filename = "%s_%d_%d.lp" % (os.path.basename(param.file), param.g1, param.g2)
    if param.unitary:
        dcj_dupindel_ilp_unitary(genomes[param.g1], genomes[param.g2], filename)
    else:
        dcj_dupindel_ilp(genomes[param.g1], genomes[param.g2], filename)

    if param.solve:
        model = solve_ilp(filename)
        model.write(filename + '.sol')
