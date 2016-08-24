#!/usr/bin/env python2
import os
import sys
import itertools
import operator
import file_ops
from gurobipy import *

class Ext:
    HEAD = 'h'
    TAIL = 't'


def ilp(genome_a, genome_b):
    # def vertex_name(genome, gene, copy, is_tail):
    #     return "%s%s_%s%s" % (genome, gene, copy, "t" if is_tail else "h")

    def vertex_name(genome, gene, copy, ext):
        return "%s%s_%s%s" % (genome, gene, copy, ext)

    # def matching_edge_name(gene, copyA, copyB, is_tail):
    #     return "x_%s_%s" % (vertex_name("A", gene, copyA, is_tail), vertex_name("B", gene, copyB, is_tail))
    def matching_edge_name(gene, copyA, copyB, ext):
        return "x_%s,%s" % (vertex_name("A", gene, copyA, ext), vertex_name("B", gene, copyB, ext))

    def balancing_edge_name(genome, gene1, copy1, ext1, gene2, copy2, ext2):
        if (gene1, copy1, ext1) > (gene2, copy2, ext2):
            gene1, copy1, ext1, gene2, copy2, ext2 = gene2, copy2, ext2, gene1, copy1, ext1
        return "w_%s,%s" % (vertex_name(genome, gene1, copy1, ext1), vertex_name(genome, gene2, copy2, ext2))

    def balancing_extremities(balancing):
        for gene, copies in balancing.iteritems():
            for copy in copies:
                yield gene, copy, Ext.HEAD
                yield gene, copy, Ext.TAIL

    # count of each gene on each genome
    gene_count_a = genome_a.gene_count()
    gene_count_b = genome_b.gene_count()

    # find all genes
    gene_count = {g: max(gene_count_a[g], gene_count_b[g]) for g in
                  set(gene_count_a.keys()).union(set(gene_count_b.keys()))}

    # build constraints:
    constraints = []
    # consistency:
    for gene, copies in gene_count.iteritems():
        constraints.append("\\ Match gene %s" % gene)
        for i, j in itertools.product(xrange(1, copies + 1), xrange(1, copies + 1)):
            constraints.append(
                "%s - %s = 0" % (matching_edge_name(gene, i, j, Ext.TAIL), matching_edge_name(gene, i, j, Ext.HEAD)))

    # unique matching: (1-to-1) (only tail, it should suffice, together with the consistency)
    constraints.append("\\ Unique matching:")
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, copies + 1):
            constraints.append(
                " + ".join([matching_edge_name(gene, i, j, Ext.TAIL) for j in xrange(1, copies + 1)]) + " = 1")

    constraints.append("\\ Telomere matching:")  # special case, "0" gene.
    n_telomeres = 2 * max(genome_a.n_chromosomes(), genome_b.n_chromosomes())
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
    label = {}
    idx = 1
    for genome in ["A", "B"]:
        # telomeres: (represented by TAIL always)
        for i in range(1, n_telomeres + 1):
            label[vertex_name(genome, 0, i, Ext.TAIL)] = idx
            idx += 1
        for gene, copies in gene_count.iteritems():
            for i in xrange(1, copies + 1):
                label[vertex_name(genome, gene, i, Ext.HEAD)] = idx
                idx += 1
                label[vertex_name(genome, gene, i, Ext.TAIL)] = idx
                idx += 1

    # for each adjacency, fix label:
    constraints.append("\\ Adjacency have the same label:")
    for genome, genome_name in [(genome_a, "A"), (genome_b, "B")]:
        copy_number = {gene: 1 for gene in gene_count.iterkeys()}
        copy_number[0] = 1
        for idx, chrom in enumerate(genome.chromosomes):
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
                constraints.append("y_%s - y_%s = 0 \\ %s <-> %s " % (label[v_i], label[v_j], v_i, v_j))

    constraints.append("\\ Matching edges with the same label:")
    for gene, copies in gene_count.iteritems():
        for i in xrange(1, copies + 1):
            for j in xrange(1, copies + 1):
                constraints.append("\\ Edges (%s%s_%s, %s%s_%s)" % ("A", gene, i, "B", gene, j))
                for ext in [Ext.HEAD, Ext.TAIL]:
                    y_i = label[vertex_name("A", gene, i, ext)]
                    y_j = label[vertex_name("B", gene, j, ext)]
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
                y_i = label[vertex_name(genome, gene_i, copy_i, ext_i)]
                y_j = label[vertex_name(genome, gene_j, copy_j, ext_j)]
                constraints.append("y_%s - y_%s + %s %s <= %d" % (
                    y_i, y_j, y_i, balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j), y_i))
                constraints.append("y_%s - y_%s + %s %s <= %d" % (
                    y_j, y_i, y_j, balancing_edge_name(genome, gene_i, copy_i, ext_i, gene_j, copy_j, ext_j), y_j))

    # z variables: since all cycles have to contains vertices from both genomes, we only add z variables
    # for genome A, that have smallest labels, so a genome B z variable will never be =1.
    constraints.append("\\ Z variables")
    for vertex, i in sorted(label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "A":
            constraints.append("%d z_%s - y_%s <= 0" % (i, i, i))

    # number of genes, to fix distance:
    constraints.append("n = %d" % (sum(gene_count.itervalues()) + n_telomeres/2))
    bounds = []
    for i in sorted(label.itervalues()):
        bounds.append("y_%d <= %d" % (i, i))

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
    for vertex, i in sorted(label.items(), key=operator.itemgetter(1)):
        if vertex[0] == "B":
            continue
        binary.append("z_%d" % i)

    # Y label are general:
    general = []
    for vertex, i in sorted(label.items(), key=operator.itemgetter(1)):
        general.append("y_%d" % i)

    # number of genes:
    general.append("n")
    # objective function:
    objective = ["obj: n - " + " - ".join(["z_%d" % i for vertex, i in sorted(label.items(), key=operator.itemgetter(1)) if vertex[0] != "B"])]


    # PRINT:
    output = "%s.lp" % os.path.basename(sys.argv[1])
    with open(output, "w") as f:
        for header, lines in [("Minimize", objective), ("Subject to", constraints), ("Bounds", bounds), ("Binary", binary), ("General", general)]:
            print >> f, header
            print >> f, "\n".join(lines)

    # pycharm complains of gurobi commands, cannot see them from the import
    model = read(output)
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        print('FINISHED: Optimal objective: %g' % model.objVal)
    else:
        print('Optimization ended with status %d' % model.status)
        exit(0)

    model.params.outputFlag = 0
    print('')
    for k in range(model.solCount):
        model.params.solutionNumber = k
        objn = 0
        for v in model.getVars():
            objn += v.obj * v.xn
        print('Solution %d has objective %g' % (k, objn))
    print('')
    model.params.outputFlag = 1

    fixed = model.fixed()
    fixed.params.presolve = 0
    fixed.optimize()

    if fixed.status != GRB.Status.OPTIMAL:
        print("Error: fixed model isn't optimal")
        exit(1)

    diff = model.objVal - fixed.objVal

    if abs(diff) > 1e-6 * (1.0 + abs(model.objVal)):
        print('Error: objective values are different')
        exit(1)

    # Print values of nonzero variables
    # for v in fixed.getVars():
    #     if v.x != 0:
    #         if v.varName[0] == "w":
    #             print('%s %g' % (v.varName, v.x))


if __name__ == '__main__':
    g = file_ops.open_genome_file(sys.argv[1])
    genomes = g.values()
    ilp(genomes[0], genomes[1])
