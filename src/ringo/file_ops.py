import glob
import os
import math

__author__ = 'pfeijao'
from model import Genome, Chromosome
from dendropy import Tree

LINEAR_END_CHR = "$"
CIRCULAR_END_CHR = ")"

def open_mgra_genomes(folder):
    """
    Opens the genomes of the output of a MGRA2 run, where each genome is a .gen file.
    """
    genomes = {}
    for file in glob.glob(os.path.join(folder, "*.gen")):
        name = os.path.basename(file)[:-4]
        genome = Genome(name)
        genomes[name] = genome
        with open(file) as f:
            for line in f:
                line = line.strip()
                if not line:  # empty
                    continue
                if line.startswith("#"):
                    # chromosome name or comment; ignore
                    continue
                else:
                    genome.add_chromosome(Chromosome(map(int, line[:-1].strip().split(" ")), circular=False))
    return genomes


def open_genome_file(filename):
    """
    Opens a genome file in GRIMM format.
    Example:

    >Genome1
    #chr1 - optional
    1 +2 -3 4 $
    #chr2
    -5 6 $
    >Genome2
    #chr1
    1 2 3 4 5 7 $

    """
    genome_list = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                genome = Genome(line[1:].strip())
                genome_list[genome.name] = genome
            elif line.startswith("#"):
                # chromosome; ignore name after, each new line is a new chromosome.
                continue
            else:
                if line.endswith(LINEAR_END_CHR):
                    circular = False
                elif line.endswith(CIRCULAR_END_CHR):
                    circular = True
                else:
                    raise RuntimeError("Invalid genome file. Unrecognized line:\n%s" % line)
                genome.add_chromosome(Chromosome(map(int, line[:-1].strip().split(" ")), circular))
    return genome_list


def write_genomes_to_file(genomes_dict, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format
    """
    with open(filename, "w") as f:
        for label, genome in genomes_dict.iteritems():
            f.write(">%s\n" % label)
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(x) for x in chromosome.gene_order]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


def open_newick_tree(filename, label_internal_nodes=False):
    """
    Open a tree file in NEWICK format
    """
    t = Tree.get_from_path(filename, schema="newick")
    for node in t.leaf_nodes():
        node.label = node.taxon.label.replace(" ", "_")
    if label_internal_nodes:
        idx = 1
        for node in t.preorder_internal_node_iter():
            if node.label is not None:
                continue
            if node == t.seed_node:
                continue
            else:
                node.label = "M%02d" % idx
                idx += 1
    return t

def write_newick_tree(tree, filename):
    tree.write_to_path(filename, schema="newick")

def write_mgra2_config(leaf_genomes, tree, filename):
    """
    Writes a MGRA2 config file, with the leaf genome names and tree topology
    """
    with open(filename, "w") as f:
        f.write("[Genomes]\n")
        for label in leaf_genomes.iterkeys():
            f.write("%s Genome_%s\n" % (label, label))
        f.write("\n[Trees]\n")
        f.write(tree.as_string(schema='newick', suppress_rooting=True, suppress_edge_lengths=True))
        f.write("\n")


ext_weight_filename = './weighted_extant_adjacencies_%f'
int_weight_filename = './weighted_internal_adjacencies_%f'
single_leaf_filename = './single_leaf_adjacencies_%f'


def write_ancestral_weights(internalWeights, filename):
    """
    Given a dictionary {label:{adj:weight}} (each label is an internal node, pointing to a dictionary of
    adjacencies pointing to weights), write a file where each line is in the format:

    >label  adj weight

    where the fields are tab-separated.
    """
    with open(filename, 'w') as file:
        for node, weights in internalWeights.iteritems():
            for adj, weight in weights.iteritems():  # for each adjacency tuple with weight
                file.write('>%s\t(%d,%d)\t%f\n' % (node, adj[0], adj[1], weight))

def open_ancestral_weights(filename, cutoff=1e-10):
    """
    Open a file with ancestral weights as defined in the function "write_ancestral_weights", store the results
    in a dictionary. Optionally, ignore weights smaller than a given cutoff value.
    """
    w = {}
    with open(filename) as f:
        for l in f:
            genome, adj, weight = l.strip().split()
            weight = float(weight)
            if weight < cutoff:
                continue
            genome = genome[1:]  # cut '>'
            if genome not in w:
                w[genome] = {}
            w[genome][eval(adj)] = weight
    return w

def write_declone_weights(singleLeafAdj, internalWeights, extantWeights, kT, folder="."):
    """
    Write the adjacency weights after running DeClone
    """
    singleLeafAdjOut = os.path.join(folder, single_leaf_filename % kT)
    internalWeightsOut = os.path.join(folder, int_weight_filename % kT)
    extantWeightsOut = os.path.join(folder, ext_weight_filename % kT)

    # ignored adjacencies are written into a special file
    f = open(singleLeafAdjOut, 'w')
    for adj in singleLeafAdj:
        f.write('(' + str(adj[0]) + ',' + str(adj[1]) + ')' + '\t' + str(singleLeafAdj[adj][0]) + '\t' + str(
            singleLeafAdj[adj][1]) + '\n')
    f.close()

    # internal node adjacencies:
    write_ancestral_weights(internalWeights, internalWeightsOut)

    # external leaves
    file = open(extantWeightsOut, 'w')
    for leaf in extantWeights:
        for adj in extantWeights[leaf]:
            file.write('>' + str(leaf) + '\t' + str(adj) + '\n')
    file.close()
