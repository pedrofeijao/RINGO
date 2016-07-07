#!/usr/bin/env python2
import pyximport; pyximport.install()
import argparse

import os
import subprocess
import tempfile
import file_ops
import ringo_config
from dendropy import Tree


def parse_NEWICK_to_nhx(t, ignore_weights, nhxOut, minimal_branch_length=0.0):
    if minimal_branch_length > 0:
        for edge in t.preorder_edge_iter():
            edge.length = max(minimal_branch_length, edge.length)
    if ignore_weights:
        for edge in t.preorder_edge_iter():
            edge.length = None
    for idx, node in enumerate(t.preorder_node_iter()):
        node.annotations.add_new(name='D', value="?")
        node.annotations.add_new(name='Ev', value="Extant" if node.is_leaf() else "Spec")
        node.annotations.add_new(name='S', value=idx)
        node.annotations.add_new(name='ND', value=node.label)
    t.write_to_path(nhxOut,schema="newick", annotations_as_nhx=True, suppress_annotations=False)


# constructs the name of the current species and gives it back.
def construct_species_name(line, nameStartPos, nameEndPos, numberUnamedNodes):
    lastSpeciesName = line[nameStartPos + 1:nameEndPos]
    if ':' in lastSpeciesName:
        lastSpeciesName = lastSpeciesName.split(':')[0]
        # if the current species is an unnamed internal node, name it node_number
        if len(lastSpeciesName) == 0:
            lastSpeciesName = "node" + str(numberUnamedNodes)
    return lastSpeciesName

# double marker id to account for orientation
def doubleMarker(marker):
    if marker < 0:
        return -marker * 2, -marker * 2 - 1
    else:
        return marker * 2 - 1, marker * 2

# double marker and find adjacencies in chromosomes
def findAdjacencies(speciesHash):
    print "Collecting extant adjacencies from marker file..."
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    for species in speciesHash:
        chromosomes = speciesHash[species]
        for chrom in chromosomes:
            gene_order = chromosomes[chrom]
            # also with telomeres: #TODO: should check for circular chromosomes.
            for adj_a, adj_b in [(0, doubleMarker(gene_order[0])[0])] + \
                    [(doubleMarker(gene_order[i])[1], doubleMarker(gene_order[i + 1])[0])
                     for i in range(0, len(gene_order) - 1)] + [(0, doubleMarker(gene_order[-1])[1])]:
                adj = tuple(sorted((adj_a, adj_b)))
                if adj in adjacencies:
                    adjacencies[adj].append((species, chrom))
                else:
                    adjacencies[adj] = [(species, chrom)]
    return adjacencies

# for each extant adjacency, use declone to compute probability that the
# adjacency is present in an adjacency forest sampled randomly from a
# Boltzmann distribution. Then assign adjacency if it is above threshold to internal node of the tree.

def deCloneProbabilities(declone, extantAdjacencies, kT, listOfInternalNodes, treefile):
    print "Computing probabilities with DeClone..."

    adjacenciesPerNode = {}
    singleLeafAdj = {}
    extantWeightedAdjacencies = {}
    for adj in extantAdjacencies:
        for species in extantAdjacencies[adj]:
            node = species[0]
            if node in extantWeightedAdjacencies:
                extantWeightedAdjacencies[node].add(adj)
            else:
                adjset = set()
                adjset.add(adj)
                extantWeightedAdjacencies[node] = adjset
    for adjacency in extantAdjacencies:
        # produce list of extant adjacencies for declone

        tmpfile = tempfile.NamedTemporaryFile(
            delete=False)  # create a temporary named file (appears with a arbitrary name in the directory)
        species = extantAdjacencies[adjacency]
        for spec in species:
            tmpfile.write(spec[0] + " " + spec[0] + "\n")
        tmpfile.seek(0)  # go to the beginning of the tmpfile
        command = "%s -t1 %s -t2 %s -a %s -i -kT %f" % (declone, treefile, treefile, tmpfile.name, kT)
        # use declone to compute probs
        output = subprocess.check_output(command, shell=True)

        # output is just matrix with probabilities
        # each line of the output should contain max one number greater 0, save for internal nodes
        lines = output.split("\n")
        for line in lines:
            if not line == "" and not line.startswith("\t") and not line.startswith(">"):
                node = line.split("\t")[
                    0]  # for each internal node,find the prob for the current adj to be at this node
                probs = line.split("\t")[1]
                probs = probs.rstrip(" ")
                probArray = probs.split(" ")
                probArrayFl = [float(x) for x in probArray]
                probability = max(probArrayFl, key=float)

                if node in listOfInternalNodes:  # if node is an internal one
                    if len(species) > 1:  # if an adjacency occurs just in one external leaf,
                    # it's ignored for internal nodes (just evolved at this leaf)
                        if node in adjacenciesPerNode:
                            adjacenciesPerNode[node][adjacency] = probability
                        else:
                            adjacenciesPerNode[node] = {adjacency: probability}
                    else:
                        # ignored adjacencies with only one leaf occuring in
                        singleLeafAdj.update({adjacency: (species[0], probability)})

        tmpfile.close()  # tmpfile is closed and immediately deleted

    return singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies

def read_Marker_file(marker):
    # read genomes
    genomes = file_ops.open_genome_file(marker)
    # convert to the format expected for this scrips: TODO: adapt script to use RINGO format directly.
    return {label:{"chr%d" % (idx+1):chrom.gene_order for idx, chrom in enumerate(genome.chromosomes)} for label, genome in genomes.iteritems()}


def main(treefile, ignore_weights, marker_file, kT=0.1, minimum=0.0, write_output=False, output_folder="."):

    # Config:
    cfg = ringo_config.RingoConfig()
    declone = cfg.declone_path()
    # nhx output:
    nhxFile = os.path.splitext(treefile)[0] + ".nhx"
    # open tree:
    t = file_ops.open_newick_tree(treefile)
    # Newick to NHX:
    parse_NEWICK_to_nhx(t, ignore_weights, nhxFile, minimum)
    listOfInternalNodes = [node.label for node in t.internal_nodes()]
    extantAdjacencies = findAdjacencies(read_Marker_file(marker_file))
    singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies = deCloneProbabilities(declone, extantAdjacencies, kT, listOfInternalNodes, nhxFile)
    if write_output:
        file_ops.write_declone_weights(singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies, kT, folder=output_folder)

    return singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies

if __name__ == '__main__':

    # parsing the input parameters
    parser = argparse.ArgumentParser(
        description='Weights given tree in nhx-format with DeClone. Also converts tree in NEWICK-format into tree in nhx-format')
    parser.add_argument("-t", '--tree', required=True, type=str, help='path to the file with NEWICK-tree')
    parser.add_argument("-i", "--ignore_weights", action='store_const', const=True,
                        help="boolean, for either ignore or consider edge length/weights, when parsing Newick Tree into nhx Tree")
    parser.add_argument("-sm", "--set_minimum", type=float,
                        help="minimal value for any edge length, when parsing Newick Tree into nhx Tree", default=0.0)
    parser.add_argument("-m", "--markers", required=True, type=str, help="path to marker-file")
    parser.add_argument("-kT", type=float, help="deClone constant", default=0.1)
    parser.add_argument("-o", "--output_folder", type=str, help="Folder for output files", default=".")
    args = parser.parse_args()

    # run main
    singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies = main(args.tree, args.ignore_weights, args.markers, kT=args.kT, minimum=args.set_minimum,
         write_output=True, output_folder=args.output_folder)
