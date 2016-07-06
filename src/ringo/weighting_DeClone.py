#!/usr/bin/env python2
import pyximport; pyximport.install()
import argparse

import os
import subprocess
import tempfile
import file_ops
import ConfigParser


nhxFileOut = './nhx_tree'


# convert NEWICK tree notation with weights into nhx notation
# parameter: file - the path to the file with the tree in NEWICK-notation
#           ignore_weights - boolean, if the parsed nhx-tree shall include weights or not
#			set_minimum - str for minimal edge length (default: 0.0)
def parse_NEWICK_to_nhx(file, ignore_weights, nhxOut, minimal_branch_length=0.0):
    print "Parsing Newick into nhx format"
    nhx_tree = ''  # empty tree, will be filled according to nhx-format
    charpos = 0  # current charpos in line
    numberUnamedNodes = 0  # current unamed internal nodes
    numberOfSpecies = -1  # numberOfSpecies is equal to the number of colons
    nameStartPos = 0  # startposition of the last species name
    nameEndPos = 0  # endposition of the last species name
    weightStartPos = 0  # startposition of the last weight
    weightEndPos = 0  # endposition of the lastweight
    lastSpeciesName = ''  # name of the current species
    isLeaf = False  # current species is a leaf or not
    newick_signs = [';', ',', '(', ')', ':']  # signs which have a distinct function in the NEWICK-Format
    isWeight = False
    includesNoWeight = False
    f = open(file, "r")  # file contains only one line
    for line in f:
        # if a : is in the Newick-Tree the tree is very likely weighted
        if ':' not in line:
            includesNoWeight = True
        # for each character
        for char in line:
            if char == '(':
                nameStartPos = charpos  # the current species name begins
                isLeaf = True  # a species directly after a ( is a leaf
                isWeight = False  # after a (, there's never a weight
            elif char == ')':
                weightEndPos = charpos
                if not includesNoWeight:
                    nhx_tree = check_edge_length(minimal_branch_length, nhx_tree, line, weightStartPos, weightEndPos)
                nameEndPos = charpos  # the current species name ends
                if line[charpos - 1] not in newick_signs:  # if there was a species id directly before this ) ...
                    # ...the nhx-block needs to be added
                    nhx_tree += '[&&NHX:D=?:Ev='
                    if isLeaf:  # was the last species a leaf or an internal node?
                        nhx_tree += 'Extant'
                    else:
                        nhx_tree += 'Spec'
                    numberOfSpecies += 1  # the number of species is increased
                    lastSpeciesName = construct_species_name(line, nameStartPos, nameEndPos, numberUnamedNodes)
                    nhx_tree += ':S=' + str(numberOfSpecies) + ':ND=' + str(lastSpeciesName) + ']'
                isLeaf = False  # the next species directly behind a ) isn't a leaf
                nameStartPos = charpos
                isWeight = False  # after a ), there's never directly a weight
            elif char == ',':
                weightEndPos = charpos
                if not includesNoWeight:
                    nhx_tree = check_edge_length(minimal_branch_length, nhx_tree, line, weightStartPos, weightEndPos)
                nameEndPos = charpos  # the current species name ended
                if line[charpos - 1] not in newick_signs:  # if there was a species id directly before this ( ...
                    # ...the nhx-block needs to be added
                    nhx_tree += '[&&NHX:D=?:Ev='
                    if isLeaf:  # was the last species a leaf or an internal node
                        nhx_tree += 'Extant'
                    else:
                        nhx_tree += 'Spec'
                    numberOfSpecies += 1  # the number of species is increased
                    lastSpeciesName = construct_species_name(line, nameStartPos, nameEndPos, numberUnamedNodes)
                    nhx_tree += ':S=' + str(numberOfSpecies) + ':ND=' + str(lastSpeciesName) + ']'
                    isLeaf = True
                    nameStartPos = charpos  # the current species name begins
                    isWeight = False
            elif char == ':':
                weightStartPos = charpos
                # setting temporarily a species name to check if it's empty
                tempSpeciesName = line[nameStartPos + 1:charpos]
                # print 'temp: '+'start: '+str(nameStartPos)+' end: '+str(charpos)+' name: '+str(tempSpeciesName)
                if len(tempSpeciesName) == 0:
                    numberUnamedNodes += 1
                    tempSpeciesName = "node" + str(numberUnamedNodes)
                    nhx_tree += tempSpeciesName
                if ignore_weights:
                    isWeight = True
            elif char == ';':  # the end of the Newick-tree is reached and a root node needs to be added
                if line[charpos - 1] != ')':  # if a identifier for root node is given
                    nameEndPos = charpos
                    lastSpeciesName = construct_species_name(line, nameStartPos, nameEndPos, numberUnamedNodes)
                    nhx_tree += "[&&NHX:D=?:Ev=Spec:S=" + str(numberOfSpecies + 1) + ":ND=" + str(lastSpeciesName) + "]"

                else:  # if no identifier for root node is given
                    nhx_tree += "root[&&NHX:D=?:Ev=Spec:S=" + str(numberOfSpecies + 1) + ":ND=root]"
                isWeight = False
            # if the given Newick -tree doesn't include :, the check for unnamed internal nodes is here
            elif includesNoWeight and line[charpos - 1] == ')':
                # checking if the current internal node is unamed:
                if line[charpos] == ',' or line[charpos] == ')':
                    numberUnamedNodes += 1
                    tempSpeciesName = "node" + str(numberUnamedNodes)
                    nhx_tree += tempSpeciesName

            # if the nhx-Tree shall not include weights, skip adding the current char to the tree, if it belongs to a weight
            if not isWeight:
                nhx_tree += char  # append the current char to the nhx_tree
            charpos += 1
        nhx_tree.strip()
        print nhx_tree
    f.close()
    # write the nhx-tree into a file
    f = open(nhxOut, "w")
    f.write(nhx_tree)
    f.close()


def check_edge_length(minimal_branch_length, nhx_tree, line, weightStartPos, weightEndPos):
    # if current edge length is smaller than demanded minimal edge length
    # replace it with minimal edge length
    weight_str = line[weightStartPos + 1:weightEndPos]
    weight = float(weight_str)
    if weight < minimal_branch_length:
        last_chars = len(weight_str)
        nhx_tree = nhx_tree[:-last_chars]  # removing last characters from nhx_tree
        nhx_tree += str(minimal_branch_length)  # and replacing it with the minimal_branch_length
    return nhx_tree


# constructs the name of the current species and gives it back.
def construct_species_name(line, nameStartPos, nameEndPos, numberUnamedNodes):
    lastSpeciesName = line[nameStartPos + 1:nameEndPos]
    if ':' in lastSpeciesName:
        lastSpeciesName = lastSpeciesName.split(':')[0]
        # if the current species is an unnamed internal node, name it node_number
        if len(lastSpeciesName) == 0:
            lastSpeciesName = "node" + str(numberUnamedNodes)
    return lastSpeciesName


# retrieve extant adjacencies
def readAdjacencyFile(file):
    print "collect extant adjacencies from provided adjacency file"
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    f = open(file, "r")
    for line in f:
        fields = line.rstrip("\n").split("#")
        adjas = fields[0].split(":")
        left = adjas[1].split(" ")[0]
        right = adjas[1].split(" ")[1]
        specieslist = []
        for extant in fields[1:]:
            species = extant.split(":")[0].split(".")[0]
            specieslist.append((species, "mockchrom"))
        specieslist = set(specieslist)
        adjacencies[(left, right)] = specieslist
    return adjacencies


# double marker id to account for orientation
def doubleMarker(marker):
    if "-" in marker:
        return int(marker[1:]) * 2, (int(marker[1:]) * 2) - 1
    else:
        return (int(marker) * 2) - 1, int(marker) * 2


# double marker and find adjacencies in chromosomes
def findAdjacencies(speciesHash):
    print "collect extant adjacencies from marker file..."
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    for species in speciesHash:
        chromosomes = speciesHash[species]
        for chrom in chromosomes:
            gene_order = chromosomes[chrom]
            # telomeres:
            # adj = (0, doubleMarker(orderList[0])[0])
            # for adj_a, adj_b in [(doubleMarker(gene_order[i])[1], doubleMarker(gene_order[i + 1])[0]) for i in range(0, len(gene_order) - 1)]:
            # print gene_order
            # print [(0, doubleMarker(gene_order[-1])[1])]
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

def deCloneProbabilities(extantAdjacencies, kT, listOfInternalNodes, treefile):
    print "Compute probabilities with DeClone..."

    # read DeClone location:
    config = ConfigParser.ConfigParser()
    path = os.path.dirname(os.path.realpath(__file__))
    config.readfp(open(os.path.join(path,'ringo.cfg')))
    declone = config.get("Paths", "declone")

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
        # print output

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
                # print probability
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

# get the internal nodes form the given nhx-treefile
#   -seperatingChar: either ':', if tree includes weights or '[', if not
#   -treefile: the path to the treefile
def get_internal_nodes_from_treefile(seperatingChar, treefile):
    print "Get internal nodes from treefile"
    listOfInternalNodes = []
    # read the nhx-tree
    print "Reading:", treefile
    file = open(treefile, 'r')
    line = file.readline()  # treefile consists only of one line
    file.close()
    # split the tree at each ')', because each internal node's name starts right after one
    line_splitted_at_closed_brackets = line.strip().split(')')
    for element in line_splitted_at_closed_brackets[1:]:
        internal_node = element.strip().split(seperatingChar)[0]
        # in case of the root node the '[&&NHX...'-block needs to be removed
        listOfInternalNodes.append(internal_node.split('[')[0])
    return listOfInternalNodes


def read_Marker_file(marker):
    # compute adjacencies, compute weights with DeClone

    # print "Number of undoubled marker: "
    # read undoubled marker, for each species
    species_marker_order = {}
    file = open(marker, "r")
    species = ""
    for line in file:
        # new species
        if line.startswith(">"):
            if not species == "":
                species_marker_order[species] = chromosomes
            species = line.strip().split("\t")[0][1:]
            chromosomes = {}
            markerCount = 0
            # new chromosome
        elif line.startswith("#"):
            chrom = line.rstrip("\n")[2:]
        elif not line == "\n":
            order = line.rstrip("\n")[:-2].split(" ")
            markerCount = markerCount + len(order)
            chromosomes[chrom] = order
    species_marker_order[species] = chromosomes
    file.close()
    return species_marker_order


def main(newick, ignore_weights, marker_file, kT=0.1, minimum=0.0, write_output=False, output_folder="."):

    nhxFile = os.path.join(output_folder, nhxFileOut)
    parse_NEWICK_to_nhx(newick, ignore_weights, nhxFile, minimum)
    sep = '[' if ignore_weights else ':'
    listOfInternalNodes = get_internal_nodes_from_treefile(sep, nhxFile)
    print "Internal:", listOfInternalNodes
    extantAdjacencies = findAdjacencies(read_Marker_file(marker_file))
    # singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies
    singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies = deCloneProbabilities(extantAdjacencies, kT, listOfInternalNodes, nhxFile)
    if write_output:
        file_ops.write_declone_weights(singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies, kT, folder=output_folder)

    return singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies

if __name__ == '__main__':

    # parsing the input parameters
    parser = argparse.ArgumentParser(
        description='Weights given tree in nhx-format with DeClone. Also converts tree in NEWICK-format into tree in nhx-format')
    parser.add_argument("-nf", '--Newick', required=True, type=str, help='path to the file with NEWICK-tree')
    parser.add_argument("-i", "--ignore_weights", action='store_const', const=True,
                        help="boolean, for either ignore or consider edge length/weights, when parsing Newick Tree into nhx Tree")
    parser.add_argument("-sm", "--set_minimum", type=float,
                        help="minimal value for any edge length, when parsing Newick Tree into nhx Tree", default=0.0)
    parser.add_argument("-m", "--markers", required=True, type=str, help="path to marker-file")
    parser.add_argument("-kT", type=float, help="deClone constant", default=0.1)
    parser.add_argument("-o", "--output_folder", type=str, help="Folder for output files", default=".")
    args = parser.parse_args()

    singleLeafAdj, adjacenciesPerNode, extantWeightedAdjacencies = main(args.Newick, args.ignore_weights, args.markers, kT=args.kT, minimum=args.set_minimum,
         write_output=True, output_folder=args.output_folder)
