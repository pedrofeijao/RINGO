__author__ = 'pfeijao'

import glob
import os
import math
from model import Genome, Chromosome
from dendropy import Tree
import ringo_config
import subprocess
import json

LINEAR_END_CHR = "$"
CIRCULAR_END_CHR = ")"

# Ringo_config:
cfg = ringo_config.RingoConfig()

def blossom5_is_available():
  return which(cfg.blossom5_path()) is not None

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def open_mgra_genomes(folder):
    """
    Opens the genomes of the output of a MGRA2 run, where each genome is a .gen file.
    """
    # MGRA creates a folder 'genomes' to store all:
    folder = os.path.join(folder, "genomes")
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


def open_genome_file(filename, as_list=False):
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
    #TODO: add option for genes as strings. I would have to transform them to numbers, and keep a name dictionary
    # to go back when necessary.
    if as_list:
        genomes = []
    else:
        genomes = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split(" ")[0]
                genome = Genome(name)
                if as_list:
                    genomes.append(genome)
                else:
                    genomes[name] = genome
            elif line.startswith("#"):
                # chromosome; ignore name after, each new line is a new chromosome.
                continue
            else:
                if line.endswith(LINEAR_END_CHR):
                    circular = False
                elif line.endswith(CIRCULAR_END_CHR):
                    circular = True
                else:
                    raise RuntimeError("Invalid genome file %s. Unrecognized line:\n%s" % (filename, line))
                genome.add_chromosome(Chromosome(map(int, line[:-1].strip().split(" ")), circular))

    return genomes


def open_copy_number_file(genomes, filename):
    current_idx = 0
    chr_idx = 0
    genome = None
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split(" ")[0]
                genome = genomes[current_idx]
                current_idx += 1
                assert genome.name == name
                chr_idx = 0
            elif line.startswith("#"):
                # chromosome; ignore name after, each new line is a new chromosome.
                continue
            else:
                if line.endswith(LINEAR_END_CHR):
                    circular = False
                elif line.endswith(CIRCULAR_END_CHR):
                    circular = True
                else:
                    raise RuntimeError("Invalid file %s. Unrecognized line:\n%s" % (filename, line))
                genome.chromosomes[chr_idx].copy_number = map(int, line[:-1].strip().split(" "))


def open_genomes_with_copy_number(genome_file, copy_file):
    genomes = open_genome_file(genome_file, as_list=True)
    open_copy_number_file(genomes, copy_file)
    return genomes

def open_adjacencies_file(filename):
    """
    Open a  file in adjacencies format. Genome identifiers are FASTA-like ">name" lines, and
    each following line is an adjacency in the format (x,y), where x and y are integers.
    """
    genome_list = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
              name = line[1:].strip().split(" ")[0]
              genome = []
              genome_list[name] = genome
              continue
            # read adjacency:
            genome.append(eval(line))
    # convert to genomes:
    return {name:Genome.from_adjacency_list(name, adj_list) for name, adj_list in genome_list.iteritems()}



def write_genomes_to_file(genomes, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format.
    """
    if isinstance(genomes, dict):
        iterator = genomes.itervalues()
    elif isinstance(genomes, list):
        iterator = iter(genomes)
    with open(filename, "w") as f:
        for genome in iterator:
            f.write(">%s\n" % genome.name)
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(gene) for gene in chromosome.gene_order]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


def write_genomes_copy_number_to_file(genomes, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format.
    """
    if isinstance(genomes, dict):
        iterator = genomes.itervalues()
    elif isinstance(genomes, list):
        iterator = iter(genomes)
    with open(filename, "w") as f:
        for genome in iterator:
            f.write(">%s\n" % genome.name)
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(gene) for gene in chromosome.copy_number]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


def open_coser_genome(file, name=None):
    chromosomes = {}
    if name is None:
        name = os.path.basename(os.path.splitext(file)[0])
    genome = Genome(name)
    with open(file) as f:
        for l in f:
            g_id, gene, chrom, circular = l.strip().split()
            if chrom not in chromosomes:
                chromosomes[chrom] = Chromosome([], circular=True if circular == 2 else False)
                genome.add_chromosome(chromosomes[chrom])
            chromosomes[chrom].gene_order.append(int(gene))
    return genome



def write_genomes_coser_format(genomes, folder):
    """
    Write genomes in a file with COSER format
    """
    if isinstance(genomes, dict):
        iterator = genomes.itervalues()
    elif isinstance(genomes, list):
        iterator = iter(genomes)

    for genome in iterator:
        with open(os.path.join(folder, "%s.coser" % genome.name), "w") as f:
            for chr_id, chrom in enumerate(genome.chromosomes):
                circular = 2 if chrom.circular else 1
                for gene, copy in zip(chrom.gene_order, chrom.copy_number):
                    f.write("%s_%s\t%d\tchr%d\t%d\n" % (abs(gene), copy, gene, chr_id+1, circular))

def open_newick_tree(filename, label_internal_nodes=True):
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


def write_mgra2_config(extant_genomes, tree, filename):
    """
    Writes a MGRA2 config file, with the leaf genome names and tree topology
    """
    with open(filename, "w") as f:
        f.write("[Genomes]\n")
        for label in extant_genomes.iterkeys():
            f.write("%s Genome_%s\n" % (label, label))
        f.write("\n[Trees]\n")
        f.write(tree.as_string(schema='newick', suppress_rooting=True, suppress_edge_lengths=True))
        f.write("\n")


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
    # read configuration and get filenames:
    singleLeafAdjOut = os.path.join(folder, cfg.declone_output_single_leaf(kT))
    internalWeightsOut = os.path.join(folder, cfg.declone_output_internal_weight(kT))
    extantWeightsOut = os.path.join(folder, cfg.declone_output_extant_weight(kT))

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

def load_ancestral_genomes(folder, method, location):
    #Treat each special case:
    if method.lower() == "mgra":
        if location == "":
            location = cfg.mgra_output_folder()
        reconstructed = open_mgra_genomes(
            os.path.join(folder, location))
    elif method.lower() == "ringo":
        reconstructed = open_genome_file(
            os.path.join(folder, location, cfg.ringo_output_genomes()))
    elif method.lower() == "physca":
        reconstructed = open_adjacencies_file(
            os.path.join(folder, location, cfg.physca_reconstructed_adjacencies()))
    elif method.lower() == "mlgo":
        reconstructed = open_genome_file(
            os.path.join(folder, location, cfg.mlgo_output_genomes()))
    else:  # passing in location the path with genome file name should also work
        reconstructed = open_genome_file(os.path.join(folder, location))

    return reconstructed

# I/O JSON parameters:
def __read_parameters(filename):
  with open(filename, 'r') as f:
      return json.load(f)

def __write_parameters(param, filename):
  with open(filename,"w") as f:
      json.dump(param.__dict__, f, sort_keys = True, indent = 4)

def read_simulation_parameters(folder):
  return __read_parameters(os.path.join(folder, cfg.sim_paramfile()))

def write_simulation_parameters(param, output):
  __write_parameters(param, os.path.join(output,cfg.sim_paramfile()))

def write_ringo_parameters(param, output):
  # filenames to a relative path to output, or absolute:
  outpath = os.path.abspath(param.__dict__["output"])
  for p in ["adj_weights_file","input_genomes","tree", "output"]:
    try:
      abs_p = os.path.abspath(param.__dict__[p])
      common_p = os.path.commonprefix([abs_p, outpath])
      param.__dict__[p] = os.path.relpath(abs_p, outpath)
    # if some path is None (adj_weights_file f.i.):
    except AttributeError:
      param.__dict__[p] = ""
  __write_parameters(param, os.path.join(output,cfg.ringo_output_parameters()))
