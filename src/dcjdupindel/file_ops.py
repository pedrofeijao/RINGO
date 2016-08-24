__author__ = 'pfeijao'

from model import Genome, Chromosome

LINEAR_END_CHR = "$"
CIRCULAR_END_CHR = "@"

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
    1 2 3 4 5 7 @

    """
    genome_list = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split(" ")[0]
                genome = Genome(name)
                genome_list[name] = genome
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
    return genome_list

