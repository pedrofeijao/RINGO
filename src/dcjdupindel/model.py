from collections import Counter

__author__ = 'pfeijao'


class DupGenome:
    def __init__(self, name="", chr_list=None):
        if not chr_list:
            chr_list = []
        self.name = name
        self.chromosomes = list()
        for gene_order in chr_list:
            if isinstance(gene_order, DupChromosome):
                self.chromosomes.append(gene_order)
            elif isinstance(gene_order, tuple):
                self.chromosomes.append(DupChromosome(list(gene_order), circular=True))
            else:
                self.chromosomes.append(DupChromosome(gene_order, circular=False))


    @staticmethod
    def identity(num_genes, num_chromosomes, name="Identity"):
        genes_per_chr = num_genes / num_chromosomes
        return DupGenome(name,
                         [range(c * genes_per_chr + 1, (c + 1) * genes_per_chr + 1) for c in range(num_chromosomes)])

    def add_chromosome(self, chromosome):
        self.chromosomes.append(chromosome)

    def gene_count(self):
        return Counter([abs(x) for chromosome in self.chromosomes for x in chromosome.gene_order])

    def n_chromosomes(self):
        return len(self.chromosomes)

    def __str__(self):
        return ">" + self.name + "\n" + "[" + ", ".join([str(c) for c in self.chromosomes]) + "]"


class DupChromosome:
    def __init__(self, gene_order, copy_order=None, circular=False):
        self.circular = circular
        self.gene_order = gene_order
        if copy_order is None:
            copy_order = [1] * len(gene_order)
        self.copy_order = copy_order


    def __iter__(self):
        return self.gene_order.__iter__()

    def next(self):
        return self.gene_order.next()

    def __str__(self):
        delimiter = ("(", ")") if self.circular else ("[", "]")
        return "%s%s%s" % (delimiter[0], ", ".join([str(x) for x in self.gene_order]), delimiter[1])

    def length(self):
        return len(self.gene_order)

    def clone(self):
        return DupChromosome(list(self.gene_order), self.circular)


class Ext:
    HEAD = 'h'
    TAIL = 't'

