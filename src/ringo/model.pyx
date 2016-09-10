import collections
import copy
from collections import Counter

import itertools

__author__ = 'pfeijao'


class Genome:
    def __init__(self, name="", chr_list=None):
        if not chr_list:
            chr_list = []
        self.name = name
        self.chromosomes = list()
        for gene_order in chr_list:
            if isinstance(gene_order, Chromosome):
                self.chromosomes.append(gene_order)
            elif isinstance(gene_order, tuple):
                self.chromosomes.append(Chromosome(list(gene_order), circular=True))
            else:
                self.chromosomes.append(Chromosome(gene_order, circular=False))

    def n_chromosomes(self):
        return len(self.chromosomes)

    def gene_count(self):
        return Counter([abs(x) for chromosome in self.chromosomes for x in chromosome.gene_order])

    @staticmethod
    #TODO: Check if this is efficient, it uses the "merge" that is not linear?
    # According to the profiler, it is ok, and it is called only about 2*n times.
    def from_adjacency_list(name, reconstructed_adjacencies):
        ext_to_chromosomes = dict()
        # Init components, as simple paths from the adjacencies:
        heads = set()
        for a, b in reconstructed_adjacencies:
            if a == 0 or b == 0:
                continue
            ext_to_chromosomes[a] = {'c': [a, b], 'type': CType.PATH}
            ext_to_chromosomes[b] = ext_to_chromosomes[a]
            if a % 2 == 0:
                heads.add(a)
            else:
                heads.add(a + 1)
            if b % 2 == 0:
                heads.add(b)
            else:
                heads.add(b + 1)
        # add gene
        for h in heads:
            if h - 1 not in ext_to_chromosomes:
                ext_to_chromosomes[h - 1] = {'c': [h - 1], 'type': CType.PATH}
            if h not in ext_to_chromosomes:
                ext_to_chromosomes[h] = {'c': [h], 'type': CType.PATH}

            BPGraph.add_adjacency(ext_to_chromosomes, (h - 1, h))

        # list of unique chromosomes:
        chromosomes = {}
        for component in ext_to_chromosomes.values():
            chromosomes[component['c'][0]] = component['c'], component['type']

        # transform to genes:
        chr_list = []
        for c, t in chromosomes.values():
            gene_c = []
            # "correct" circular chromosome:
            if t == CType.CYCLE:
                c = c[1:] + c[:1]
            for i in xrange(0, len(c), 2):
                if c[i] > c[i + 1]:
                    gene_c.append(-c[i] / 2)
                else:
                    gene_c.append(c[i + 1] / 2)
            if t == CType.CYCLE:
                gene_c = tuple(gene_c)
            chr_list.append(gene_c)

        return Genome(name=name, chr_list=chr_list)

    @staticmethod
    def identity(num_genes, num_chromosomes, name="Identity", circular=False):
        genes_per_chr = num_genes / num_chromosomes
        return Genome(name, [Chromosome(range(c * genes_per_chr + 1, (c + 1) * genes_per_chr + 1), circular=circular)
                      for c in range(num_chromosomes)])


    def add_chromosome(self, chromosome):
        self.chromosomes.append(chromosome)

    def gene_set(self):
        return set([abs(x) for chromosome in self.chromosomes for x in chromosome.gene_order])

    def adjacency_set(self):
        return set.union(*[chromosome.adjacency_set() for chromosome in self.chromosomes])

    def clone(self, name=None):
        # return copy.deepcopy(self)
        g = Genome(name=name if name is not None else self.name)
        for chromosome in self.chromosomes:
            g.add_chromosome(chromosome.clone())
        return g

    def common_adjacencies(self, g):
        return self.adjacency_set().intersection(g.adjacency_set())

    def __str__(self):
        return ">" + self.name + "\n" + "[" + ", ".join([str(c) for c in self.chromosomes]) + "]"

    # Duplicated genes:
    def adjacency_iter_with_copies(self):
        # using the tuplet list generated from "build_extremity_order_lists", outputs a
        # list of pairs of tuplets, represeting the adjacencies.
        ext_order_list = self.build_extremity_order_lists()
        for ext_order in ext_order_list:
            # rotate 1 to make the adjacencies:
            a = iter(ext_order[1:] + ext_order[:1])
            # yield
            for i, j in itertools.izip(a, a):
                yield i, j    # i=(gene_i, copy_i, ext_i), j=(gene_j, copy_j, ext_j)

    def build_extremity_order_lists(self):
        ext_order_list = []
        for idx, chrom in enumerate(self.chromosomes):
            ext_order_list.append(chrom.build_chromosome_ext_order())
        return ext_order_list

    def gene_copies(self):
        copies = collections.defaultdict(set)
        for chrom in self.chromosomes:
            for gene_i, copy_i in zip(chrom.gene_order, chrom.copy_number):
                copies[abs(gene_i)].add(copy_i)
        return copies

class Chromosome:
    def __init__(self, gene_order, circular=False, copy_number=None):
        self.circular = circular
        self.gene_order = gene_order
        self.copy_number = copy_number  # only used in simulations if needed for duplications.

    def __iter__(self):
        return self.gene_order.__iter__()

    def next(self):
        return self.gene_order.next()

    def __str__(self):
        delimiter = ("(", ")") if self.circular else ("[", "]")
        return "%s%s%s" % (delimiter[0], ", ".join([str(x) for x in self.gene_order]), delimiter[1])

    def length(self):
        return len(self.gene_order)

    def adjacency_set(self):
        def adjacency(g1, g2):
            ext1 = 2 * abs(g1) - 1 if g1 < 0 else 2 * abs(g1)
            ext2 = 2 * abs(g2) if g2 < 0 else 2 * abs(g2) - 1
            return (ext1, ext2) if ext1 < ext2 else (ext2, ext1)

        s = {adjacency(gene1, gene2) for gene1, gene2 in zip(self.gene_order, self.gene_order[1:])}
        if self.circular:
            s.add(adjacency(self.gene_order[-1], self.gene_order[0]))
        return s

    def clone(self):
        return Chromosome(list(self.gene_order), self.circular,
                          copy_number=list(self.copy_number) if self.copy_number is not None else None)


    # Duplicated genes:
    def build_chromosome_ext_order(self):
        # returns a list of tuplets (gene, copy, extremity) for the extremities
        ext_list = []
        for gene, copy_n in zip(self.gene_order, self.copy_number):
            if gene >= 0:
                orientation = [Ext.TAIL, Ext.HEAD]
            else:
                orientation = [Ext.HEAD, Ext.TAIL]
            ext_list.extend([(abs(gene), copy_n, ext) for ext in orientation])
        return ext_list

class Ext:
    HEAD = 'h'
    TAIL = 't'


# Breakpoint graph Classes.
class CType:
    PATH = "P"
    CYCLE = "C"
    A_PATH = "P_a"
    B_PATH = "P_b"
    AB_PATH = "P_ab"
    AA_PATH = "P_aa"
    BB_PATH = "P_bb"

    @staticmethod
    def all():
        return [CType.PATH, CType.CYCLE, CType.A_PATH, CType.B_PATH, CType.AB_PATH, CType.AA_PATH, CType.BB_PATH]

    @staticmethod
    def merge(type1, type2):
        # PATH + PATH = PATH
        # PATH + A_PATH = A_PATH
        # PATH + B_PATH = B_PATH
        # A_PATH + A_PATH = AA_PATH
        # A_PATH + B_PATH = AB_PATH
        # B_PATH + B_PATH = BB_PATH
        closed = [CType.AB_PATH, CType.AA_PATH, CType.BB_PATH]
        if type1 in closed or type2 in closed:
            raise RuntimeError("Invalid component merge of type %s and %s" % (type1, type2))
        if type1 == CType.PATH:
            return type2
        if type2 == CType.PATH:
            return type1
        if type1 == type2:
            if type1 == CType.A_PATH:
                return CType.AA_PATH
            if type1 == CType.B_PATH:
                return CType.BB_PATH
            raise RuntimeError("Invalid component merge of type %s and %s" % (type1, type2))
        return CType.AB_PATH

    @staticmethod
    def str(c_type):
        if c_type == CType.PATH:
            return "Path"
        if c_type == CType.A_PATH:
            return "A-Path"
        if c_type == CType.B_PATH:
            return "B-Path"
        if c_type == CType.CYCLE:
            return "Cycle"
        if c_type == CType.AA_PATH:
            return "AA-Path"
        if c_type == CType.BB_PATH:
            return "BB-Path"
        if c_type == CType.AB_PATH:
            return "AB-Path"

class BPGraph:
    # helper functions:
    @staticmethod
    def add_adjacency(components, adjacency):
        ext1, ext2 = adjacency
        c1 = components[ext1]
        c2 = components[ext2]

        # if same component, close a cycle:
        if c1 == c2:
            c1["type"] = CType.CYCLE
            return
        # if not, merge:
        if c1["c"][-1] != ext1:
            c1["c"].reverse()
        if c2["c"][0] != ext2:
            c2["c"].reverse()
        c1["c"].extend(c2["c"])

        # New Type:
        c1["type"] = CType.merge(c1["type"], c2["type"])

        # Combine: TODO not linear, should improve.
        for adj in components[ext2]['c']:
            components[adj] = c1

    @staticmethod
    def other_ext(ext):
        if ext % 2 == 0:
            return ext - 1
        return ext + 1

    def __init__(self, genome_a, genome_b):
        self.genomeA = genome_a
        self.genomeB = genome_b
        self.gene_set_A = genome_a.gene_set()
        self.gene_set_B = genome_b.gene_set()
        self.common_AB = self.gene_set_A.intersection(self.gene_set_B)
        self.unique_A = self.gene_set_A - self.common_AB
        self.unique_B = self.gene_set_B - self.common_AB
        self.type_dict = {t: [] for t in CType.all()}

        all_genes = self.gene_set_A.union(self.gene_set_B)
        # adjacency dict for fast access:
        adj_a = {}
        adj_b = {}
        for ext1, ext2 in self.genomeA.adjacency_set():
            adj_a[ext1] = ext2
            adj_a[ext2] = ext1
        for ext1, ext2 in self.genomeB.adjacency_set():
            adj_b[ext1] = ext2
            adj_b[ext2] = ext1

        # all extremities:
        extremities = set([ext for g in all_genes for ext in (2 * g - 1, 2 * g)])

        # open extremities:
        a_open_extremities = set([ext for g in self.unique_B for ext in (2 * g - 1, 2 * g)])
        b_open_extremities = set([ext for g in self.unique_A for ext in (2 * g - 1, 2 * g)])

        # telomeres:
        a_extremities = set([ext for g in self.gene_set_A for ext in (2 * g - 1, 2 * g)])
        b_extremities = set([ext for g in self.gene_set_B for ext in (2 * g - 1, 2 * g)])
        a_telomeres = {ext for ext in a_extremities if ext not in adj_a}
        b_telomeres = {ext for ext in b_extremities if ext not in adj_b}

        # A-open PATHS:
        while len(a_open_extremities) > 0:
            curr_ext = a_open_extremities.pop()
            comp = [curr_ext]
            while True:
                extremities.remove(curr_ext)
                # Test extremity in B:
                if curr_ext in b_open_extremities:
                    self.type_dict[CType.AB_PATH].append(comp)
                    b_open_extremities.remove(curr_ext)
                    break
                elif curr_ext in b_telomeres:
                    self.type_dict[CType.A_PATH].append(comp)
                    b_telomeres.remove(curr_ext)
                    break

                curr_ext = adj_b[curr_ext]
                comp.append(curr_ext)
                extremities.remove(curr_ext)

                # Now test in A:
                if curr_ext in a_open_extremities:
                    self.type_dict[CType.AA_PATH].append(comp)
                    a_open_extremities.remove(curr_ext)
                    break
                elif curr_ext in a_telomeres:
                    self.type_dict[CType.A_PATH].append(comp)
                    a_telomeres.remove(curr_ext)
                    break
                curr_ext = adj_a[curr_ext]
                comp.append(curr_ext)

        # B-open PATHS:
        while len(b_open_extremities) > 0:
            curr_ext = b_open_extremities.pop()
            comp = [curr_ext]
            while True:
                extremities.remove(curr_ext)
                # test if it is adj in A: (no need to test A-open, all treated by now).
                if curr_ext in a_telomeres:
                    self.type_dict[CType.B_PATH].append(comp)
                    a_telomeres.remove(curr_ext)
                    break
                curr_ext = adj_a[curr_ext]
                comp.append(curr_ext)
                extremities.remove(curr_ext)
                # Test extremity in B:
                if curr_ext in b_open_extremities:
                    self.type_dict[CType.BB_PATH].append(comp)
                    b_open_extremities.remove(curr_ext)
                    break
                elif curr_ext in b_telomeres:
                    self.type_dict[CType.B_PATH].append(comp)
                    b_telomeres.remove(curr_ext)
                    break

                curr_ext = adj_b[curr_ext]
                comp.append(curr_ext)

        # A-telomere paths:

        while len(a_telomeres) > 0:
            curr_ext = a_telomeres.pop()
            comp = [curr_ext]
            while True:
                extremities.remove(curr_ext)
                if curr_ext in b_telomeres:
                    self.type_dict[CType.PATH].append(comp)
                    b_telomeres.remove(curr_ext)
                    break
                else:
                    curr_ext = adj_b[curr_ext]
                    comp.append(curr_ext)

                extremities.remove(curr_ext)
                if curr_ext in a_telomeres:
                    self.type_dict[CType.PATH].append(comp)
                    a_telomeres.remove(curr_ext)
                    break
                else:
                    curr_ext = adj_a[curr_ext]
                    comp.append(curr_ext)

        # B-telomere paths:
        while len(b_telomeres) > 0:
            curr_ext = b_telomeres.pop()
            comp = [curr_ext]
            while True:
                extremities.remove(curr_ext)
                curr_ext = adj_a[curr_ext]
                comp.append(curr_ext)
                extremities.remove(curr_ext)
                if curr_ext in b_telomeres:
                    self.type_dict[CType.PATH].append(comp)
                    b_telomeres.remove(curr_ext)
                    break
                curr_ext = adj_b[curr_ext]
                comp.append(curr_ext)

        # cycles:
        while len(extremities) > 0:
            st = extremities.pop()
            curr_ext = st
            comp = [curr_ext]
            while True:
                curr_ext = adj_a[curr_ext]
                comp.append(curr_ext)
                extremities.remove(curr_ext)
                curr_ext = adj_b[curr_ext]
                if curr_ext == st:
                    self.type_dict[CType.CYCLE].append(comp)
                    break
                comp.append(curr_ext)
                extremities.remove(curr_ext)

    @staticmethod
    def int2ext(x):
        if x % 2 == 0:
            return "%dh" % (x / 2)
        return "%dt" % ((x + 1) / 2)

    def __str__(self):

        s = ""
        for t, l in self.type_dict.items():
            if len(l) == 0:
                continue
            s += t + ":\n"
            for c in l:
                # s += str([int2ext(x) for x in c]) + ", "
                s += str(["%d (%s)" % (x, BPGraph.int2ext(x)) for x in c]) + ", "
            s = s[:-2] + "\n"
        return s[:-1]
