#!/usr/bin/env python2
import collections

import math
from dendropy import Tree, Taxon
from model import CType, BPGraph, Genome, Chromosome
import sys
import networkx as nx
import mwmatching
import dcj


# Int Gen:
def ig_indel_small_phylogeny(leaf_genomes, tree, ancestral_adj, solve_median=True):
    # helper functions:
    def add_match_nodes(c_a, c_b, label_a, label_b, ancestral_weight, graph, c_type, path_parity=None):
        # build combined cycle:
        component = collections.deque(reversed(c_a))
        component.extend(c_b)
        # TOdo: fix parity of paths, it' wrong now. should add TWO for even path?
        if c_type == CType.PATH and path_parity == 1:
            # pass
            component.append(0)

        # build adjacency guide:
        component_set = set(component)
        component_weights = {}
        for adj, w in ancestral_weight.iteritems():
            if adj[0] in component_set and adj[1] in component_set:
                component_weights[adj] = w
                # TODO: also check parity of adjacency?
        # solve Max-IS:
        if len(component_weights) > 0:
            max_adj, w = max_weight_ind_set(component, component_weights)

            # create edge, if there is some weight:
            if w > 0:
                graph.add_edge(label_a, label_b, weight=w, max_adj=max_adj, component=component, c_type=c_type)

    def add_path_match_nodes(bp, path_type, label, adj_weight_per_comp, graph):
        odd_path = []
        even_path = []
        for idx, c in enumerate(bp.type_dict[path_type]):
            if len(c) % 2 == 0:
                even_path.append((c, idx))
            else:
                odd_path.append((c, idx))
        # inverted parity results in odd path:
        for i, (path_i, idx_i) in enumerate(even_path):
            for j, (path_j, idx_j) in enumerate(odd_path):
                adj_weights = dict(adj_weight_per_comp[path_type][idx_i])
                adj_weights.update((adj_weight_per_comp[path_type][idx_j]))
                add_match_nodes(path_i, path_j, label + "e%d" % i, label + "o%d" % j, adj_weights, graph,
                                c_type=CType.PATH, path_parity=1)

        # complete the matching if one parity is larger:
        if len(odd_path) == len(even_path):
            comp_list = []
            char = "-"
        elif len(odd_path) > len(even_path):
            comp_list = odd_path
            char = "o"
        else:
            comp_list = even_path
            char = "e"
        # same parity matching, even path:
        for i, (path_i, idx_i) in enumerate(comp_list):
            for j, (path_j, idx_j) in enumerate(comp_list[(i + 1):]):
                j += i + 1
                adj_weights = dict(adj_weight_per_comp[path_type][idx_i])
                adj_weights.update((adj_weight_per_comp[path_type][idx_j]))
                add_match_nodes(path_i, path_j, label + char + "%d" % i, label + char + "%d" % j, adj_weights,
                                graph, c_type=CType.PATH, path_parity=0)
        return comp_list, char

    def reconstruct_node(node):
        # Assuming we have both children of this node are leafs (either originally, or reconstructed).
        node1, node2 = node.child_nodes()

        label = node.label
        print >> sys.stderr, "\nReconstructing %s ..." % label

        # TODO: Leaf adj weights -> DeClone does not use "singleton adj", that are present in just
        # one leaf. Can I use this info?

        # TODO: Logging system.

        # Rebuild node:
        # 1 - find BP graph:
        bp = BPGraph(genomes[node1.label], genomes[node2.label])

        # 2 - adjacencies for this node,  from DeClone:
        ancestral_weight = ancestral_adj[label]

        # Todo: reconstructed_adjacencies is a set. It was original a list, but I was getting
        # repeated adjacencies. HOw can this happen? check;
        reconstructed_adjacencies = set()
        ambiguous_components = list()

        ext_in_comp = {c_type: {} for c_type in CType.all()}

        # Distribute adjacencies by components:
        # first, point each extremity to its component:
        for c_type in CType.all():
            for idx, comp in enumerate(bp.type_dict[c_type]):
                for ext in comp:
                    ext_in_comp[c_type][ext] = idx

        # now "distribute adjacencies:
        adj_weight_per_comp = {c_type: {idx: {} for idx in xrange(len(bp.type_dict[c_type]))} for c_type in CType.all()}
        for adj, w in ancestral_weight.iteritems():
            ext1, ext2 = adj
            for c_type in CType.all():
                if ext1 in ext_in_comp[c_type]:
                    adj_weight_per_comp[c_type][ext_in_comp[c_type][ext1]][adj] = w

        # 3 - Solve the "easy" cases:
        max_adjacencies_from_cycles(bp, reconstructed_adjacencies, ambiguous_components, adj_weight_per_comp)
        # print reconstructed_adjacencies

        # 4 - Find the matching weights and solve max matching:
        comp_match_graph = nx.Graph()

        # AB to AB matching:
        for i, ab_i in enumerate(bp.type_dict[CType.AB_PATH]):
            for j, ab_j in enumerate(bp.type_dict[CType.AB_PATH][(i + 1):]):
                j += i + 1
                adj_weights = dict(adj_weight_per_comp[CType.AB_PATH][i])
                adj_weights.update((adj_weight_per_comp[CType.AB_PATH][j]))
                add_match_nodes(ab_i, ab_j, "AB%d" % i, "AB%d" % j, adj_weights, comp_match_graph,
                                c_type=CType.CYCLE)

        # Ae to Ao matching:
        larger_parity_A = add_path_match_nodes(bp, CType.A_PATH, "A", adj_weight_per_comp, comp_match_graph)

        # Be to Bo matching:
        larger_parity_B = add_path_match_nodes(bp, CType.B_PATH, "B", adj_weight_per_comp, comp_match_graph)

        # SOLVE MAX MATCHING, depending on parity of P_AB:
        # TODO: right now, easy way: just take the unmatched p_a, p_b and p_ab after the regular matching;
        # Todo: because of weight 0, there might be more than one unmatched node on each type;

        # Build edges for the mwmatching script:
        # nodes must be intergers, 0 to n.
        idx_to_v = comp_match_graph.nodes()
        v_to_idx = {v: idx for idx, v in enumerate(idx_to_v)}
        # build edges:
        edges = [(v_to_idx[a], v_to_idx[b], comp_match_graph.get_edge_data(a, b)['weight']) for a, b in
                 comp_match_graph.edges_iter()]
        # solve:
        mw_mate = mwmatching.maxWeightMatching(edges)
        # transform in a dictionary where {a:b} for each (a,b) vertex.
        mate = {idx_to_v[a]: idx_to_v[b] for a, b in enumerate(mw_mate) if b != -1}

        # if odd, we have unmatched AB path('s?)
        # TODO: because of 0 weight edges, we can do this even for even AB.
        if len(bp.type_dict[CType.AB_PATH]) % 2 != 0:
            find_best_triple(bp, mate, adj_weight_per_comp, reconstructed_adjacencies, ambiguous_components)

        # add max weight adjacencies from the best matching:
        seen = set()
        for a, b in mate.iteritems():
            if a in seen:
                continue
            seen.add(b)
            data = comp_match_graph.get_edge_data(a, b)
            component = data['component']
            max_adj = data['max_adj']
            c_type = data['c_type']
            uniq_adj, ambiguous = find_unique_adj(list(component), list(max_adj), c_type)
            max_adj.extend(uniq_adj)
            ambiguous_components.extend(ambiguous)
            reconstructed_adjacencies.update(max_adj)

        s = 1
        for n in [len(comp['c']) / 2 for comp in ambiguous_components]:
            s *= math.factorial(2 * n) // math.factorial(n + 1) // math.factorial(n)
        print "Possible genomes:", s

        # 4.5 = solve IG median to add more adjacencies:
        if solve_median and label != "Root":
            reconstructed_adjacencies, ambiguous_components = adjacencies_from_median(
                label, reconstructed_adjacencies, ambiguous_components, dynamic_tree, genomes)

        # 5 - create genome by using the BP algorithms:
        reconstructed_genome = Genome.from_adjacency_list(label, reconstructed_adjacencies)

        # 6 - add or remove InDel genes:
        # add missing genes:
        gene_set = reconstructed_genome.gene_set()
        for gene in bp.common_AB:
            if gene not in gene_set:
                reconstructed_genome.add_chromosome(Chromosome([gene]))

        # remove indel genes:
        # TODO: this way is "dynamic"; should we do static, from the initial tree?
        ancestral_gene_set = build_ancestral_gene_set(Tree(dynamic_tree), genomes, one_node=node)[node.label]

        no_indel = list(reconstructed_genome.chromosomes)
        for chrom in reconstructed_genome.chromosomes:
            if all([abs(x) in bp.unique_A for x in chrom]) or all([abs(x) in bp.unique_B for x in chrom]):
                # If it is circular, always remove:
                if chrom.circular:
                    no_indel.remove(chrom)
                    continue
                # if linear, use the "guess":
                guess = [abs(g) in ancestral_gene_set for g in chrom]

                # gene_in = [abs(g) in perfect[label].gene_set() for g in chrom]
                # if cmp(guess, gene_in) != 0:
                #     print "Circular:", chrom.circular, "Weights:", guess, "avg:", sum(guess) / len(guess)
                #     print "Gene In ancestral:", gene_in
                #     print "Adj in ancestral:", [adj in perfect[label].adjacency_set() for adj in chrom.adjacency_set()]
                #     print

                # if all are good, continue without removing; else, break the chromosome into
                # parts that are good.
                if all(guess):
                    continue
                no_indel.remove(chrom)
                if any(guess):
                    current_chrom = []
                    for gene_in, gene in zip(guess, chrom):
                        if gene_in:
                            current_chrom.append(gene)
                        else:
                            if len(current_chrom) > 0:
                                no_indel.append(Chromosome(current_chrom))
                                current_chrom = []
                    if len(current_chrom) > 0:
                        no_indel.append(Chromosome(current_chrom))

        # build the new genome after inserting/deleting chromosomes.
        reconstructed_genome = Genome(name=label, chr_list=no_indel)
        # add to the list of genomes:
        genomes[label] = reconstructed_genome

        # Final step - prune dynamic tree:
        # prune the leaves, internal node becomes leaf;
        new_leaf = dynamic_tree.find_node_with_label(label)
        new_leaf.set_child_nodes([])
        new_leaf.taxon = Taxon(label=label)

    def reconstruct_closest_first():
        while True:
            leafs = set(dynamic_tree.leaf_nodes())
            # if only 2 leafs, stop, since we do not reconstruct the root.
            if len(leafs) == 2:
                break
            closest_parent = None
            smallest_dist = 1e6
            while len(leafs) > 0:
                try:
                    node1 = leafs.pop()
                    node2 = node1.sister_nodes()[0]
                    if not node2.is_leaf():
                        continue
                    leafs.remove(node2)
                    if node1.edge.length is None or node2.edge.length is None:
                        dist = dcj.dcj_distance(genomes[node1.label], genomes[node2.label])
                    else:
                        dist = node1.edge.length + node2.edge.length
                    if dist < smallest_dist:
                        smallest_dist = dist
                        closest_parent = node1.parent_node
                except ValueError or TypeError:
                    pass
            reconstruct_node(closest_parent)

    # ==== STARTING MAIN FUNCTION: ig_indel_small_phylogeny

    # dict to store genomes; starts with the leaves, will get the ancestors.
    genomes = dict(leaf_genomes)

    # dynamic tree (gets updated during the process, cutting the leafs of reconstructed nodes.
    dynamic_tree = Tree(tree)

    # bottom up: reconstruct closest leaves first;
    reconstruct_closest_first()
    #
    # at end, return genomes
    return {node.label: genomes[node.label] for node in tree.internal_nodes() if node != tree.seed_node}




def max_adjacencies_from_cycles(bp, reconstructed_adjacencies, ambiguous_components, adj_weight_per_comp):
    # AA and BB paths are self-closed, so basically a regular path:
    for c_type in [CType.AA_PATH, CType.BB_PATH, CType.CYCLE, CType.PATH]:
        for idx, component in enumerate(bp.type_dict[c_type]):
            if len(component) < 2:
                continue
            # only cycle == 2 can be accepted here; paths only in the find_unique recursion.
            if c_type == CType.CYCLE and len(component) == 2:
                reconstructed_adjacencies.add(tuple(sorted(component)))
                continue

            if c_type == CType.PATH:
                # paths need extra "telomere(s)"
                if len(component) % 2 == 0:
                    # component.extend([0, -1])
                    # c_type = CType.CYCLE
                    pass  # TODO: later add TWO telomeres!
                else:
                    pass
                    # c_type = CType.CYCLE
                    # component.append(0)

            component_weights = {}
            add_zero = False
            comp_set = set(component)
            for adj, w in adj_weight_per_comp[c_type][idx].iteritems():
                if c_type == CType.PATH and adj[0] == 0 and adj[1] in comp_set and len(component) % 2 == 1:
                    add_zero = True
                    component_weights[adj] = w
                elif adj[0] in comp_set and adj[1] in comp_set:
                    # if 0 in adj:
                    #     print "WE GOT ONE!", adj, w
                    component_weights[adj] = w
            if add_zero:
                component.append(0)
                c_type = CType.CYCLE
            if len(component_weights) > 0:
                # print "Components weights:", component_weights
                # print "COMP:", component
                max_adj, w = max_weight_ind_set(component, component_weights)

                # for a in max_adj:
                #     if 0 in a:
                #         print "AND IT IS ON THE MAX W SOLUTION!!!!"
                #         tel = a

                # print "Max IS:", max_adj
                # Now i have to check if there are uniquely defined adjacencies
                find_c_type = c_type
                if find_c_type in [CType.AA_PATH, CType.BB_PATH]:
                    find_c_type = CType.CYCLE
                uniq_adj, ambiguous = find_unique_adj(list(component), list(max_adj), find_c_type)
                max_adj.extend(uniq_adj)
                ambiguous_components.extend(ambiguous)

                # print "Max IS extended:", max_adj
                # print "C:", component
                # print "Max_adj:", max_adj
                reconstructed_adjacencies.update(max_adj)
            else:
                ambiguous_components.append({'c': component, 'type': c_type})


def adjacencies_from_median(label, reconstructed_adjacencies, ambiguous_components, current_tree, genomes):
    """
    Find the closest leafs, and complete the ambiguous components by trying
    to minimise distance to this leaf.
    """

    # build current genome to then do a BP:
    reconstructed_genome = Genome.from_adjacency_list(label, reconstructed_adjacencies)
    # closest leafs by traversal:
    search_tree = Tree(current_tree)
    median_node = search_tree.find_node_with_label(label)
    search_tree.reroot_at_node(median_node)
    leaf_distance = {}
    for node in search_tree.preorder_node_iter():
        if node.label == label:
            leaf_distance[node] = 0
        else:
            leaf_distance[node] = node.edge_length + leaf_distance[node.parent_node]
    leafs_by_traversal = sorted([(leaf.label, leaf_distance[leaf]) for leaf in search_tree.leaf_node_iter()],
                                key=lambda x: x[1])

    # closest leafs by DCJ distance:
    # It seems that by traversal works better, specially for higher nodes where DCJ distance starts to be
    # not so precise due to fragmentation and saturation. I will leave here, since this has to be used in
    # case the tree does not have weights.
    # leafs_by_distance = sorted([(leaf.label, dcj.dcj_distance(reconstructed_genome, genomes[leaf.label])) for leaf in
    #                             search_tree.leaf_node_iter()], key=lambda x: x[1])

    # TODO: distance threshold? or just take 3,4 closest genomes:
    # for leaf_label, dcj_dist in leafs_by_distance[:3]:
    for leaf_label, dcj_dist in leafs_by_traversal[:3]:
        # n = len(genomes[leaf.label].gene_set())
        # if dcj_dist > n/2:
        #     break
        bp_outgroup = BPGraph(reconstructed_genome, genomes[leaf_label])
        for p in bp_outgroup.type_dict[CType.PATH]:
            if len(p) % 2 == 0:
                adj = tuple(sorted([p[0], p[-1]]))
                new_ambiguous = []
                for comp in ambiguous_components:
                    new_adj = False
                    if adj[0] in comp['c'] and adj[1] in comp['c']:
                        i = list(comp['c']).index(adj[0])
                        j = list(comp['c']).index(adj[1])
                        if (i - j) % 2 == 1:
                            # ADD ADJ:
                            rec = set([item for sublist in reconstructed_adjacencies for item in sublist])
                            if adj[0] not in rec and adj[1] not in rec:
                                new_adj = True
                                uniq_adj, ambig = find_unique_adj(list(comp['c']), [adj], comp['type'])
                                uniq_adj.append(adj)
                                reconstructed_adjacencies.update(uniq_adj)
                                new_ambiguous.extend(ambig)
                            else:
                                print >> sys.stderr, "CONFLICT!", adj
                                if adj in reconstructed_adjacencies:
                                    print >> sys.stderr, "FULL ADJ."
                    if not new_adj:
                        new_ambiguous.append(comp)
                ambiguous_components = new_ambiguous
    return reconstructed_adjacencies, ambiguous_components

# Tree helpers:
def set_distance_from_node_to_leaves(in_tree):
    for node in in_tree.postorder_node_iter():
        if node.is_leaf():
            node.d_leaf = 0
        else:
            node.d_leaf = min([n.d_leaf + n.edge.length for n in node.child_nodes()])


def set_distance_from_node_to_root(in_tree):
    for node in in_tree.preorder_node_iter():
        if node == in_tree.seed_node:
            node.d_root = 0
        else:
            node.d_root = node.edge.length + node.parent_node.d_root


# best triple:
def find_best_triple(bp, mate, adj_weight_per_comp, reconstructed_adjacencies, ambiguous_components):
    cdef int idx_a, idx_ab, idx_b
    cdef double best_w, w

    def find_unmatched_paths():
      # find unmatched paths:
      cdef int idx_a_e, idx_a_o, idx_b_e, idx_b_o
      idx_a_e = idx_a_o = idx_b_e = idx_b_o = 0
      a_paths = []
      ab_paths = []
      b_paths = []
      for idx_a, p_a in enumerate(bp.type_dict[CType.A_PATH]):
          # find the label of A_path, and if even or odd:
          if len(p_a) % 2 == 0:
              label_a = "Ae%d" % idx_a_e
              idx_a_e += 1
          else:
              label_a = "Ao%d" % idx_a_o
              idx_a_o += 1
          if label_a in mate:
              continue
          a_paths.append((idx_a, p_a))

      for idx_ab, p_ab in enumerate(bp.type_dict[CType.AB_PATH]):
          if "AB%d" % idx_ab in mate:
              continue
          ab_paths.append((idx_ab, p_ab))

      for idx_b, p_b in enumerate(bp.type_dict[CType.B_PATH]):
          # find the label of B_path, and if even or odd:
          if len(p_b) % 2 == 0:
              label_b = "Be%d" % idx_b_e
              idx_b_e += 1
          else:
              label_b = "Bo%d" % idx_b_o
              idx_b_o += 1
          if label_b in mate:
              continue
          b_paths.append((idx_b,p_b))
      return a_paths, ab_paths, b_paths

    best_w = -1
    best_adj = None
    best_component = None

    # unmatched paths:
    a_paths, ab_paths, b_paths = find_unmatched_paths()
    # now make all triples and solve MWIS:
    for idx_a, p_a in a_paths:
        for idx_ab, p_ab in ab_paths:
            for idx_b, p_b in b_paths:
                # calculate the MaxIS for this triple component:

                # build adjacency guide:
                component_weights = {}
                comp_set = set.union(set(p_a),set(p_ab),set(p_b))
                for adj_dict in [adj_weight_per_comp[CType.A_PATH][idx_a], adj_weight_per_comp[CType.B_PATH][idx_b], adj_weight_per_comp[CType.AB_PATH][idx_ab]]:
                    for adj, w in adj_dict.iteritems():
                        if adj[0] in comp_set and adj[1] in comp_set:
                            component_weights[adj] = w

                # solve Max-IS:
                if len(component_weights) > 0:
                    component = list(reversed(p_a))
                    component.extend(p_ab)
                    component.extend(p_b)
                    max_adj, w = max_weight_ind_set(component, component_weights)
                    if w > best_w:
                        best_w = w
                        best_adj = max_adj
                        best_component = component

    # now add the best triple:
    if best_adj is not None:
        uniq_adj, ambiguous = find_unique_adj(list(best_component), list(best_adj), CType.PATH)
        best_adj.extend(uniq_adj)
        ambiguous_components.extend(ambiguous)
        reconstructed_adjacencies.update(best_adj)


# find adjacencies: given
def find_unique_adj(comp, adj_list, c_type):
    # if c_type == CType.PATH and len(comp) % 2 == 1:
    #     c_type = CType.CYCLE
    return find_unique_adj_rec(comp, adj_list, c_type)


def find_unique_adj_rec(comp, adj_list, c_type):
    cdef int idx1, idx2
    if len(comp) < 2:
        return [], []
    if len(comp) == 2:
        if c_type == CType.CYCLE:
            return [tuple(sorted(comp))], []
    while len(adj_list) > 0:
        try:
            adj = adj_list.pop()
            idx1 = comp.index(adj[0])
            idx2 = comp.index(adj[1])
            if idx1 > idx2:
                idx1, idx2 = idx2, idx1

            comp1 = comp[:idx1] + comp[(idx2 + 1):]
            comp2 = comp[(idx1 + 1):idx2]
            unique_adj = []
            ambiguous = []
            if len(comp1) > 0:
                find_adj, find_ambiguous = find_unique_adj_rec(comp1, list(adj_list), c_type)
                unique_adj += find_adj
                ambiguous += find_ambiguous
            if len(comp2) > 0:
                find_adj, find_ambiguous = find_unique_adj_rec(comp2, list(adj_list), CType.CYCLE)
                unique_adj += find_adj
                ambiguous += find_ambiguous

            return unique_adj, ambiguous
        except ValueError:
            pass
    # no more adjacencies; return "ambiguous components"
    return [], [{'c': comp, 'type': c_type}]


def max_weight_ind_set(cycle, weight, parity_filter=True):
    """
    Given a cycle of extremities (deque, list) and weighted adjacencies/edges (dict (tuple):weight),
    find the maximum weight set of noncrossing edges.
    :param cycle:
    :param weight:
    """

    def vertex_idx_and_orient_edges(c, w, p_filter=False):
        length = len(c)
        # dictionary with the index of vertices:
        idx = {c[i]: i for i in range(length)}
        # orient edges and create
        oriented_weight = {}
        for (i, j), w in w.items():
            if p_filter and (idx[i]-idx[j]) % 2 == 0:
                continue
            if idx[i] < idx[j]:
                oriented_weight[(i, j)] = w
            else:
                oriented_weight[(j, i)] = w
        out_edges = collections.defaultdict(list)
        # create out_edges, sorted in order in the cycles
        for edge in sorted(oriented_weight.keys(), key=lambda s: (idx[s[0]], idx[s[1]])):
            out_edges[edge[0]].append(edge[1])
        return idx, out_edges, oriented_weight

    # MAIN:
    cdef int index, i, j, k, s, n, curr_node
    cdef double opt1, opt2
    # idx and orient edges:
    v, out_edges, weight = vertex_idx_and_orient_edges(cycle, weight, parity_filter)

    # create new nodes, for nodes with more than one edge:
    curr_node = max(cycle)+1
    name = {}
    new_cycle = []
    new_cycle_app = new_cycle.append
    for s in cycle:
        if len(out_edges[s]) >= 1:

            for index, j in enumerate(out_edges[s]):
                # create new node
                k = curr_node
                # save the original
                name[k] = s
                curr_node += 1
                out_edges[j].append(k)
                weight[(j, k)] = weight[(s, j)]
                del weight[(s, j)]
                new_cycle_app(k)
            del out_edges[s]
        else:
            new_cycle_app(s)
    cycle = new_cycle
    n = len(cycle)
    # redo: idx and orient edges:
    v, out_edges, weight = vertex_idx_and_orient_edges(cycle, weight, False)

    # structs
    MIS = {}
    T = [None] * (n+1)
    W = {}
    W_MIS = {}
    for edge in sorted(weight.keys(), key=lambda v_s: v[v_s[1]] - v[v_s[0]]):
        MIS[edge] = {edge}
        W_MIS[edge] = weight[edge]
        i, j = v[edge[0]], v[edge[1]]
        if j - i >= 3:
            T[j] = set()
            W[j] = 0
            k = j - 1
            while i < k:
                T[k] = T[k + 1]
                W[k] = W[k+1]
                if len(out_edges[cycle[k]]) > 0:
                    l = v[out_edges[cycle[k]][0]]
                    e_prime = (cycle[k], cycle[l])

                    if l < j and W[k] < W_MIS[e_prime] + W[l + 1]:
                        T[k] = MIS[e_prime].union(T[l + 1])
                        W[k] = W_MIS[e_prime] + W[l + 1]
                k -= 1
            MIS[edge].update(T[i + 1])
            W_MIS[edge] += W[i+1]

    T[n] = set()
    W = [0] * (n + 1)
    for vertex in reversed(cycle):
        i = v[vertex]
        T[i] = set()
        if len(out_edges[vertex]) == 0:
            T[i] = T[i + 1]
            W[i] = W[i + 1]

        else:
            j = v[out_edges[vertex][0]]
            e = (cycle[i], cycle[j])
            opt_1 = W_MIS[e] + W[j + 1]
            opt_2 = W[i + 1]
            if opt_1 > opt_2:
                T[i] = MIS[e].union(T[j + 1])
                W[i] = opt_1
            else:
                T[i] = T[i + 1]
                W[i] = opt_2

    # T[0] has the answer; but we need to go back to "original" vertex names: (because
    # new vertices were added for degree>1 vertices.
    edges2 = [tuple(sorted((name[i], name[j]))) for i, j in T[0]]
    return edges2, W[0]


def build_ancestral_gene_set(tree, leaf_genomes, one_node=None):
    gene_sets = {}
    all_nodes = tree.internal_nodes() if one_node is None else [one_node]
    for ancestral in all_nodes:
        label = ancestral.label
        t = Tree(tree)
        if label != "Root":
            current = t.find_node_with_label(label)
            current.set_child_nodes([])
            t.reroot_at_node(current)
        gene_set_intermediate = {}
        for node in t.postorder_node_iter():
            if node.is_leaf():
                gene_set_intermediate[node.label] = leaf_genomes[node.label].gene_set()
            else:
                gene_set_intermediate[node.label] = set.union(*[gene_set_intermediate[child.label] for child in node.child_nodes()])
        # The set for this current node is on the root:
        gene_sets[label] = gene_set_intermediate[t.seed_node.label]
    return gene_sets


def ancestral_adjacency_weights(tree, leaf_genomes, node=None):
    return ancestral_weights(tree, leaf_genomes, node, genes=False)


def ancestral_gene_weights(tree, leaf_genomes, node=None):
    return ancestral_weights(tree, leaf_genomes, node, genes=True)


def ancestral_weights(tree, leaf_genomes, one_node, genes):
    anc_weights = {}
    all_nodes = tree.internal_nodes() if one_node is None else [one_node]
    for ancestral in all_nodes:
        label = ancestral.label
        t = Tree(tree)
        if ancestral != tree.seed_node:
            current = t.find_node_with_label(label)
            current.set_child_nodes([])
            t.reroot_at_node(current)
        weights = {}
        for node in t.postorder_node_iter():
            if node.is_leaf():
                if genes:
                    weights[node.label] = {adj: 1 for adj in leaf_genomes[node.label].gene_set()}
                else:
                    weights[node.label] = {adj: 1 for adj in leaf_genomes[node.label].adjacency_set()}
                continue
            # weighted sum with the child nodes:
            d = 1e-5 # avoid division by zero in the presence of zero weight edges
            all_adj = set()
            children = list(node.child_nodes())
            for child in children:
                all_adj.update(weights[child.label].iterkeys())
                d += child.edge.length
            # define weight for each adj from children:
            node_adj_weights = {}
            for adj in all_adj:
                children_w = [weights[child.label][adj] * (d - child.edge.length) if adj in weights[child.label] else 0 for child in children]
                node_adj_weights[adj] = sum(children_w)/(d * (len(children)-1))
            weights[node.label] = node_adj_weights

        # The weights for this current node are on the root:
        anc_weights[label] = weights[t.seed_node.label]
    return anc_weights
