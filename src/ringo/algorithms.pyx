#!/usr/bin/env python2
import collections

from dendropy import Tree, Taxon
from model import CType, BPGraph, Genome, Chromosome
import sys
import networkx as nx
import mwmatching
import blossom5_perfect_matching
import dcj
import numpy as np
from scipy.optimize import linprog
from numpy.linalg import lstsq
import random
import copy

## Sampling random genomes:

def random_adjacencies_from_cycle(cycle):
    # choose [0,i], where i is odd; cycle length should be even; (odd paths get +1, even get +2 nodes)
    if len(cycle) == 0:
        return []
    assert len(cycle) % 2 == 0
    pos = random.randint(0,(len(cycle)/2)-1) * 2 + 1
    adj = cycle[0], cycle[pos]
    return random_adjacencies_from_cycle(cycle[1:pos])+ [adj] + random_adjacencies_from_cycle(cycle[(pos+1):])

def random_adjacencies_from_cycles(cycle_list):
    return [adj for adjacencies in [random_adjacencies_from_cycle(cycle) for cycle in cycle_list] for adj in adjacencies]

def random_adjacencies_from_list(adj_list):
  assert len(adj_list) % 2 ==0
  random_adj = []
  while len(adj_list)>0:
    # adj = (pos,n-1); pos is random:
    pos = random.randint(0,len(adj_list)-2)
    ex1 = adj_list.pop()
    ex2 = adj_list.pop(pos)
    if ex1 > ex2:
      ex1, ex2 = ex2, ex1
    random_adj.append((ex1,ex2))
  return random_adj



# Int Gen:
def ig_indel_small_phylogeny(leaf_genomes, tree, ancestral_adj, solve_median=True, perfect_matching=False,
    random_repeat=0,add_open_2_cycles=False):
    # helper functions:
    def add_match_nodes(c_a, c_b, label_a, label_b, ancestral_weight, graph, c_type, path_parity=None):
        # build combined cycle:
        component = list(reversed(c_a)) + c_b
        # TOdo: fix parity of paths, Should add TWO zeros for even path and
        # treat that.
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
        else:
            max_adj, w = [], 0
        # create edge, if there is some weight, or if the matching is perfect.
        if w > 0 or perfect_matching:
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
            parity = "-"
        elif len(odd_path) > len(even_path):
            comp_list = odd_path
            parity = "o"
        else:
            comp_list = even_path
            parity = "e"
        # same parity matching, even path:
        for i, (path_i, idx_i) in enumerate(comp_list):
            for j, (path_j, idx_j) in enumerate(comp_list[(i + 1):]):
                j += i + 1
                adj_weights = dict(adj_weight_per_comp[path_type][idx_i])
                adj_weights.update((adj_weight_per_comp[path_type][idx_j]))
                add_match_nodes(path_i, path_j, label + parity + "%d" % i, label + parity + "%d" % j, adj_weights,
                                graph, c_type=CType.PATH, path_parity=0)
        return comp_list, parity

    def build_vertex_edges_for_mwm(comp_match_graph, avoid=None, add_bogus=False, invert_weight=False, scale=1):
        # Build nodes and edges for the mwmatching script/blossum5:
        # output nodes are integers, 0 to n-1.
        # avoid any nodes in a given set 'avoid';
        # 'add bogus' adds an extra vertex with index n, connected to all vertices
        # with weight 0. Usually it is done when n is odd, to allow a perfect matching.
        if avoid is None:
          avoid = set()
        idx_to_v = comp_match_graph.nodes()
        if add_bogus:
          idx_to_v.append("bogus")
        # idx_to_v = sorted(comp_match_graph.nodes())
        v_to_idx = {v: idx for idx, v in enumerate(idx_to_v)}
        # build edges:
        edges = [(v_to_idx[a], v_to_idx[b], -comp_match_graph.get_edge_data(a, b)['weight']*scale if invert_weight else comp_match_graph.get_edge_data(a, b)['weight']*scale) for a, b in comp_match_graph.edges_iter()]
        if add_bogus:
          edges.extend([(v_to_idx[v], v_to_idx["bogus"], 0) for v in idx_to_v[:-1]])
        # debug: sort:
        # edges = sorted([(a,b,w) if a<b else (b,a,w) for a,b,w in edges])
        return idx_to_v, v_to_idx, edges


    def add_adjacencies_from_matching(comp_match_graph, a, b, reconstructed_adjacencies, ambiguous_components):
        data = comp_match_graph.get_edge_data(a, b)
        component = data['component']
        max_adj = data['max_adj']
        c_type = data['c_type']
        uniq_adj, ambiguous = find_unique_adj(
            list(component), list(max_adj), c_type)
        max_adj.extend(uniq_adj)
        ambiguous_components.extend(ambiguous)
        reconstructed_adjacencies.update(max_adj)

    def solve_max_weight_matching(bp, adj_weight_per_comp, reconstructed_adjacencies, ambiguous_components):
        comp_match_graph = {ctype: nx.Graph() for ctype in [CType.AB_PATH, CType.A_PATH, CType.B_PATH]}

        # AB to AB matching:
        for i, ab_i in enumerate(bp.type_dict[CType.AB_PATH]):
            for j, ab_j in enumerate(bp.type_dict[CType.AB_PATH][(i + 1):]):
                j += i + 1
                adj_weights = dict(adj_weight_per_comp[CType.AB_PATH][i])
                adj_weights.update((adj_weight_per_comp[CType.AB_PATH][j]))
                add_match_nodes(ab_i, ab_j, "AB%d" % i, "AB%d" % j, adj_weights, comp_match_graph[CType.AB_PATH],
                                c_type=CType.CYCLE)

        # A matching:
        larger_parity_A_components, parity_A = add_path_match_nodes(
            bp, CType.A_PATH, "A", adj_weight_per_comp, comp_match_graph[CType.A_PATH])

        # B matching:
        larger_parity_B_components, parity_B = add_path_match_nodes(
            bp, CType.B_PATH, "B", adj_weight_per_comp, comp_match_graph[CType.B_PATH])

        # SOLVE MAX MATCHING, depending on parity of P_AB:
        # TODO: right now, easy way: build a matching like it is the even P_AB case, and then
        # just take the unmatched p_a, p_b and p_ab after the regular matching;
        # Because of weight 0, there might be more than one unmatched node on each type; take the
        # best; or maybe I should I try to add more triplets than just the
        # best?

        if perfect_matching:
          # If p_AB is odd, we would have remove every triplet and solve a perfect MWM
          # problem for each triplet, getting the best combination. That is a TODO.
          # As a quicker approach, just add three bogus A_, B_, and AB_ vertices that connects to
          # all triplet candidates;

          for ctype in [CType.AB_PATH, CType.A_PATH, CType.B_PATH]:
            add_bogus = len(bp.type_dict[CType.AB_PATH]) % 2 != 0
            idx_to_v, v_to_idx, edges = build_vertex_edges_for_mwm(comp_match_graph[ctype], add_bogus=add_bogus, invert_weight=True, scale=1000)
            bogus = len(idx_to_v)-1
            if len(edges)>0:
              solution = blossom5_perfect_matching.run_blossom5(len(idx_to_v), len(edges), edges)
              # solution is already
              for a,b in solution:
                # skip bogus edge:
                if add_bogus and b == bogus:
                  continue
                add_adjacencies_from_matching(comp_match_graph[ctype], idx_to_v[a], idx_to_v[b], reconstructed_adjacencies, ambiguous_components)

        # not perfect
        else:
          # Maximum Weight matching, without the restriction of being perfect. In practice,
          # this means that the completion might not be optimal. But, we can get higher weights,
          # which might lead to better reconstruction.
          matched = set()
          for ctype in [CType.AB_PATH, CType.A_PATH, CType.B_PATH]:
            idx_to_v, v_to_idx, edges = build_vertex_edges_for_mwm(comp_match_graph[ctype])
            if len(edges)>0:
              solution = mwmatching.maxWeightMatching(edges)
              # for each edge of the matching, get the adjacencies for the mathed component:
              for a,b in solution:
                matched.add(idx_to_v[a]);matched.add(idx_to_v[b])
                add_adjacencies_from_matching(comp_match_graph[ctype], idx_to_v[a], idx_to_v[b], reconstructed_adjacencies, ambiguous_components)

            # if odd, we have unmatched AB path('s?)
            # TODO: because of 0 weight edges, we can potentially do this also for
            # even AB. For now, only in the odd case we find the best triplet:
            # UPDATE: might have a low rate of TP here, I got around FP/TP = 2,
            # in a simple test, not good. Therefore is it disabled for now.
          # if False:
          if len(bp.type_dict[CType.AB_PATH]) % 2 != 0:
            find_best_triple(bp, matched, adj_weight_per_comp,
                                 reconstructed_adjacencies, ambiguous_components)


    def complete_with_random_adjacencies(genomes, label, tree, reconstructed_adjacencies, ambiguous_components, random_repeats, bp, node):
        # find closest leaf:
        sorted_leafs = closest_leafs_from_node(tree, label)
        closest_label, closest_dist = sorted_leafs[0]
        closest_genome = genomes[closest_label]

        # best:
        best_dist = 1e10
        best_adjs = None
        for i in range(random_repeats):
          cycles = []
          for elem in ambiguous_components:
            if elem['type'] in [CType.PATH, CType.A_PATH, CType.B_PATH]:
              if len(elem['c']) % 2 == 0:
                cycles.append(elem['c'] + [0,0])
              else:
                cycles.append(elem['c'] + [0])
            else:
              cycles.append(elem['c'])

          # cycles = copy.deepcopy(ambiguous_components)
          recon = set(reconstructed_adjacencies)
          # from IntGen only:
          # adjs = random_adjacencies_from_cycles(cycles)
          # or totally random:
          adjs = random_adjacencies_from_list([adj for component in ambiguous_components for adj in component['c']])

          recon.update(adjs)
          g = Genome.from_adjacency_list("Random", recon)
          g = add_remove_singletons(g, bp, tree, genomes, node)

          dcj_distance = dcj.dcj_distance(g, closest_genome)
          if dcj_distance < best_dist:
              best_adjs = adjs
              best_dist = dcj_distance

        return best_adjs

    def add_remove_singletons(reconstructed_genome, bp, tree, genomes, node):
      # add missing genes:
      gene_set = reconstructed_genome.gene_set()
      for gene in bp.common_AB:
          if gene not in gene_set:
              reconstructed_genome.add_chromosome(Chromosome([gene]))

      # remove indel genes:
      # TODO: this way is "dynamic"; should we do static, from the initial
      # tree?
      ancestral_gene_set = build_ancestral_gene_set(
          Tree(tree), genomes, one_node=node)[node.label]

      no_indel = list(reconstructed_genome.chromosomes)
      for chrom in reconstructed_genome.chromosomes:
          if all([abs(x) in bp.unique_A for x in chrom]) or all([abs(x) in bp.unique_B for x in chrom]):
              # If it is circular, always remove:
              if chrom.circular:
                  no_indel.remove(chrom)
                  continue
              # if linear, use the "guess":
              guess = [abs(g) in ancestral_gene_set for g in chrom]

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
      return Genome(name=node.label, chr_list=no_indel)


    def reconstruct_node(node):
        # Assuming we have both children of this node are leafs (either
        # originally, or reconstructed).
        node1, node2 = node.child_nodes()

        label = node.label
        print >> sys.stderr, "Reconstructing %s ..." % label

        # Rebuild node:
        # == find BP graph:
        bp = BPGraph(genomes[node1.label], genomes[node2.label])

        # == adjacency weight for this node:
        ancestral_weight = ancestral_adj[label]

        # Todo: Sometimes I get conflicting adjacencies, how can this happen?
        # check;
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
        adj_weight_per_comp = {c_type: {idx: {} for idx in xrange(
            len(bp.type_dict[c_type]))} for c_type in CType.all()}
        for adj, w in ancestral_weight.iteritems():
            ext1, ext2 = adj
            for c_type in CType.all():
                if ext1 in ext_in_comp[c_type]:
                    adj_weight_per_comp[c_type][
                        ext_in_comp[c_type][ext1]][adj] = w

        # == Solve the "easy" cases:
        max_adjacencies_from_cycles(
            bp, reconstructed_adjacencies, ambiguous_components, adj_weight_per_comp,
             add_open_2_cycles=add_open_2_cycles)

        # == Find the matching weights and solve max matching to find the maximum weight completion,
        # also adding to the reconstructed adjacencies.
        if any([len(bp.type_dict[CType.AB_PATH])>0,len(bp.type_dict[CType.A_PATH])>0,len(bp.type_dict[CType.B_PATH])>0]):
          solve_max_weight_matching(bp, adj_weight_per_comp, reconstructed_adjacencies, ambiguous_components)

        # Count all possible ways of "completing" the genome:
        # s = 1
        # for n in [len(comp['c']) / 2 for comp in ambiguous_components]:
        #     s *= math.factorial(2 * n) // math.factorial(n + 1) // math.factorial(n)
        # print "Possible genomes:", s

        # === solve IG median to add more adjacencies, if not at the root:
        if solve_median and node != dynamic_tree.seed_node:
            reconstructed_adjacencies, ambiguous_components = adjacencies_from_median(
                label, reconstructed_adjacencies, ambiguous_components, dynamic_tree, genomes)

        # === randomly complete:
        if random_repeat > 0:
          new_adj = complete_with_random_adjacencies(genomes, label, dynamic_tree, reconstructed_adjacencies,
                ambiguous_components, random_repeat, bp, node)
          reconstructed_adjacencies.update(new_adj)

        # == solve a Max W Matching with the remaining weights for unmatched adjacencies:
        # TODO: can be a good idea for a more agressive reconstruction, but might get
        # lots of FP.
        if False:
          unmatched = set([ext for component in ambiguous_components for ext in component['c']])
          graph = nx.Graph()
          for adj, w in ancestral_weight.iteritems():
            if w>0.9 and all([(ext in unmatched) for ext in adj]):
              graph.add_edge(adj[0], adj[1], weight=w)

          idx_to_v, v_to_idx, edges = build_vertex_edges_for_mwm(graph)
          if len(edges)>0:
              solution = mwmatching.maxWeightMatching(edges)
              # for each edge of the matching, get the adjacencies for the mathed component:
              for a,b in solution:
                reconstructed_adjacencies.add((idx_to_v[a],idx_to_v[b]))




        # add_adjacencies_from_matching(comp_match_graph[ctype], idx_to_v[a], idx_to_v[b], reconstructed_adjacencies, ambiguous_components)

        # 5 - create genome by using the BP algorithms:
        reconstructed_genome = Genome.from_adjacency_list(label, reconstructed_adjacencies)

        # 6 - add or remove InDel genes:
        # add missing genes, remove singletons:
        reconstructed_genome = add_remove_singletons(reconstructed_genome, bp, dynamic_tree, genomes, node)

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
                        dist = dcj.dcj_distance(
                            genomes[node1.label], genomes[node2.label])
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

    # dynamic tree (gets updated during the process, cutting the leafs of
    # reconstructed nodes.
    dynamic_tree = Tree(tree)

    # bottom up: reconstruct closest leaves first;
    reconstruct_closest_first()
    #
    # at end, return genomes
    return {node.label: genomes[node.label] for node in tree.internal_nodes() if node != tree.seed_node}


def max_adjacencies_from_cycles(bp, reconstructed_adjacencies, ambiguous_components, adj_weight_per_comp, add_open_2_cycles=False):
    # AA and BB paths are self-closed, so basically a regular path:
    for c_type in [CType.AA_PATH, CType.BB_PATH, CType.CYCLE, CType.PATH]:
        for idx, component in enumerate(bp.type_dict[c_type]):
            # if c_type == CType.AA_PATH or c_type == CType.BB_PATH:
            #     c_type = CType.CYCLE
            if len(component) < 2:
                continue
            # only cycle == 2 can be accepted here; paths only in the
            # find_unique recursion.
            if c_type == CType.CYCLE and len(component) == 2:
            # if c_type != CType.PATH and len(component) == 2:
                # print "ADD cycle:",sorted(component)
                reconstructed_adjacencies.add(tuple(sorted(component)))
                continue
            if add_open_2_cycles:
              if c_type in [CType.AA_PATH, CType.BB_PATH] and len(component) == 2:
                  reconstructed_adjacencies.add(tuple(sorted(component)))
                  adj = tuple(sorted(component))
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
                    component_weights[adj] = w
            if add_zero:
                component.append(0)
                c_type = CType.CYCLE
            if len(component_weights) > 0:
                max_adj, w = max_weight_ind_set(component, component_weights)

                # add uniquely defined adjacencies:
                find_c_type = c_type
                if find_c_type in [CType.AA_PATH, CType.BB_PATH]:
                    find_c_type = CType.CYCLE
                uniq_adj, ambiguous = find_unique_adj(
                    list(component), list(max_adj), find_c_type)
                max_adj.extend(uniq_adj)
                ambiguous_components.extend(ambiguous)
                reconstructed_adjacencies.update(max_adj)
            else:
                ambiguous_components.append({'c': component, 'type': c_type})

# Find closest leafs to current internal node:
def closest_leafs_from_node(tree, node_label):
  search_tree = Tree(tree)
  median_node = search_tree.find_node_with_label(node_label)
  search_tree.reroot_at_node(median_node)
  leaf_distance = {}

  # try to use edge lengths; if not available, fall back to DCJ distance:
  # It seems that traversing edge lengths works better, specially for higher nodes where DCJ distance starts to be
  # not so precise due to fragmentation and saturation.
  # UPDATE: now that we have estimated branch lengts at least, just BL then.
  for node in search_tree.preorder_node_iter():
      if node.label == node_label:
          leaf_distance[node] = 0
      else:
          leaf_distance[node] = node.edge_length + \
              leaf_distance[node.parent_node]
  return sorted([(leaf.label, leaf_distance[leaf]) for leaf in search_tree.leaf_node_iter()],
                        key=lambda x: x[1])

def adjacencies_from_median(label, reconstructed_adjacencies, ambiguous_components, current_tree, genomes):
    """
    Find the closest leafs, and complete the ambiguous components by trying
    to minimise distance to these leafs.
    """

    # find closest leafs:
    sorted_leafs = closest_leafs_from_node(current_tree, label)

    # build current genome to then do a BP for each leaf:
    reconstructed_genome = Genome.from_adjacency_list(
        label, reconstructed_adjacencies)

    # TODO: distance threshold? or just take 3,4 closest genomes:
    for leaf_label, dcj_dist in sorted_leafs[:3]:
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
                            rec = set(
                                [item for sublist in reconstructed_adjacencies for item in sublist])
                            if adj[0] not in rec and adj[1] not in rec:
                                new_adj = True
                                uniq_adj, ambig = find_unique_adj(
                                    list(comp['c']), [adj], comp['type'])
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

# Tree functions:
def estimate_branch_lengths_lp(tree, extant_genomes):
  return __estimate_branch_lengths(tree, extant_genomes, "lp")

def estimate_branch_lengths_least_squares(tree, extant_genomes):
  return __estimate_branch_lengths(tree, extant_genomes, "least_squares")


def __estimate_branch_lengths(tree, extant_genomes, method):
    """
    Estimates the branch lenghts of a tree by two availavble
    Solve a min evolution Linear Program based on the DCJ-Indel distances
    between the leaves to
    """
    # index the edges:
    e_idx = 0
    # pairwise paths between leafs:
    pw_paths = {}
    # list of edges in the order of labeling 0,1,...,n-1.
    edges = []
    # paths from the leafs to any internal node. If the internal node is the LCA of two leafs,
    # merging the paths of both leafs gives the PW path.
    leaf_paths = {node.label:{} for node in tree.internal_nodes()}
    # Iterate between all leafs, to create the leaf_paths. If a LCA with previous leafs
    # is found, create the pw_path.
    for node in tree.leaf_node_iter():
        current_node = node
        path = []
        # traverse from leaf to root.
        while current_node != tree.seed_node:
            # if current edge is not labelled, create a label:
            if current_node.edge.label is None:
                current_node.edge.label = e_idx
                e_idx += 1
                edges.append(current_node.edge)
            # current path from leaf to internal node:
            path.append(current_node.edge)
            # update current node
            current_node = current_node.parent_node
            # check if there are LCAs in this internal node with the current leaf:
            for k in leaf_paths[current_node.label].iterkeys():
                # if a PW path exists, this is not LCA, but higher up the tree, so ignore.
                # if the PW does not exist yet, this is a LCA:
                if (node.label, k) not in pw_paths:
                    pw_paths[(node.label, k)] = [p.label for p in leaf_paths[current_node.label][k] + path]
            # store in this internal node the path from the leaf:
            leaf_paths[current_node.label][node.label] = list(path)
    # Now, with all the PW paths, create the LP:
    # min cx
    # s.a. Ax >= b   =>  -Ax <= -b   (the linprog package needs <= type of ineq.)
    # where A is a n x m matrix, n is the number of paths (ineqs.), m the number of edges.
    n = len(extant_genomes)*(len(extant_genomes)-1)/2
    m = e_idx
    A = np.zeros((n,m))
    b = np.zeros(n)
    c = np.ones(m)

    # for each path, fill the correct row of the A matrix and b vector:
    for i, (l, path) in enumerate((pw_paths.iteritems())):
        b[i] = -dcj.dcj_distance(extant_genomes[l[0]], extant_genomes[l[1]])
        for j in path:
            A[i][j] = -1
    if method == "lp":
      # solve the LP:
      result = linprog(c, A_ub=A, b_ub=b).x
    elif method == "least_squares":
      # alternatively: least squares:
      result = lstsq(A, b)[0]
    else:
      print >> sys.stderr, "Unknown method for branch length estimation, skipping..."
      return
    # Apply the lengths in the tree:
    for e, x in zip(edges, result):
        e.length = x
    # the edges from the root are "ambiguous", so each gets the average of the two;
    # from the solution, usually one gets zero and the other the full length;
    node_1, node_2 = tree.seed_node.child_nodes()
    avg = (node_1.edge.length + node_2.edge.length)/2.0
    node_1.edge.length = node_2.edge.length = avg
    # that's it.

def tree_diameter(tree):
  max_leaf_distance = {}
  for node in tree.postorder_node_iter():

    if node == tree.seed_node:
      c1, c2 = node.child_nodes()
      return max_leaf_distance[c1] + c1.edge.length + max_leaf_distance[c2] + c2.edge.length
    if node.is_leaf():
      max_leaf_distance[node] = 0
    else:
      c1, c2 = node.child_nodes()
      max_leaf_distance[node] = max(max_leaf_distance[c1] + c1.edge.length, max_leaf_distance[c2] + c2.edge.length)

def set_all_tree_distances(tree):
  set_min_distance_from_node_to_leaves(tree)
  set_max_distance_from_node_to_leaves(tree)
  set_avg_distance_from_node_to_leaves(tree)
  set_distance_from_node_to_root(tree)
  set_total_evolution_from_node_to_leaves(tree)


def set_min_distance_from_node_to_leaves(in_tree):
    for node in in_tree.postorder_node_iter():
        if node.is_leaf():
            node.min_d_leaf = 0
        else:
            node.min_d_leaf = min([n.min_d_leaf + n.edge.length for n in node.child_nodes()])

def set_max_distance_from_node_to_leaves(in_tree):
    for node in in_tree.postorder_node_iter():
        if node.is_leaf():
            node.max_d_leaf = 0
        else:
            node.max_d_leaf = max([n.max_d_leaf + n.edge.length for n in node.child_nodes()])

def set_avg_distance_from_node_to_leaves(in_tree):
    for node in in_tree.postorder_node_iter():
        if node.is_leaf():
            node.avg_d_leaf = 0
        else:
            node.avg_d_leaf = np.mean([n.max_d_leaf + n.edge.length for n in node.child_nodes()])

def set_total_evolution_from_node_to_leaves(in_tree):
    for node in in_tree.postorder_node_iter():
        if node.is_leaf():
            node.total_ev = 0
        else:
            node.total_ev = sum([n.total_ev + n.edge.length for n in node.child_nodes()])


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
            b_paths.append((idx_b, p_b))
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
                comp_set = set.union(set(p_a), set(p_ab), set(p_b))
                for adj_dict in [adj_weight_per_comp[CType.A_PATH][idx_a], adj_weight_per_comp[CType.B_PATH][idx_b], adj_weight_per_comp[CType.AB_PATH][idx_ab]]:
                    for adj, w in adj_dict.iteritems():
                        if adj[0] in comp_set and adj[1] in comp_set:
                            component_weights[adj] = w

                # solve Max-IS:
                if len(component_weights) > 0:
                    component = list(reversed(p_a))
                    component.extend(p_ab)
                    component.extend(p_b)
                    max_adj, w = max_weight_ind_set(
                        component, component_weights)
                    if w > best_w:
                        best_w = w
                        best_adj = max_adj
                        best_component = component

    # now add the best triple:
    if best_adj is not None:
        uniq_adj, ambiguous = find_unique_adj(
            list(best_component), list(best_adj), CType.PATH)
        best_adj.extend(uniq_adj)
        ambiguous_components.extend(ambiguous)
        reconstructed_adjacencies.update(best_adj)


# find unique adjacencies
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
                find_adj, find_ambiguous = find_unique_adj_rec(
                    comp1, list(adj_list), c_type)
                unique_adj += find_adj
                ambiguous += find_ambiguous
            if len(comp2) > 0:
                find_adj, find_ambiguous = find_unique_adj_rec(
                    comp2, list(adj_list), CType.CYCLE)
                unique_adj += find_adj
                ambiguous += find_ambiguous

            return unique_adj, ambiguous
        except ValueError:
            pass
    # no more adjacencies; return "ambiguous components"
    return [], [{'c': comp, 'type': c_type}]


def max_weight_ind_set(cycle, weight, parity_filter=True):
    """
    Given a cycle of extremities  (list) and weighted adjacencies/edges (dict (tuple):weight),
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
            if p_filter and (idx[i] - idx[j]) % 2 == 0:
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
    v, out_edges, weight = vertex_idx_and_orient_edges(
        cycle, weight, parity_filter)

    # create new nodes, for nodes with more than one edge:
    curr_node = max(cycle) + 1
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
    T = [None] * (n + 1)
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
                W[k] = W[k + 1]
                if len(out_edges[cycle[k]]) > 0:
                    l = v[out_edges[cycle[k]][0]]
                    e_prime = (cycle[k], cycle[l])

                    if l < j and W[k] < W_MIS[e_prime] + W[l + 1]:
                        T[k] = MIS[e_prime].union(T[l + 1])
                        W[k] = W_MIS[e_prime] + W[l + 1]
                k -= 1
            MIS[edge].update(T[i + 1])
            W_MIS[edge] += W[i + 1]

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
    edges = [tuple(sorted((name[i], name[j]))) for i, j in T[0]]
    return edges, W[0]


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
                gene_set_intermediate[node.label] = leaf_genomes[
                    node.label].gene_set()
            else:
                gene_set_intermediate[node.label] = set.union(
                    *[gene_set_intermediate[child.label] for child in node.child_nodes()])
        # The set for this current node is on the root:
        gene_sets[label] = gene_set_intermediate[t.seed_node.label]
    return gene_sets


def ancestral_adjacency_weights(tree, leaf_genomes, node=None):
    return ancestral_weights(tree, leaf_genomes, node, genes=False)


def ancestral_gene_weights(tree, leaf_genomes, node=None):
    return ancestral_weights(tree, leaf_genomes, node, genes=True)


def ancestral_weights(tree, leaf_genomes, one_node, genes):

    # set edge lenghts if not present:
    for edge in tree.preorder_edge_iter():
        if edge.length is None:
            edge.length = 1

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
                    weights[node.label] = {
                        adj: 1 for adj in leaf_genomes[node.label].gene_set()}
                else:
                    weights[node.label] = {
                        adj: 1 for adj in leaf_genomes[node.label].adjacency_set()}
                continue
            # weighted sum with the child nodes:
            d = 0
            all_adj = set()
            children = list(node.child_nodes())
            for child in children:
                all_adj.update(weights[child.label].iterkeys())
                d += child.edge.length
            if d == 0:
                d = 0.1
            # define weight for each adj from children:
            node_adj_weights = {}
            for adj in all_adj:
                children_w = [weights[child.label][
                    adj] * (d - child.edge.length) if adj in weights[child.label] else 0 for child in children]
                node_adj_weights[adj] = sum(
                    children_w) / (d * (len(children) - 1))
            weights[node.label] = node_adj_weights

        # The weights for this current node are on the root:
        anc_weights[label] = weights[t.seed_node.label]
    return anc_weights
