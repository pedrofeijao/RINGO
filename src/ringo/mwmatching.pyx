"""Weighted maximum matching in general graphs.

The algorithm is taken from "Efficient Algorithms for Finding Maximum
Matching in Graphs" by Zvi Galil, ACM Computing Surveys, 1986.
It is based on the "blossom" method for finding augmenting paths and
the "primal-dual" method for finding a matching of maximum weight, both
due to Jack Edmonds.
Some ideas came from "Implementation of algorithms for maximum matching
on non-bipartite graphs" by H.J. Gabow, Standford Ph.D. thesis, 1973.

A C program for maximum weight matching by Ed Rothberg was used extensively
to validate this new code.
"""

#
# Changes:
#
# 2013-04-07
#   * Added Python 3 compatibility with contributions from Daniel Saunders.
#
# 2008-06-08
#   * First release.
#


def maxWeightMatching(edges, maxcardinality=False):
    """Compute a maximum-weighted matching in the general undirected
    weighted graph given by "edges".  If "maxcardinality" is true,
    only maximum-cardinality matchings are considered as solutions.

    Edges is a sequence of tuples (i, j, wt) describing an undirected
    edge between vertex i and vertex j with weight wt.  There is at most
    one edge between any two vertices; no vertex has an edge to itself.
    Vertices are identified by consecutive, non-negative integers.

    Return a list "mate", such that mate[i] == j if vertex i is
    matched to vertex j, and mate[i] == -1 if vertex i is not matched.

    This function takes time O(n ** 3)."""

    #
    # Vertices are numbered 0 .. (nvertex-1).
    # Non-trivial blossoms are numbered nvertex .. (2*nvertex-1)
    #
    # Edges are numbered 0 .. (nedge-1).
    # Edge endpoints are numbered 0 .. (2*nedge-1), such that endpoints
    # (2*k) and (2*k+1) both belong to edge k.
    #
    # Many terms used in the comments (sub-blossom, T-vertex) come from
    # the paper by Galil; read the paper before reading this code.
    #

    integer_types = (int, long)

    # Deal swiftly with empty graphs.
    if not edges:
        return []

    # cdefs:
    cdef int nedge, n_vertex, i, j, w, v, k, p, t, b
    cdef double wt, maxweight, kslack

    # Count vertices.
    nedge = len(edges)
    n_vertex = 0

    for (i, j, wt) in edges:
        if i >= n_vertex:
            n_vertex = i + 1
        if j >= n_vertex:
            n_vertex = j + 1

    # Find the maximum edge weight.
    maxweight = max(0, max([wt for (i, j, wt) in edges]))

    # If p is an edge endpoint,
    # endpoint[p] is the vertex to which endpoint p is attached.
    # Not modified by the algorithm.
    # cdef int[:] endpoint
    endpoint = [edges[p // 2][p % 2] for p in range(2 * nedge)]

    # If v is a vertex,
    # neighbend[v] is the list of remote endpoints of the edges attached to v.
    # Not modified by the algorithm.
    neighbend = [[] for i in range(n_vertex)]

    for k in range(len(edges)):
        (i, j, w) = edges[k]
        neighbend[i].append(2*k + 1)
        neighbend[j].append(2*k)

    # If v is a vertex,
    # mate[v] is the remote endpoint of its matched edge, or -1 if it is single
    # (i.e. endpoint[mate[v]] is v's partner vertex).
    # Initially all vertices are single; updated during augmentation.
    # cdef array.array ar_mate = array.array(nvertex * [-1])
    # cdef int[:] mate = ar_mate
    mate = n_vertex * [-1]

    # If b is a top-level blossom,
    # label[b] is 0 if b is unlabeled (free);
    #             1 if b is an S-vertex/blossom;
    #             2 if b is a T-vertex/blossom.
    # The label of a vertex is found by looking at the label of its
    # top-level containing blossom.
    # If v is a vertex inside a T-blossom,
    # label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
    # Labels are assigned during a stage and reset after each augmentation.
    label = (2 * n_vertex) * [0]

    # If b is a labeled top-level blossom,
    # labelend[b] is the remote endpoint of the edge through which b obtained
    # its label, or -1 if b's base vertex is single.
    # If v is a vertex inside a T-blossom and label[v] == 2,
    # labelend[v] is the remote endpoint of the edge through which v is
    # reachable from outside the blossom.
    labelend = (2 * n_vertex) * [-1]

    # If v is a vertex,
    # inblossom[v] is the top-level blossom to which v belongs.
    # If v is a top-level vertex, v is itself a blossom (a trivial blossom)
    # and inblossom[v] == v.
    # Initially all vertices are top-level trivial blossoms.
    inblossom = list(xrange(n_vertex))

    # If b is a sub-blossom,
    # blossomparent[b] is its immediate parent (sub-)blossom.
    # If b is a top-level blossom, blossomparent[b] is -1.
    blossomparent = (2 * n_vertex) * [-1]

    # If b is a non-trivial (sub-)blossom,
    # blossomchilds[b] is an ordered list of its sub-blossoms, starting with
    # the base and going round the blossom.
    blossomchilds = (2 * n_vertex) * [None]

    # If b is a (sub-)blossom,
    # blossombase[b] is its base VERTEX (i.e. recursive sub-blossom).
    blossombase = list(range(n_vertex)) + n_vertex * [-1]

    # If b is a non-trivial (sub-)blossom,
    # blossomendps[b] is a list of endpoints on its connecting edges,
    # such that blossomendps[b][i] is the local endpoint of blossomchilds[b][i]
    # on the edge that connects it to blossomchilds[b][wrap(i+1)].
    blossomendps = (2 * n_vertex) * [None]

    # If v is a free vertex (or an unreached vertex inside a T-blossom),
    # bestedge[v] is the edge to an S-vertex with least slack,
    # or -1 if there is no such edge.
    # If b is a (possibly trivial) top-level S-blossom,
    # bestedge[b] is the least-slack edge to a different S-blossom,
    # or -1 if there is no such edge.
    # This is used for efficient computation of delta2 and delta3.
    bestedge = (2 * n_vertex) * [-1]

    # If b is a non-trivial top-level S-blossom,
    # blossombestedges[b] is a list of least-slack edges to neighbouring
    # S-blossoms, or None if no such list has been computed yet.
    # This is used for efficient computation of delta3.
    blossombestedges = (2 * n_vertex) * [None]

    # List of currently unused blossom numbers.
    unusedblossoms = list(xrange(n_vertex, 2 * n_vertex))

    # If v is a vertex,
    # dualvar[v] = 2 * u(v) where u(v) is the v's variable in the dual
    # optimization problem (multiplication by two ensures integer values
    # throughout the algorithm if all edge weights are integers).
    # If b is a non-trivial blossom,
    # dualvar[b] = z(b) where z(b) is b's variable in the dual optimization
    # problem.
    dualvar = n_vertex * [maxweight] + n_vertex * [0]

    # If allowedge[k] is true, edge k has zero slack in the optimization
    # problem; if allowedge[k] is false, the edge's slack may or may not
    # be zero.
    allowedge = nedge * [False]

    # Queue of newly discovered S-vertices.
    queue = []

    # Return 2 * slack of edge k (does not work inside blossoms).
    def slack(int v_k):
        cdef int v_i, v_j
        cdef double w_t
        (v_i, v_j, w_t) = edges[v_k]
        return dualvar[v_i] + dualvar[v_j] - 2 * w_t

    # Generate the leaf vertices of a blossom.
    def blossom_leaves(int bl):
        cdef int v_t, v_v
        if bl < n_vertex:
            yield bl
        else:
            for v_t in blossomchilds[bl]:
                if v_t < n_vertex:
                    yield v_t
                else:
                    for v_v in blossom_leaves(v_t):
                        yield v_v

    # Assign label t to the top-level blossom containing vertex w
    # and record the fact that w was reached through the edge with
    # remote endpoint p.
    # def assignLabel(int v_w, t, int p):
    def assign_label(int v_w, int v_t, int v_p):
        cdef int v_v, bl, v_base
        bl = inblossom[v_w]
        assert label[v_w] == 0 and label[bl] == 0
        label[v_w] = label[bl] = v_t
        labelend[v_w] = labelend[bl] = v_p
        bestedge[v_w] = bestedge[bl] = -1
        if v_t == 1:
            # b became an S-vertex/blossom; add it(s vertices) to the queue.
            queue.extend(blossom_leaves(bl))
        elif v_t == 2:
            # b became a T-vertex/blossom; assign label S to its mate.
            # (If b is a non-trivial blossom, its base is the only vertex
            # with an external mate.)
            v_base = blossombase[bl]
            assert mate[v_base] >= 0
            assign_label(endpoint[mate[v_base]], 1, mate[v_base] ^ 1)

    # Trace back from vertices v and w to discover either a new blossom
    # or an augmenting path. Return the base vertex of the new blossom or -1.
    def scan_blossom(int v_v, int v_w):
        # Trace back from v and w, placing breadcrumbs as we go.
        path = []
        cdef int v_base, bl
        v_base = -1
        while v_v != -1 or v_w != -1:
            # Look for a breadcrumb in v's blossom or put a new breadcrumb.
            bl = inblossom[v_v]
            if label[bl] & 4:
                v_base = blossombase[bl]
                break
            assert label[bl] == 1
            path.append(bl)
            label[bl] = 5
            # Trace one step back.
            assert labelend[bl] == mate[blossombase[bl]]
            if labelend[bl] == -1:
                # The base of blossom b is single; stop tracing this path.
                v_v = -1
            else:
                v_v = endpoint[labelend[bl]]
                bl = inblossom[v_v]
                assert label[bl] == 2
                # b is a T-blossom; trace one more step back.
                assert labelend[bl] >= 0
                v_v = endpoint[labelend[bl]]
            # Swap v and w so that we alternate between both paths.
            if v_w != -1:
                v_v, v_w = v_w, v_v
                # Remove breadcrumbs.
        for bl in path:
            label[bl] = 1
        # Return base vertex, if we found one.
        return v_base

    # Construct a new blossom with given base, containing edge k which
    # connects a pair of S vertices. Label the new blossom as S; set its dual
    # variable to zero; relabel its T-vertices to S and add them to the queue.
    def add_blossom(int v_base, int v_k):
        cdef int v_i, v_j, v_v, v_w, bb,bv,bw,bl
        cdef double v_wt
        (v_v, v_w, v_wt) = edges[v_k]
        bb = inblossom[v_base]
        bv = inblossom[v_v]
        bw = inblossom[v_w]
        # Create blossom.
        bl = unusedblossoms.pop()
        blossombase[bl] = v_base
        blossomparent[bl] = -1
        blossomparent[bb] = bl
        # Make list of sub-blossoms and their interconnecting edge endpoints.
        blossomchilds[bl] = path = []
        blossomendps[bl] = end_p = []
        # Trace back from v to base.
        while bv != bb:
            # Add bv to the new blossom.
            blossomparent[bv] = bl
            path.append(bv)
            end_p.append(labelend[bv])
            assert (label[bv] == 2 or
                    (label[bv] == 1 and labelend[bv] == mate[blossombase[bv]]))
            # Trace one step back.
            assert labelend[bv] >= 0
            v_v = endpoint[labelend[bv]]
            bv = inblossom[v_v]
        # Reverse lists, add endpoint that connects the pair of S vertices.
        path.append(bb)
        path.reverse()
        end_p.reverse()
        end_p.append(2 * v_k)
        # Trace back from w to base.
        while bw != bb:
            # Add bw to the new blossom.
            blossomparent[bw] = bl
            path.append(bw)
            end_p.append(labelend[bw] ^ 1)
            assert (label[bw] == 2 or
                    (label[bw] == 1 and labelend[bw] == mate[blossombase[bw]]))
            # Trace one step back.
            assert labelend[bw] >= 0
            v_w = endpoint[labelend[bw]]
            bw = inblossom[v_w]
        # Set label to S.
        assert label[bb] == 1
        label[bl] = 1
        labelend[bl] = labelend[bb]
        # Set dual variable to zero.
        dualvar[bl] = 0
        # Relabel vertices.
        for v_v in blossom_leaves(bl):
            if label[inblossom[v_v]] == 2:
                # This T-vertex now turns into an S-vertex because it becomes
                # part of an S-blossom; add it to the queue.
                queue.append(v_v)
            inblossom[v_v] = bl
        # Compute blossombestedges[b].
        bestedgeto = (2 * n_vertex) * [-1]
        for bv in path:
            if blossombestedges[bv] is None:
                # This subblossom does not have a list of least-slack edges;
                # get the information from the vertices.
                nblists = [[v_p // 2 for v_p in neighbend[v_v]]
                           for v_v in blossom_leaves(bv)]
            else:
                # Walk this subblossom's least-slack edges.
                nblists = [blossombestedges[bv]]
            for nblist in nblists:
                for v_k in nblist:
                    (v_i, v_j, v_wt) = edges[v_k]
                    if inblossom[v_j] == bl:
                        v_i, v_j = v_j, v_i
                    bj = inblossom[v_j]
                    if (bj != bl and label[bj] == 1 and
                            (bestedgeto[bj] == -1 or
                                     slack(v_k) < slack(bestedgeto[bj]))):
                        bestedgeto[bj] = v_k
            # Forget about least-slack edges of the subblossom.
            blossombestedges[bv] = None
            bestedge[bv] = -1
        blossombestedges[bl] = [v_k for v_k in bestedgeto if v_k != -1]
        # Select bestedge[b].
        bestedge[bl] = -1
        for v_k in blossombestedges[bl]:
            if bestedge[bl] == -1 or slack(v_k) < slack(bestedge[bl]):
                bestedge[bl] = v_k

    # Expand the given top-level blossom.
    def expand_blossom(int bl, end_stage):
        cdef int s, v_v, v_p, v_j
        # Convert sub-blossoms into top-level blossoms.
        for s in blossomchilds[bl]:
            blossomparent[s] = -1
            if s < n_vertex:
                inblossom[s] = s
            elif end_stage and dualvar[s] == 0:
                # Recursively expand this sub-blossom.
                expand_blossom(s, end_stage)
            else:
                for v_v in blossom_leaves(s):
                    inblossom[v_v] = s
        # If we expand a T-blossom during a stage, its sub-blossoms must be
        # relabeled.
        if (not end_stage) and label[bl] == 2:
            # Start at the sub-blossom through which the expanding
            # blossom obtained its label, and relabel sub-blossoms untili
            # we reach the base.
            # Figure out through which sub-blossom the expanding blossom
            # obtained its label initially.
            assert labelend[bl] >= 0
            entry_child = inblossom[endpoint[labelend[bl] ^ 1]]
            # Decide in which direction we will go round the blossom.
            v_j = blossomchilds[bl].index(entry_child)
            if v_j & 1:
                # Start index is odd; go forward and wrap.
                v_j -= len(blossomchilds[bl])
                j_step = 1
                endptrick = 0
            else:
                # Start index is even; go backward.
                j_step = -1
                endptrick = 1
            # Move along the blossom until we get to the base.
            v_p = labelend[bl]
            while v_j != 0:
                # Relabel the T-sub-blossom.
                label[endpoint[v_p ^ 1]] = 0
                label[endpoint[blossomendps[bl][v_j - endptrick] ^ endptrick ^ 1]] = 0
                assign_label(endpoint[v_p ^ 1], 2, v_p)
                # Step to the next S-sub-blossom and note its forward endpoint.
                allowedge[blossomendps[bl][v_j - endptrick] // 2] = True
                v_j += j_step
                v_p = blossomendps[bl][v_j - endptrick] ^ endptrick
                # Step to the next T-sub-blossom.
                allowedge[v_p // 2] = True
                v_j += j_step
            # Relabel the base T-sub-blossom WITHOUT stepping through to
            # its mate (so don't call assignLabel).
            bv = blossomchilds[bl][v_j]
            label[endpoint[v_p ^ 1]] = label[bv] = 2
            labelend[endpoint[v_p ^ 1]] = labelend[bv] = v_p
            bestedge[bv] = -1
            # Continue along the blossom until we get back to entrychild.
            v_j += j_step
            while blossomchilds[bl][v_j] != entry_child:
                # Examine the vertices of the sub-blossom to see whether
                # it is reachable from a neighbouring S-vertex outside the
                # expanding blossom.
                bv = blossomchilds[bl][v_j]
                if label[bv] == 1:
                    # This sub-blossom just got label S through one of its
                    # neighbours; leave it.
                    v_j += j_step
                    continue
                for v_v in blossom_leaves(bv):
                    if label[v_v] != 0:
                        break
                # If the sub-blossom contains a reachable vertex, assign
                # label T to the sub-blossom.
                if label[v_v] != 0:
                    assert label[v_v] == 2
                    assert inblossom[v_v] == bv
                    label[v_v] = 0
                    label[endpoint[mate[blossombase[bv]]]] = 0
                    assign_label(v_v, 2, labelend[v_v])
                v_j += j_step
        # Recycle the blossom number.
        label[bl] = labelend[bl] = -1
        blossomchilds[bl] = blossomendps[bl] = None
        blossombase[bl] = -1
        blossombestedges[bl] = None
        bestedge[bl] = -1
        unusedblossoms.append(bl)

    # Swap matched/unmatched edges over an alternating path through blossom b
    # between vertex v and the base vertex. Keep blossom bookkeeping consistent.
    def augment_blossom(int v_b, int v_v):
        cdef int v_t, v_i, v_j, jstep, endptrick, p
        # Bubble up through the blossom tree from vertex v to an immediate
        # sub-blossom of b.
        v_t = v_v
        while blossomparent[v_t] != v_b:
            v_t = blossomparent[v_t]
        # Recursively deal with the first sub-blossom.
        if v_t >= n_vertex:
            augment_blossom(v_t, v_v)
        # Decide in which direction we will go round the blossom.
        v_i = v_j = blossomchilds[v_b].index(v_t)
        if v_i & 1:
            # Start index is odd; go forward and wrap.
            v_j -= len(blossomchilds[v_b])
            jstep = 1
            endptrick = 0
        else:
            # Start index is even; go backward.
            jstep = -1
            endptrick = 1
        # Move along the blossom until we get to the base.
        while v_j != 0:
            # Step to the next sub-blossom and augment it recursively.
            v_j += jstep
            v_t = blossomchilds[v_b][v_j]
            p = blossomendps[v_b][v_j - endptrick] ^ endptrick
            if v_t >= n_vertex:
                augment_blossom(v_t, endpoint[p])
            # Step to the next sub-blossom and augment it recursively.
            v_j += jstep
            v_t = blossomchilds[v_b][v_j]
            if v_t >= n_vertex:
                augment_blossom(v_t, endpoint[p ^ 1])
            # Match the edge connecting those sub-blossoms.
            mate[endpoint[p]] = p ^ 1
            mate[endpoint[p ^ 1]] = p
        # Rotate the list of sub-blossoms to put the new base at the front.
        blossomchilds[v_b] = blossomchilds[v_b][v_i:] + blossomchilds[v_b][:v_i]
        blossomendps[v_b] = blossomendps[v_b][v_i:] + blossomendps[v_b][:v_i]
        blossombase[v_b] = blossombase[blossomchilds[v_b][0]]
        assert blossombase[v_b] == v_v

    # Swap matched/unmatched edges over an alternating path between two
    # single vertices. The augmenting path runs through edge k, which
    # connects a pair of S vertices.
    def augmentMatching(int v_k):
        cdef int v_v,v_w,v_wt, s, v_p, bs
        v_v, v_w, v_wt = edges[v_k]
        for (s, v_p) in ((v_v, 2 * v_k + 1), (v_w, 2 * v_k)):
            # Match vertex s to remote endpoint p. Then trace back from s
            # until we find a single vertex, swapping matched and unmatched
            # edges as we go.
            while 1:
                bs = inblossom[s]
                assert label[bs] == 1
                assert labelend[bs] == mate[blossombase[bs]]
                # Augment through the S-blossom from s to base.
                if bs >= n_vertex:
                    augment_blossom(bs, s)
                # Update mate[s]
                mate[s] = v_p
                # Trace one step back.
                if labelend[bs] == -1:
                    # Reached single vertex; stop.
                    break
                v_t = endpoint[labelend[bs]]
                bt = inblossom[v_t]
                assert label[bt] == 2
                # Trace one step back.
                assert labelend[bt] >= 0
                s = endpoint[labelend[bt]]
                v_j = endpoint[labelend[bt] ^ 1]
                # Augment through the T-blossom from j to base.
                assert blossombase[bt] == v_t
                if bt >= n_vertex:
                    augment_blossom(bt, v_j)
                # Update mate[j]
                mate[v_j] = labelend[bt]
                # Keep the opposite endpoint;
                # it will be assigned to mate[s] in the next step.
                v_p = labelend[bt] ^ 1


    # Main loop: continue until no further improvement is possible.
    for t in range(n_vertex):

        # Each iteration of this loop is a "stage".
        # A stage finds an augmenting path and uses that to improve
        # the matching.

        # Remove labels from top-level blossoms/vertices.
        label[:] = (2 * n_vertex) * [0]

        # Forget all about least-slack edges.
        bestedge[:] = (2 * n_vertex) * [-1]
        blossombestedges[n_vertex:] = n_vertex * [None]

        # Loss of labeling means that we can not be sure that currently
        # allowable edges remain allowable througout this stage.
        allowedge[:] = nedge * [False]

        # Make queue empty.
        queue[:] = []

        # Label single blossoms/vertices with S and put them in the queue.
        for v in range(n_vertex):
            if mate[v] == -1 and label[inblossom[v]] == 0:
                assign_label(v, 1, -1)

        # Loop until we succeed in augmenting the matching.
        augmented = 0
        while 1:

            # Each iteration of this loop is a "substage".
            # A substage tries to find an augmenting path;
            # if found, the path is used to improve the matching and
            # the stage ends. If there is no augmenting path, the
            # primal-dual method is used to pump some slack out of
            # the dual variables.

            # Continue labeling until all vertices which are reachable
            # through an alternating path have got a label.
            while queue and not augmented:

                # Take an S vertex from the queue.
                v = queue.pop()
                assert label[inblossom[v]] == 1

                # Scan its neighbours:
                for p in neighbend[v]:
                    k = p // 2
                    w = endpoint[p]
                    # w is a neighbour to v
                    if inblossom[v] == inblossom[w]:
                        # this edge is internal to a blossom; ignore it
                        continue
                    if not allowedge[k]:
                        kslack = slack(k)
                        if kslack <= 0:
                            # edge k has zero slack => it is allowable
                            allowedge[k] = True
                    if allowedge[k]:
                        if label[inblossom[w]] == 0:
                            # (C1) w is a free vertex;
                            # label w with T and label its mate with S (R12).
                            assign_label(w, 2, p ^ 1)
                        elif label[inblossom[w]] == 1:
                            # (C2) w is an S-vertex (not in the same blossom);
                            # follow back-links to discover either an
                            # augmenting path or a new blossom.
                            base = scan_blossom(v, w)
                            if base >= 0:
                                # Found a new blossom; add it to the blossom
                                # bookkeeping and turn it into an S-blossom.
                                add_blossom(base, k)
                            else:
                                # Found an augmenting path; augment the
                                # matching and end this stage.
                                augmentMatching(k)
                                augmented = 1
                                break
                        elif label[w] == 0:
                            # w is inside a T-blossom, but w itself has not
                            # yet been reached from outside the blossom;
                            # mark it as reached (we need this to relabel
                            # during T-blossom expansion).
                            assert label[inblossom[w]] == 2
                            label[w] = 2
                            labelend[w] = p ^ 1
                    elif label[inblossom[w]] == 1:
                        # keep track of the least-slack non-allowable edge to
                        # a different S-blossom.
                        b = inblossom[v]
                        if bestedge[b] == -1 or kslack < slack(bestedge[b]):
                            bestedge[b] = k
                    elif label[w] == 0:
                        # w is a free vertex (or an unreached vertex inside
                        # a T-blossom) but we can not reach it yet;
                        # keep track of the least-slack edge that reaches w.
                        if bestedge[w] == -1 or kslack < slack(bestedge[w]):
                            bestedge[w] = k

            if augmented:
                break

            # There is no augmenting path under these constraints;
            # compute delta and reduce slack in the optimization problem.
            # (Note that our vertex dual variables, edge slacks and delta's
            # are pre-multiplied by two.)
            deltatype = -1
            delta = deltaedge = deltablossom = None

            # Compute delta1: the minumum value of any vertex dual.
            if not maxcardinality:
                deltatype = 1
                delta = min(dualvar[:n_vertex])

            # Compute delta2: the minimum slack on any edge between
            # an S-vertex and a free vertex.
            for v in range(n_vertex):
                if label[inblossom[v]] == 0 and bestedge[v] != -1:
                    d = slack(bestedge[v])
                    if deltatype == -1 or d < delta:
                        delta = d
                        deltatype = 2
                        deltaedge = bestedge[v]

            # Compute delta3: half the minimum slack on any edge between
            # a pair of S-blossoms.
            for b in range(2 * n_vertex):
                if (blossomparent[b] == -1 and label[b] == 1 and
                            bestedge[b] != -1):
                    kslack = slack(bestedge[b])
                    if isinstance(kslack, integer_types):
                        assert (kslack % 2) == 0
                        d = kslack // 2
                    else:
                        d = kslack / 2
                    if deltatype == -1 or d < delta:
                        delta = d
                        deltatype = 3
                        deltaedge = bestedge[b]

            # Compute delta4: minimum z variable of any T-blossom.
            for b in range(n_vertex, 2 * n_vertex):
                if (blossombase[b] >= 0 and blossomparent[b] == -1 and
                            label[b] == 2 and
                        (deltatype == -1 or dualvar[b] < delta)):
                    delta = dualvar[b]
                    deltatype = 4
                    deltablossom = b

            if deltatype == -1:
                # No further improvement possible; max-cardinality optimum
                # reached. Do a final delta update to make the optimum
                # verifyable.
                assert maxcardinality
                deltatype = 1
                delta = max(0, min(dualvar[:n_vertex]))

            # Update dual variables according to delta.
            for v in range(n_vertex):
                if label[inblossom[v]] == 1:
                    # S-vertex: 2*u = 2*u - 2*delta
                    dualvar[v] -= delta
                elif label[inblossom[v]] == 2:
                    # T-vertex: 2*u = 2*u + 2*delta
                    dualvar[v] += delta
            for b in range(n_vertex, 2 * n_vertex):
                if blossombase[b] >= 0 and blossomparent[b] == -1:
                    if label[b] == 1:
                        # top-level S-blossom: z = z + 2*delta
                        dualvar[b] += delta
                    elif label[b] == 2:
                        # top-level T-blossom: z = z - 2*delta
                        dualvar[b] -= delta

            # Take action at the point where minimum delta occurred.
            if deltatype == 1:
                # No further improvement possible; optimum reached.
                break
            elif deltatype == 2:
                # Use the least-slack edge to continue the search.
                allowedge[deltaedge] = True
                (i, j, wt) = edges[deltaedge]
                if label[inblossom[i]] == 0:
                    i, j = j, i
                assert label[inblossom[i]] == 1
                queue.append(i)
            elif deltatype == 3:
                # Use the least-slack edge to continue the search.
                allowedge[deltaedge] = True
                (i, j, wt) = edges[deltaedge]
                assert label[inblossom[i]] == 1
                queue.append(i)
            elif deltatype == 4:
                # Expand the least-z blossom.
                expand_blossom(deltablossom, False)

                # End of a this substage.

        # Stop when no more augmenting path can be found.
        if not augmented:
            break

        # End of a stage; expand all S-blossoms which have dualvar = 0.
        for b in range(n_vertex, 2 * n_vertex):
            if (blossomparent[b] == -1 and blossombase[b] >= 0 and
                        label[b] == 1 and dualvar[b] == 0):
                expand_blossom(b, True)

    # Transform mate[] such that mate[v] is the vertex to which v is paired.
    for v in range(n_vertex):
        if mate[v] >= 0:
            mate[v] = endpoint[mate[v]]

    for a,b in enumerate(mate):
      if a<b:
        yield a,b


# Unit tests
if __name__ == '__main__':
    import unittest, math


    class MaxWeightMatchingTests(unittest.TestCase):
        def test10_empty(self):
            # empty input graph
            self.assertEqual(maxWeightMatching([]), [])

        def test11_singleedge(self):
            # single edge
            self.assertEqual(maxWeightMatching([(0, 1, 1)]), [1, 0])

        def test12(self):
            self.assertEqual(maxWeightMatching([(1, 2, 10), (2, 3, 11)]), [-1, -1, 3, 2])

        def test13(self):
            self.assertEqual(maxWeightMatching([(1, 2, 5), (2, 3, 11), (3, 4, 5)]), [-1, -1, 3, 2, -1])

        def test14_maxcard(self):
            # maximum cardinality
            self.assertEqual(maxWeightMatching([(1, 2, 5), (2, 3, 11), (3, 4, 5)], True), [-1, 2, 1, 4, 3])

        def test15_float(self):
            # floating point weigths
            self.assertEqual(
                maxWeightMatching([(1, 2, math.pi), (2, 3, math.exp(1)), (1, 3, 3.0), (1, 4, math.sqrt(2.0))]),
                [-1, 4, 3, 2, 1])

        def test16_negative(self):
            # negative weights
            self.assertEqual(maxWeightMatching([(1, 2, 2), (1, 3, -2), (2, 3, 1), (2, 4, -1), (3, 4, -6)], False),
                             [-1, 2, 1, -1, -1])
            self.assertEqual(maxWeightMatching([(1, 2, 2), (1, 3, -2), (2, 3, 1), (2, 4, -1), (3, 4, -6)], True),
                             [-1, 3, 4, 1, 2])

        def test20_sblossom(self):
            # create S-blossom and use it for augmentation
            self.assertEqual(maxWeightMatching([(1, 2, 8), (1, 3, 9), (2, 3, 10), (3, 4, 7)]), [-1, 2, 1, 4, 3])
            self.assertEqual(maxWeightMatching([(1, 2, 8), (1, 3, 9), (2, 3, 10), (3, 4, 7), (1, 6, 5), (4, 5, 6)]),
                             [-1, 6, 3, 2, 5, 4, 1])

        def test21_tblossom(self):
            # create S-blossom, relabel as T-blossom, use for augmentation
            self.assertEqual(maxWeightMatching([(1, 2, 9), (1, 3, 8), (2, 3, 10), (1, 4, 5), (4, 5, 4), (1, 6, 3)]),
                             [-1, 6, 3, 2, 5, 4, 1])
            self.assertEqual(maxWeightMatching([(1, 2, 9), (1, 3, 8), (2, 3, 10), (1, 4, 5), (4, 5, 3), (1, 6, 4)]),
                             [-1, 6, 3, 2, 5, 4, 1])
            self.assertEqual(maxWeightMatching([(1, 2, 9), (1, 3, 8), (2, 3, 10), (1, 4, 5), (4, 5, 3), (3, 6, 4)]),
                             [-1, 2, 1, 6, 5, 4, 3])

        def test22_s_nest(self):
            # create nested S-blossom, use for augmentation
            self.assertEqual(
                maxWeightMatching([(1, 2, 9), (1, 3, 9), (2, 3, 10), (2, 4, 8), (3, 5, 8), (4, 5, 10), (5, 6, 6)]),
                [-1, 3, 4, 1, 2, 6, 5])

        def test23_s_relabel_nest(self):
            # create S-blossom, relabel as S, include in nested S-blossom
            self.assertEqual(maxWeightMatching(
                [(1, 2, 10), (1, 7, 10), (2, 3, 12), (3, 4, 20), (3, 5, 20), (4, 5, 25), (5, 6, 10), (6, 7, 10),
                 (7, 8, 8)]), [-1, 2, 1, 4, 3, 6, 5, 8, 7])

        def test24_s_nest_expand(self):
            # create nested S-blossom, augment, expand recursively
            self.assertEqual(maxWeightMatching(
                [(1, 2, 8), (1, 3, 8), (2, 3, 10), (2, 4, 12), (3, 5, 12), (4, 5, 14), (4, 6, 12), (5, 7, 12),
                 (6, 7, 14), (7, 8, 12)]), [-1, 2, 1, 5, 6, 3, 4, 8, 7])

        def test25_s_t_expand(self):
            # create S-blossom, relabel as T, expand
            self.assertEqual(maxWeightMatching(
                [(1, 2, 23), (1, 5, 22), (1, 6, 15), (2, 3, 25), (3, 4, 22), (4, 5, 25), (4, 8, 14), (5, 7, 13)]),
                             [-1, 6, 3, 2, 8, 7, 1, 5, 4])

        def test26_s_nest_t_expand(self):
            # create nested S-blossom, relabel as T, expand
            self.assertEqual(maxWeightMatching(
                [(1, 2, 19), (1, 3, 20), (1, 8, 8), (2, 3, 25), (2, 4, 18), (3, 5, 18), (4, 5, 13), (4, 7, 7),
                 (5, 6, 7)]), [-1, 8, 3, 2, 7, 6, 5, 4, 1])

        def test30_tnasty_expand(self):
            # create blossom, relabel as T in more than one way, expand, augment
            self.assertEqual(maxWeightMatching(
                [(1, 2, 45), (1, 5, 45), (2, 3, 50), (3, 4, 45), (4, 5, 50), (1, 6, 30), (3, 9, 35), (4, 8, 35),
                 (5, 7, 26), (9, 10, 5)]), [-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9])

        def test31_tnasty2_expand(self):
            # again but slightly different
            self.assertEqual(maxWeightMatching(
                [(1, 2, 45), (1, 5, 45), (2, 3, 50), (3, 4, 45), (4, 5, 50), (1, 6, 30), (3, 9, 35), (4, 8, 26),
                 (5, 7, 40), (9, 10, 5)]), [-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9])

        def test32_t_expand_leastslack(self):
            # create blossom, relabel as T, expand such that a new least-slack S-to-free edge is produced, augment
            self.assertEqual(maxWeightMatching(
                [(1, 2, 45), (1, 5, 45), (2, 3, 50), (3, 4, 45), (4, 5, 50), (1, 6, 30), (3, 9, 35), (4, 8, 28),
                 (5, 7, 26), (9, 10, 5)]), [-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9])

        def test33_nest_tnasty_expand(self):
            # create nested blossom, relabel as T in more than one way, expand outer blossom such that inner blossom ends up on an augmenting path
            self.assertEqual(maxWeightMatching(
                [(1, 2, 45), (1, 7, 45), (2, 3, 50), (3, 4, 45), (4, 5, 95), (4, 6, 94), (5, 6, 94), (6, 7, 50),
                 (1, 8, 30), (3, 11, 35), (5, 9, 36), (7, 10, 26), (11, 12, 5)]),
                             [-1, 8, 3, 2, 6, 9, 4, 10, 1, 5, 7, 12, 11])

        def test34_nest_relabel_expand(self):
            # create nested S-blossom, relabel as S, expand recursively
            self.assertEqual(maxWeightMatching(
                [(1, 2, 40), (1, 3, 40), (2, 3, 60), (2, 4, 55), (3, 5, 55), (4, 5, 50), (1, 8, 15), (5, 7, 30),
                 (7, 6, 10), (8, 10, 10), (4, 9, 30)]), [-1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8])


    CHECK_DELTA = True
    unittest.main()

# end
