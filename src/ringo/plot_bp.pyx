#!/usr/bin/env python2
import math

import sys

import algorithms
import os
import simulation
import file_ops
import scj
from model import BPGraph, CType
from pyx import *
from dendropy import Tree

# Drawing parameters:
COMP_RADIUS = 2  # X_SEP

# styles:
g1_style = [style.linewidth(0.05), color.cmyk.Black]
g2_style = [style.linewidth(0.05), color.cmyk.Gray]
ancestral_style = [style.linewidth(0.2), color.cmyk.Blue, style.linecap.round]
reconstructed_style = [style.linewidth(0.12), color.cmyk.Orange, style.linecap.round]
weight_color = color.cmyk.LimeGreen

def draw_genome(c, adjacency_set, pos, style, draw_only_in_given_components=None):
    for a, b in adjacency_set:
        draw = False
        if draw_only_in_given_components is not None:
            for cycle in draw_only_in_given_components:
                if a in cycle['c'] and b in cycle['c']:
                    draw = True
                    break
        else:
            draw = True
        if draw:
            try:
                c.stroke(path.line(pos[a][0], pos[a][1], pos[b][0], pos[b][1]), style)
            except KeyError:
                # if no key, two cycle that was not drawn
                pass


def draw_bp_graph(c, bp, g1, g2, component_list=None, x_offset=0, y_offset=0, include_2comp=True):
    pos = {}
    if component_list is None:
        component_list = [{'c':val, 'type':key} for key, c_list in bp.type_dict.iteritems() for val in c_list]

    for component in sorted(component_list, key=lambda x: len(x['c']), reverse=True):
        n = len(component['c'])
        if not include_2comp and n <= 2:
            continue
        if component['type'] in [CType.AB_PATH]:
            n += 1
        if component['type'] in [CType.PATH, CType.A_PATH, CType.B_PATH]:
            if n % 2 == 0:
                n += 2
            else:
                n += 1
        # print n
        # if n == 2:
        #     continue
        for x, ext in enumerate(component['c']):
            pos[ext] = (math.cos(2 * math.pi / n * x) * COMP_RADIUS + x_offset,
                         math.sin(2 * math.pi / n * x) * COMP_RADIUS + y_offset)
            # txt = r"%d" % gene
            txt = BPGraph.int2ext(ext)
            c.text(math.cos(2 * math.pi / n * x) * 1.17 * COMP_RADIUS + x_offset - len(str(ext))*0.08,
                   math.sin(2 * math.pi / n * x) * 1.17 * COMP_RADIUS + y_offset - 0.08, txt, [text.size.tiny])
        x_offset += 2.8 * COMP_RADIUS

    # input genomes:
    draw_genome(c, g1.adjacency_set(), pos, g1_style)
    draw_genome(c, g2.adjacency_set(), pos, g2_style)

    # A and B-open vertices:
    a_open = []
    b_open = []
    for component in bp.type_dict[CType.A_PATH]:
        if include_2comp or len(component)>2:
            a_open.append(component[0])
    for component in bp.type_dict[CType.B_PATH]:
        if include_2comp or len(component)>2:
            b_open.append(component[0])
    for component in bp.type_dict[CType.AB_PATH]:
        if include_2comp or len(component)>2:
            a_open.append(component[0])
            b_open.append(component[-1])
    for component in bp.type_dict[CType.AA_PATH]:
        if include_2comp or len(component)>2:
            a_open.append(component[0])
            a_open.append(component[-1])
    for component in bp.type_dict[CType.BB_PATH]:
        if include_2comp or len(component)>2:
            b_open.append(component[0])
            b_open.append(component[-1])

    for px, py in [pos[a] for a in a_open]:
        c.stroke(path.circle(px, py, 0.15), g1_style)
    for px, py in [pos[b] for b in b_open]:
        c.stroke(path.circle(px, py, 0.15), g2_style)

    return pos


def draw_all_bp(reconstructed, tree, folder, internalAdjWeight, ancestral=None):
    # TODO: draw the completion that was used to build the ancestral genome.
    y_offset = 0
    c = canvas.canvas()
    for node in tree.preorder_internal_node_iter():
        if node == tree.seed_node:
            continue
        g1, g2 = reconstructed[node.child_nodes()[0].label], reconstructed[node.child_nodes()[1].label]
        bp = BPGraph(g1, g2)

        c.text(-4*COMP_RADIUS,y_offset+COMP_RADIUS/2.0, "BP(%s,%s). Anc:%s" % (g1.name, g2.name, node.label))
        comp_list = [{'c': val, 'type': key} for key, c_list in bp.type_dict.iteritems() for val in c_list]
        pos = draw_bp_graph(c, bp, g1, g2, y_offset=y_offset)
        y_offset += 2.5 * COMP_RADIUS
        # true ancestor:
        if ancestral is not None:
            draw_genome(c, ancestral[node.label].adjacency_set(), pos, ancestral_style,
                                draw_only_in_given_components=comp_list)
        # reconstructed:
        draw_genome(c, reconstructed[node.label].adjacency_set(), pos,
                            reconstructed_style,
                            draw_only_in_given_components=comp_list)
        # weights:
        for max_w in [x / 10.0 for x in range(1, 11)]:
            draw_genome(c, [adj for adj, w in internalAdjWeight[node.label].iteritems() if w > max_w], pos,
                                [style.linewidth(max_w / 10), weight_color],
                                draw_only_in_given_components=comp_list)

    # legend:
    y_offset -= COMP_RADIUS/2.0
    for name, st in [("Ancestral adjacencies", ancestral_style), ("Adjacency weights",[weight_color]),
            ("Reconstructed Genome", reconstructed_style), ("Genome B", g2_style), ("Genome A", g1_style)]:
        c.stroke(path.line(0, y_offset, COMP_RADIUS, y_offset), st)
        c.text(COMP_RADIUS+0.4, y_offset-0.1, name)
        y_offset += COMP_RADIUS/3.0
    c.text(COMP_RADIUS+0.4, y_offset-0.2, " ")
    #
    #
    # c.stroke(path.line(0, y_offset, COMP_RADIUS, y_offset), ancestral_style)
    # c.text(COMP_RADIUS+0.4, y_offset-0.2, "Ancestral adjacencies")
    #
    # y_offset += COMP_RADIUS/3
    # c.stroke(path.line(0, y_offset, COMP_RADIUS, y_offset), [weight_color])
    # y_offset += COMP_RADIUS/3
    # c.stroke(path.line(0, y_offset, COMP_RADIUS, y_offset), g1_style)
    # y_offset += COMP_RADIUS/3
    # c.stroke(path.line(0, y_offset, COMP_RADIUS, y_offset), g2_style)
    # y_offset += COMP_RADIUS/3

    c.writePDFfile(os.path.join(folder, "bp_plot.pdf"))
