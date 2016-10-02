import math

from model import BPGraph, CType

def dcj_distance(genome1, genome2):
    # find BP graph:
    bp = BPGraph(genome1, genome2)
    # cycles:
    c = len(bp.type_dict[CType.CYCLE]) + len(bp.type_dict[CType.AA_PATH]) + len(bp.type_dict[CType.BB_PATH])
    # odd paths:
    # pure:
    p_o = len([comp for comp in bp.type_dict[CType.PATH] if len(comp) % 2 == 1])
    # from A and B, by merging opposing parity:
    p_a_o = p_b_o = p_b_e = p_a_e = 0
    for comp in bp.type_dict[CType.A_PATH]:
        if len(comp) % 2 == 0:
            p_a_e += 1
        else:
            p_a_o += 1
    for comp in bp.type_dict[CType.B_PATH]:
        if len(comp) % 2 == 0:
            p_b_e += 1
        else:
            p_b_o += 1
    n = len(bp.common_AB) + len(bp.unique_A) + len(bp.unique_B)
    delta = 0
    if (p_a_e < p_a_o and p_b_e < p_b_o) or (p_a_e > p_a_o and p_b_e > p_b_o):
        delta = 1
    return n - c - (p_o + min([p_a_o, p_a_e]) + min([p_b_o, p_b_e]) + delta)/2
