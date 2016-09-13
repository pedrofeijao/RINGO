#!/usr/bin/env python2
import argparse
import collections
import os
import re

# Best objective 4.937500000000e+02, best bound 4.937500000000e+02, gap 0.0%
sol_regexp = re.compile('Best objective (.+), best bound (.+), gap (.+)%')

# genome2 total      25 duplicons of length   1
dup_regexp = re.compile('genome(\d) total\s+(\d+) duplicons')


def parse_filetime(file):
    with open(file) as f:
        t = 0
        for l in f:
            if l.startswith("real"):
                s, time = l.strip().split()
                m = re.search("(\d+)m(\d+)\.(\d+)s", time)
                t += 60 * int(m.group(1)) + int(m.group(2)) + float(m.group(3)) / 1000
        return t


def parse_coser_sol(folder, correct_matching):
    duplications = {1: 0, 2: 0}
    obj = gap = 0
    filename = os.path.join(folder, "coser.out")
    with open(filename) as f:
        for l in f:
            l = l.strip()
            m = sol_regexp.match(l)
            if m is not None:
                obj, bound, gap = map(float, m.groups())
            m = dup_regexp.match(l)
            if m is not None:
                duplications[int(m.group(1))] = m.group(2)

    # Log for time:
    time = parse_filetime(filename.replace("out", "log"))

    # orthology:
    from parse_orthology import parse_orthology_quality

    # get solution matching:
    # 984_1 984_1
    solution_matching = collections.defaultdict(list)
    with open(os.path.join(folder, "mapping")) as f:
        for l in f:
            pair_a, pair_b = l.strip().split()
            g_a, c_a, g_b, c_b = map(int, pair_a.split("_") + pair_b.split("_"))
            solution_matching[g_a].append((c_a, c_b))

    # compare:
    tp, fp, fn = parse_orthology_quality(solution_matching, correct_matching)
    return {"dcj_distance": obj, "dup_a": duplications[1], "dup_b": duplications[2],
            "time": time, "gap": gap, "ortho_TP": len(tp), "ortho_FP": len(fp), "ortho_FN": len(fn)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of COSER output files.")
    parser.add_argument("files", type=str, nargs="+", help="Output file(s).")
    param = parser.parse_args()
    print "\t".join(["File", "Dist", "Dup1", "Dup2", "RunningTime", "Gap"])
    for sol_file in param.files:
        r = parse_coser_sol(sol_file)
        title = sol_file.replace(".out", "").replace("coser.", "").replace(".", "_")
        print "\t".join([title] + [str(x) for x in r])
