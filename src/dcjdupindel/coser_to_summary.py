#!/usr/bin/env python2
import argparse
import re

# Best objective 4.937500000000e+02, best bound 4.937500000000e+02, gap 0.0%
sol_regexp = re.compile('Best objective (.+), best bound (.+), gap (.+)%')

# genome2 total      25 duplicons of length   1
dup_regexp = re.compile('genome(\d) total\s+(\d+) duplicons')


def parse_coser_sol(filename):
    duplications = {1: 0, 2: 0}
    obj = gap = 0
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
    logfile = filename.replace("out", "log")
    with open(logfile) as f:
        time = f.readline().strip()

    return obj, duplications[1], duplications[2], time, gap


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a summary file from a list of COSER output files.")
    parser.add_argument("files", type=str, nargs="+", help="Output file(s).")
    param = parser.parse_args()
    print "\t".join(["File", "Dist", "Dup1", "Dup2", "RunningTime", "Gap"])
    for sol_file in param.files:
        r = parse_coser_sol(sol_file)
        title = sol_file.replace(".out", "").replace("coser.", "").replace("G_","").replace(".","_")
        print "\t".join([title] + [str(x) for x in r])
