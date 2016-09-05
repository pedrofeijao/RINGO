#!/usr/bin/env python2
import argparse
import collections
import os
import pyximport;
pyximport.install()

import matplotlib.pyplot as plt

# TODO: this is a hack to import from other directory; should use packages
ringo_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../ringo"))
os.sys.path.insert(0, ringo_path)
import file_ops

def try_float(c):
    try:
        return float(c)
    except ValueError:
        return c


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parse results from random walk folders.")
    parser.add_argument("folders", type=str, nargs="+", help="folder(s).")
    param = parser.parse_args()

    sim_group = "(%(duplication_length)d,%(duplication_p).2f)"

    results = collections.defaultdict(lambda: collections.defaultdict(list))
    correct = collections.defaultdict(lambda: collections.defaultdict(list))
    for folder in param.folders:
        sim_parameters = file_ops.read_simulation_parameters(folder)
        sim_label = sim_group % sim_parameters
        #results[sim_label].append(sim_parameters)
        distances = {}
        with open(os.path.join(folder, "distances.txt")) as f:
            cols = f.readline().strip().split("\t")
            for l in f:
                tabs = l.strip().split("\t")
                distances[tabs[0]] = {a:b for a, b in zip(cols, map(int, tabs[1:]))}
        with open(os.path.join(folder, "sol_summary.txt")) as f:
            cols = f.readline().strip().split("\t")

            for l in f:
                tabs = l.strip().split("\t")
                genome = "G_%s" % tabs[0].split("_")[1]
                # for col, result in zip(cols, tabs):
                group_results = {a:b for a, b in zip(cols[1:], map(try_float, tabs[1:]))}
                results[sim_label][tabs[0]].append((group_results, distances[genome]))

    # print results
    # Plot distance difference:
    for key, sim in results.iteritems():
        # for each sim group:
        x = []
        y = []
        for pair, group_results in sorted(sim.items(), key=lambda x: x[1][0][1]["total"]):
            x_i = 0
            y_i = 0
            for result in group_results:
                y_i += result[1]["total"] - result[0]["Dist"]
                # y_i += result[1]["duplication"] - result[0]["Indels"]
                x_i += result[1]["total"]
            x.append(x_i/len(group_results))
            y.append(y_i/len(group_results))
        print key, x,y
        plt.plot(x, y, label=key)
    plt.legend(fontsize=12, loc='upper left')
    plt.grid(True)
    # # TITLE:
    # plt.title(title)
    # plt.xlabel(x_axis)
    plt.savefig("parse.pdf", bbox_inches='tight')





