#!/usr/bin/env python2
import argparse

import pyximport;

pyximport.install()

import matplotlib
matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt

import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a plot from tabular file(s).")
    parser.add_argument("x", type=str, help="Column for x")
    parser.add_argument("y", type=str, help="Column for y")
    parser.add_argument("files", type=str, nargs="+", help="Tabular files.")
    parser.add_argument("-xlim", type=float, nargs=2)
    parser.add_argument("-ylim", type=float, nargs=2)
    param = parser.parse_args()

    cl = plt.cm.Paired([x - (1.0 / 24) for x in np.linspace(0, 1, 12)])[1::2]
    alpha = 1.0
    for idx, tabfile in enumerate(param.files):
        with open(tabfile) as f:
            l = f.readline()
            header = l.strip().split()
            data = {tab: [] for tab in header}
            for l in f:
                for tab, val in zip(header, map(float, l.strip().split())):
                    data[tab].append(val)
                data["dcj_distance"][-1] = data["real_distance"][-1] - data["dcj_distance"][-1]
                data["ortho_TP"][-1] /= data["ortho_TOTAL"][-1]
                if data["real_distance"][-1] >= 600:
                    break

        plt.plot(data[param.x], data[param.y], label=tabfile)
        # plt.scatter(data[param.x], data[param.y], label=tabfile, c=cl[idx % 6], alpha=alpha,
        #             edgecolors="none" if idx < 6 else "black")
    plt.legend(fontsize=8, loc='upper left')
    plt.grid(True)
    if param.xlim is not None:
        plt.xlim(param.xlim)
    if param.ylim is not None:
        plt.ylim(param.ylim)
    plt.savefig('tab_%s_vs_%s.pdf' % (param.y, param.x), bbox_inches='tight')
