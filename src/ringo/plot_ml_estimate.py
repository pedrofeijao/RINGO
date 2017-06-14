#!/usr/bin/env python2
import argparse

import pandas as pd
import matplotlib

matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plots a boxplot of the DCJ distances and estimates from random walks.")
    parser.add_argument("-n", type=int, required=True, help="Genome size")
    parser.add_argument("-m", type=int, required=True, help="Max size of walk")
    parser.add_argument("-s", type=int, default=4, help="Step in walk lengths")
    parser.add_argument("-r", type=int, default=20, help="Number of repeats")
    param = parser.parse_args()

    data = {"DCJ": [], "ML": [], "eDCJ": [], 'Real Distance': []}
    step_range = range(param.s, param.m, param.s)
    for step in step_range:

        ml_rep_data = []
        dcj_rep_data = []
        est_rep_data = []
        for rep in range(1, 11):
            with open("rw.n%d.step%d.rep%d/genomes.txt.ml" % (param.n, step, rep)) as f:
                l = f.readline()
                dcj, ml, est = map(float, f.readline().strip().split())
                #         dcj_rep_data.append(dcj - step)
                #         ml_rep_data.append(ml - step)
                #         est_rep_data.append(est - step)
                data["DCJ"].append(dcj - step)
                data["ML"].append(ml - step)
                data["eDCJ"].append(est - step)
                data["Real Distance"].append(step)

    # PANDA:
    del data["DCJ"]
    df = pd.DataFrame.from_dict(data)
    df.boxplot(by='Real Distance', layout=(3,1), figsize=(14,10))
    plt.savefig('rw_%s_results.pdf' % param.n, bbox_inches='tight')

