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

    data = {"DCJ": [], "ML": [], "eDCJ": [], "Avg": [], 'step': []}
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
                data["Avg"].append((est+ml)/2 - step)
                data["step"].append(step)
                # data["dcj"].append(dcj_rep_data)
                # data["ml"].append(ml_rep_data)
                # data["est"].append(est_rep_data)

    # PANDA:
    df = pd.DataFrame.from_dict(data)
    # pd.options.display.mpl_style = 'default'
    # import ipdb; ipdb.set_trace()
    df.boxplot(by='step')
    plt.savefig('rw_%s_results.pdf' % param.n, bbox_inches='tight')

