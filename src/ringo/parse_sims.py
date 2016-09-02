#! /usr/bin/env python2
import ringo_config
cfg = ringo_config.RingoConfig()
import pyximport; pyximport.install(build_dir=cfg.pyximport_build())
import argparse
import re
import datetime
import math
from tabulate import tabulate
import file_ops
from simulation import Simulation
import os
from ringo_config import RingoConfig
import sys
import collections
import algorithms
from matplotlib.font_manager import FontProperties
import matplotlib
import numpy as np
import dcj
matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font', family='serif')
## for Palatino and other serif fonts use:
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


def td_format(td_object):
    seconds = int(td_object.total_seconds())
    periods = [
        ('year', 60 * 60 * 24 * 365),
        ('month', 60 * 60 * 24 * 30),
        ('d', 60 * 60 * 24),
        ('h', 60 * 60),
        ('m', 60),
        ('s', 1)
    ]

    strings = []
    for period_name, period_seconds in periods:
        if seconds > period_seconds:
            period_value, seconds = divmod(seconds, period_seconds)
            if period_name == "s":
                pass
                strings.append("%.0fs" % (period_value + math.modf(seconds)[0]))
            else:
                strings.append("%s%s" % (period_value, period_name))

    return "".join(strings)


def parse_filetime(file):
    with open(file) as f:
        t = 0
        for l in f:
            if l.startswith("user"):
                s, time = l.strip().split()
                m = re.search("(\d+)m(\d+)\.(\d+)s", time)
                t += 60 * int(m.group(1)) + int(m.group(2)) + float(m.group(3)) / 1000
        return t




def get_stats(ancestral_g, reconstructed_g):
    a_set = ancestral_g.adjacency_set()
    r_set = reconstructed_g.adjacency_set()
    return {'tp': len(a_set.intersection(r_set)), 'fp': len(r_set.difference(a_set)),
            'fn': len(a_set.difference(r_set)), 'ancestral': len(a_set), 'reconstructed': len(r_set)}


def get_ancestral_stats(ancestral_genomes, reconstructed_genomes, tree, exclude_root=True,
                        calc_tree_distances=False, calc_dcj_distance=False):
    stats = []
    for node in tree.preorder_internal_node_iter():
        if node == tree.seed_node and exclude_root:
            continue
        r = get_stats(ancestral_genomes[node.label], reconstructed_genomes[node.label])
        # distances:
        if calc_tree_distances:
            r.update({k:node.__dict__[k] for k in ['min_d_leaf','max_d_leaf','avg_d_leaf','d_root','total_ev']})
            # r.update({'min_d_leaf':node.min_d_leaf, 'max_d_leaf':node.max_d_leaf, 'avg_d_leaf':node.avg_d_leaf,
            #           'd_root':node.d_root, 'total_ev':node.total_ev})
        if calc_dcj_distance:
            r.update({'dcj_distance':dcj.dcj_distance(ancestral_genomes[node.label], reconstructed_genomes[node.label])})
        r['label'] = node.label
        stats.append(r)
    return stats

def sum_ds_result(ds, stats=None):
    if stats is None:
        stats = ds[0].keys()
    total = {k:0 for k in stats}
    for sim_result in ds:
        for key in stats:
            total[key] += sim_result[key]
    return total

def average_ds_result(ds, stats=None):
    if stats is None:
        stats = ds[0].keys()
    total = {k:[] for k in stats}
    for sim_result in ds:
        for key in stats:
            total[key].append(sim_result[key])
    return {key:np.mean(v) for key,v in total.iteritems()}

def show_time_tabular(results, ds_per_row=6, fmt="rst"):
    datasets = sorted(results[algs[0]].keys())
    for ds_slice in [datasets[i:i+ds_per_row] for i in range(0,len(datasets),ds_per_row)]:
        table = []
        for alg in algs:
            table.append([alg] + [td_format(datetime.timedelta(seconds=np.mean(results[alg][ds]))) for ds in ds_slice])
        print
        print "DataSets:", ds_slice
        # print tabulate(table, headers=["algorithms"] + ([""]*len(ds_slice)), tablefmt=fmt)
        print tabulate(table, headers=["algorithms"] + ds_slice, tablefmt=fmt)
        print tabulate(table, headers=["algorithms"] + ds_slice, tablefmt="latex_booktabs")


def show_tabular(results, algs, stats, ds_per_row=None, fmt="rst"):
    if ds_per_row is None:
        ds_per_row = 12 / len(stats)
    datasets = sorted(results[algs[0]].keys())
    for ds_slice in [datasets[i:i+ds_per_row] for i in range(0,len(datasets),ds_per_row)]:
        table = []

        for alg in algs:
            row = [alg]
            for ds in ds_slice:
                r = sum_ds_result(results[alg][ds], stats=stats+["ancestral","reconstructed"])
                row.extend([100 * float(r[s])/r["ancestral"] if r["ancestral"]>0 else 0 for s in stats])
            table.append(row)
        print
        print "DataSets:", ds_slice
        print tabulate(table, headers=["algorithms"] + (stats*len(ds_slice)), tablefmt=fmt, floatfmt=".1f")

def output_csv(filename, results, algs, stats):
    datasets = sorted(results[algs[0]].keys())
    with open(filename,"w") as f:
        # header:
        # import ipdb; ipdb.set_trace()
        f.write("Algorithms\t%s\n" % ("\t".join([ds+"\t"*(len(stats)-1)  for ds in datasets])))
        f.write("\t"+"\t".join(["\t".join(stats)] * len(datasets))+"\n")
        table = []
        for alg in algs:
            row = [alg]
            for ds in datasets:
                # r = sum_ds_result(results[alg][ds], stats=stats+['ancestral'])
                # row.extend(["%.5f" % (float(r[s])/r["ancestral"]) if r["ancestral"]>0 else "-" for s in stats])
                for r in results[alg][ds]:
                    row.append("%s" % r["label"])
                    row.extend(["%d" % float(r[s]) for s in stats])


            # table.append(row)
            f.write("\t".join(row)+"\n")
        # f.write(tabulate(table, headers=["algorithms"] + (stats*len(datasets)), tablefmt="plain", floatfmt=".3f"))
    # f.write(alg)
    # f.write("%s\t%s\n" % (alg, "\t".join(["\t".join([str(float(r[s])/r['ancestral']) if r['ancestral']>0 else "-" for s in stats])  for ds in datasets]))_




def plot_scatter(results, x_stats="total_ev", y_stats="tp"):
    cl = plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[1::2]

    ALPHA=0.8
    algs = results.keys()
    datasets = sorted(results[algs[0]].keys())
    for ds in datasets:
        plt.clf()
        plt.cla()
        for color_idx, alg in enumerate(algs):
            x_vector = [r[x_stats] for r in results[alg][ds]]
            y_vector = [float(r[y_stats])/r['ancestral'] for r in results[alg][ds]]
            plt.scatter(x_vector, y_vector, c=cl[color_idx], s=6, alpha=ALPHA, edgecolors='none', label=alg)
        plt.legend(fontsize=6, loc='upper right')
        plt.grid(True)
        plt.savefig('scatter_%s.pdf' % ds, bbox_inches='tight')



def plot_distance(filename, results, algs, title=None,x_axis=None, plot_size=None):
    cl = plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[1::2]


    n_algs = len(algs)
    bar_width = 1.0 / (n_algs + 1)
    datasets = sorted(results[algs[0]].keys())
    n_datasets = len(datasets)
    index = np.arange(n_datasets)

    plt.cla()
    plt.clf()
    if plot_size is not None:
        plt.figure(figsize=tuple(plot_size))
    for idx, alg in enumerate(algs):
        label = "%s" % (alg)
        y = list()
        for ds in datasets:
            r = average_ds_result(results[alg][ds], stats=['dcj_distance'])
            y.append(r['dcj_distance'])


        plt.bar(index + (bar_width * idx), y, bar_width, linewidth=0.1, alpha=1, label=label, color=cl[idx])
    # plt.ylim([0,1.05])
    plt.xlim([-2 * bar_width, n_datasets]) # + bar_width * (n_algs + 1)])
    plt.xticks(index + (bar_width/2*n_algs), datasets, rotation=40, fontsize=9)
    ax = plt.gca()
    plt.legend(fontsize=12, loc='upper left')
    # plt.grid(True)
    # TITLE:
    plt.title(title)
    plt.xlabel(x_axis)

    plt.savefig(filename, bbox_inches='tight')



def plot_bar(filename, results, algs, stats=["tp","fp"], title=None,x_axis=None,plot_size=None):
    # cl = plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[1::2]
    cl = plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[1::2]
    # cl += plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[0::2]
    # print "C:",cl
    alpha={"tp":1, "fp":0.5, "fn":0.3}

    n_algs = len(algs)
    bar_width = 1.0 / (n_algs + 1)
    datasets = sorted(results[algs[0]].keys())
    n_datasets = len(datasets)
    index = np.arange(n_datasets, dtype=float)
    offset = 0.0
    for i in range(n_datasets):
        index[i] += offset
        if (i+1) % 4 == 0:
            offset += bar_width

    plt.cla()
    plt.clf()
    if plot_size is not None:
        plt.figure(figsize=tuple(plot_size))
    for idx, alg in enumerate(algs):
        y = collections.defaultdict(list)
        for ds in datasets:
            r = sum_ds_result(results[alg][ds], stats=stats+['ancestral'])
            for s in stats:
                y[s].append(float(r[s])/r['ancestral'] if r['ancestral']>0 else 0)

        bottom = np.zeros(n_datasets)
        label = "%s (%s)" % (alg, "+".join(map(str.upper,stats)))
        for s in stats:
            plt.bar(index + (bar_width * idx), y[s], bar_width, bottom=bottom,  linewidth=0.1, alpha=alpha[s], label=label, color=cl[idx])
            bottom += y[s]
            label = ""
        # plt.bar(index + (bar_width * idx), y['fp'], bar_width, bottom=y['tp'],  linewidth=0.1, alpha=0.4, label=FP, color="gray")
        # plt.bar(index + (bar_width * idx), y['fn'], bar_width, bottom=[a+b for a,b in zip(y['tp'],y['fp'])],  linewidth=0.1, alpha=0.4, label=FN, color="darkgray")
        # plt.bar(index + (bar_width * idx), y['tp'], bar_width, alpha=ALPHA, linewidth=0.1, color=cl[idx], label=alg)
        # FP=None
        # FN=None

    plt.ylim([0,1.05])
    plt.xlim([-2 * bar_width, max(index)+1+bar_width]) # + bar_width * (n_algs + 1)])
    plt.xticks(index + (bar_width/2*n_algs), datasets, rotation=40, fontsize=9)
    # plt.xticks(index + (bar_width/2*n_algs), datasets, fontsize=9)
    ax = plt.gca()
    # reorder legend:
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[1:] + handles[:1], labels[1:] + labels[:1], loc='lower right', fontsize=8)
    # normal legend:
    plt.legend(fontsize=12, loc='lower right')
    # plt.grid(True)
    # TITLE:
    plt.title(title)
    plt.xlabel(x_axis)

    plt.savefig(filename, bbox_inches='tight')


def parse_sims(folders, algorithm_list, parameters, sim_group, calc_dcj_distance=False, parse_time=False):

    # if algorithm_list is None:
    #     # TODO: auto-detect?
    #     pass
    all_results = collections.defaultdict(lambda: collections.defaultdict(list))
    time_results = collections.defaultdict(lambda: collections.defaultdict(list))
    for idx,folder in enumerate(folders):
        sim = Simulation.open_folder(folder)
        # set tree distances:
        algorithms.set_all_tree_distances(sim.sim_tree)
        # define simulation parameter label:
        sim.sim_parameters.folder = folder
        sim.sim_parameters.idx = idx
        # indel perc
        sim.sim_parameters.indel_p = sim.sim_parameters.insertion_p + sim.sim_parameters.deletion_p
        sim_label = sim_group % sim.sim_parameters.__dict__

        try:
            for alg in algorithm_list:

                label, method, location = alg.split(",")
                # Get ancestral genomes:
                reconstructed = file_ops.load_ancestral_genomes(folder, method, location)

                algo_res = get_ancestral_stats(sim.ancestral_genomes, reconstructed, sim.sim_tree,
                                               calc_tree_distances=True, calc_dcj_distance=calc_dcj_distance)
                if parse_time:
                    # convention: timefile is method.time, for instance, "ringo.time" or "mgra.time"
                    time_results[label][sim_label].append(parse_filetime(os.path.join(folder,location,"%s.time" % method.lower())))
                all_results[label][sim_label].extend(algo_res)
        except (RuntimeError,KeyError,IOError):
            print >> sys.stderr, "Results not present for all methods on folder %s, skipping..." % folder

    return all_results, time_results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parses a list of simulation folders.")
    parser.add_argument("-b", "--bar_plot", type=str,
                        help="PDF filename for a bar plot with adjacency results.")
    parser.add_argument("--bar_plot_title", type=str, default="", help="Optionally give the plot title.")

    parser.add_argument("-d", "--distance_plot", type=str,
                        help="PDF filename for a bar plot with DCJ distance results.")
    parser.add_argument("--distance_plot_title", type=str, default="", help="Optionally give the distance plot title.")

    parser.add_argument("--plot_xaxis", type=str, default="", help="Optionally give the X axis legend.")
    parser.add_argument("--plot_size", type=float, nargs=2, help="Optionally give the plot size, (height and width).")

    parser.add_argument("-a", "--algorithms", type=str, nargs='+',
                        help="Algorithms to compare.")
    parser.add_argument("-t", "--table", nargs='?', const=4, type=int, help="Prints a table(s) with the results. You can optionally give an integer, meaning how many datasets are shown per row of the table output")
    parser.add_argument("-f", "--folders", required=True, type=str, nargs='+',
                        help="Simulation folder(s).")

    parser.add_argument("--parse_time", action="store_true", default=False, help="Outputs a table with running times.")
    parser.add_argument("--latex", action="store_true", default=False, help="Tables in latex format.")
    param = parser.parse_args()


    algorithm_list = param.algorithms
    folders = sorted(param.folders)
    algs = [t.split(",")[0] for t in algorithm_list]
    parameters = ["tp", "fn", "fp"]
    sim_group = "(%(scale).1f,%(indel_p).1f)"
    # sim_group = "%(scale).1f"
    # sim_group = "All"
    # sim_group = "%(folder)s"
    # sim_group = "%(idx)s"
    # sim_group = "(%(scale).1f,%(indel_perc).1f,%(idx)d)"

    calc_dcj_distance = param.distance_plot is not None

    all_results, time_results = parse_sims(folders, algorithm_list, parameters, sim_group,
                             calc_dcj_distance=calc_dcj_distance, parse_time=param.parse_time)

    if param.parse_time:
        show_time_tabular(time_results)

    # print all_results
    output_csv("output.csv", all_results, algs, ["tp","fp"])
    # sys.exit()
    if param.table is not None:
        fmt = "latex_booktabs" if param.latex else "rst"
        show_tabular(all_results, algs, stats=['tp', 'fp'], ds_per_row=param.table, fmt=fmt)

    # plot_scatter(all_results)
    if param.bar_plot is not None:
        plot_bar(param.bar_plot, all_results, algs, title=param.bar_plot_title, x_axis=param.plot_xaxis, plot_size=param.plot_size)

    if param.distance_plot is not None:
        plot_distance(param.distance_plot, all_results, algs, title=param.distance_plot_title, x_axis=param.plot_xaxis, plot_size=param.plot_size)
