#! /usr/bin/env python2
import ringo_config
cfg = ringo_config.RingoConfig()
import pyximport; pyximport.install(build_dir=cfg.pyximport_build())
import argparse
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

matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font', family='serif')
## for Palatino and other serif fonts use:
# rc('font', **{'family': 'serif', 'serif': ['Palatino']})
# rc('text', usetex=True)


def get_stats(ancestral_g, reconstructed_g):
    a_set = ancestral_g.adjacency_set()
    r_set = reconstructed_g.adjacency_set()
    return {'tp': len(a_set.intersection(r_set)), 'fp': len(r_set.difference(a_set)),
            'fn': len(a_set.difference(r_set)), 'ancestral': len(a_set), 'reconstructed': len(r_set)}


def get_ancestral_stats(ancestral_genomes, reconstructed_genomes, tree, exclude_root=True, calc_distances=False):
    stats = []
    for node in tree.preorder_internal_node_iter():
        if node == tree.seed_node and exclude_root:
            continue
        r = get_stats(ancestral_genomes[node.label], reconstructed_genomes[node.label])
        # distances:
        if calc_distances:
            r.update({k:node.__dict__[k] for k in ['min_d_leaf','max_d_leaf','avg_d_leaf','d_root','total_ev']})
            # r.update({'min_d_leaf':node.min_d_leaf, 'max_d_leaf':node.max_d_leaf, 'avg_d_leaf':node.avg_d_leaf,
            #           'd_root':node.d_root, 'total_ev':node.total_ev})
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

def show_tabular(results, algs, stats=None, ds_per_row=None):
    if ds_per_row is None:
        ds_per_row = 12 / len(stats)
    datasets = sorted(results[algs[0]].keys())
    for ds_slice in [datasets[i:i+ds_per_row] for i in range(0,len(datasets),ds_per_row)]:
        table = []
        for alg in algs:
            row = [alg]
            for ds in ds_slice:
                r = sum_ds_result(results[alg][ds], stats=stats+["ancestral","reconstructed"])
                row.extend([float(r[s])/r["ancestral"] if r["ancestral"]>0 else 0 for s in stats])
            table.append(row)
        print
        print "DataSets:", ds_slice
        print tabulate(table, headers=["algorithms"] + (stats*len(ds_slice)), tablefmt="rst", floatfmt=".3f")

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


def plot_bar(results, algs, stats=["tp","fp"]):
    cl = plt.cm.Paired([x-(1.0/24) for x in np.linspace(0, 1, 12)])[1::2]

    alpha={"tp":1, "fp":0.7, "fn":0.4}

    n_algs = len(algs)
    bar_width = 1.0 / (n_algs + 1)
    FP = "FP (all)"
    FN = "FN (all)"
    datasets = sorted(results[algs[0]].keys())
    n_datasets = len(datasets)
    index = np.arange(n_datasets)

    colors = {}
    # plt.figure(figsize=(12,8))
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
    plt.xlim([-2 * bar_width, n_datasets]) # + bar_width * (n_algs + 1)])
    plt.xticks(index + (bar_width/2*n_algs), datasets, rotation=40, fontsize=9)
    ax = plt.gca()
    # reorder legend:
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[1:] + handles[:1], labels[1:] + labels[:1], loc='lower right', fontsize=8)
    # normal legend:
    plt.legend(fontsize=8, loc='lower right')
    # plt.grid(True)
    plt.savefig('bar_plot.pdf', bbox_inches='tight')


def parse_sims(folders, algorithm_list, parameters, sim_group):

    # if algorithm_list is None:
    #     # TODO: auto-detect?
    #     pass
    all_results = collections.defaultdict(lambda: collections.defaultdict(list))
    for folder in folders:
        sim = Simulation.open_folder(folder)
        # set tree distances:
        algorithms.set_all_tree_distances(sim.sim_tree)
        # define simulation parameter label:
        sim_label = sim_group % sim.sim_parameters.__dict__

        try:
            for alg in algorithm_list:
                label, method, location = alg.split(",")
                # Treat each special case:
                if method.lower() == "mgra":
                    if location == "":
                        location = cfg.mgra_output_folder()
                    reconstructed = file_ops.open_mgra_genomes(
                        os.path.join(folder, location))
                elif method.lower() == "ringo":
                    reconstructed = file_ops.open_genome_file(
                        os.path.join(folder, location, cfg.ringo_output_genomes()))
                elif method.lower() == "physca":
                    reconstructed = file_ops.open_adjacencies_file(
                        os.path.join(folder, location, cfg.physca_reconstructed_adjacencies()))

                else:  # passing in location the path with genome file name should also work
                    reconstructed = file_ops.open_genome_file(
                        os.path.join(folder, location))

                algo_res = get_ancestral_stats(sim.ancestral_genomes, reconstructed, sim.sim_tree, calc_distances=True)
                all_results[label][sim_label].extend(algo_res)
        except (RuntimeError,KeyError,IOError):
            print >> sys.stderr, "Results not present for all methods on folder %s, skipping..." % folder

    return all_results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parses a list of simulation folders.")
    parser.add_argument("-b", "--bar_plot", type=str,
                        help="PDF filename for a bar plot with adjacency results.")
    parser.add_argument("-a", "--algorithms", type=str, nargs='+',
                        help="Algorithms to compare.")
    parser.add_argument("-t", "--table", nargs='?', const=4, type=int, help="Prints a table(s) with the results. You can optionally give an integer, meaning how many datasets are shown per row of the table output")
    parser.add_argument("-f", "--folders", required=True, type=str, nargs='+',
                        help="Simulation folder(s).")
    param = parser.parse_args()


    algorithm_list = param.algorithms
    algs = [t.split(",")[0] for t in algorithm_list]
    parameters = ["tp", "fn", "fp"]
    sim_group = "(%(scale).1f,%(indel_perc).1f)"
    # sim_group = "%(scale).1f"
    # sim_group = "All"

    all_results = parse_sims(param.folders, algorithm_list, parameters, sim_group)

    if param.table is not None:
        show_tabular(all_results, algs, stats=['tp', 'fp'], ds_per_row=param.table)

    # plot_scatter(all_results)
    if param.bar_plot:
        plot_bar(all_results, algs)
