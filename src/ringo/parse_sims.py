#! /usr/bin/env python2
import pyximport
pyximport.install()
import argparse
from tabulate import tabulate
import file_ops
from simulation import Simulation
import os
from ringo_config import RingoConfig
import sys

def get_total_stats(ancestral_genomes, reconstructed_genomes, average=False, exclude_root=True):
    stats = ['tp', 'fp', 'fn', 'ancestral', 'reconstructed']
    count = {st: 0 for st in stats}
    ancestral_adj = 0
    total = 0
    for label, anc in ancestral_genomes.iteritems():
        if label == "Root" and exclude_root:
            continue
        total += 1
        r = get_stats(anc, reconstructed_genomes[label])
        for st in stats:
            count[st] += r[st]
    if average:
        count = {st: float(val)/total for st, val in count.iteritems()}
    # count['tp'] = float(count['tp'])/count['ancestral']
    # count['fn'] = float(count['fn'])/count['ancestral']
    # count['fp'] = float(count['fp'])/count['reconstructed']
    return count


def get_stats(ancestral_g, reconstructed_g):
    a_set = ancestral_g.adjacency_set()
    r_set = reconstructed_g.adjacency_set()
    return {'tp': len(a_set.intersection(r_set)), 'fp': len(r_set.difference(a_set)),
            'fn': len(a_set.difference(r_set)), 'ancestral': len(a_set), 'reconstructed': len(r_set)}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parses a list of simulation folders.")
    parser.add_argument("-o", "--output", type=str,
                        default="adj_reconstruction.pdf",
                        help="Output PDF file for plot.")
    parser.add_argument("-a", "--algorithms", type=str, nargs='+',
                        help="Algorithms to compare.")
    parser.add_argument("-f", "--folders", required=True, type=str, nargs='+',
                        help="Simulation folder(s).")
    parser.add_argument("-c", "--combine", action="store_true", default=False, help="Combine all folders into one result.")
    param = parser.parse_args()

    # res = [["M1",696000,1989100000],["M2",6371,5973.6], ["M3",1737,73.5],["M4",3390,641.85]]

    cfg = RingoConfig()
    algorithms = param.algorithms
    parameters = ["tp", "fn", "fp"]
    if algorithms is None:
        # TODO: auto-detect?
        pass
    all_results = {}
    for folder in param.folders:
        results = []

        sim = Simulation.open_folder(folder)
        try:
            for alg in algorithms:
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

                algo_res = get_total_stats(sim.ancestral_genomes, reconstructed)
                # results.append([label] + ["%.1f %%" % (100*algo_res[par]) for par in parameters])
                results.append([label] + [algo_res[par] for par in parameters])
            all_results[folder] = results
        except:
            print >> sys.stderr, "Results not present for all methods on folder %s, skipping..." % folder

    # print all_results
    if param.combine:
        combined = [[alg.split(",")[0]] + [0]*len(parameters) for alg in algorithms]
        for results in all_results.itervalues():
            for idx, result in enumerate(results):
                for i in xrange(len(parameters)):
                    combined[idx][i+1] += result[i+1]
        print "Combined results:"
        print tabulate(combined, headers=["algorithms"] + parameters, tablefmt="rst", floatfmt=".1f")
    else:
        for folder, results in sorted(all_results.iteritems()):
            print "FOLDER:", folder
            print tabulate(results, headers=["algorithms"] + parameters, tablefmt="rst", floatfmt=".1f")
            print
