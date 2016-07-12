import tempfile
import os
from ringo_config import RingoConfig
import subprocess
import sys

__author__ = 'pfeijao'


def run_blossom5(v, e, edges):
    # config:
    cfg = RingoConfig()
    blossom5 = cfg.blossom5_path()
    # create temporary files for input/output:
    tmpfile_in = tempfile.NamedTemporaryFile(delete=False)
    tmpfile_out = tempfile.NamedTemporaryFile(delete=False)
    solution = []
    try:
        # write size of problem, and then edges:
        tmpfile_in.write("%d %d\n" % (v, e))
        for e in edges:
            tmpfile_in.write("%d %d %f\n" % tuple(e))
        tmpfile_in.close()
        FNULL = open(os.devnull, 'w')
        subprocess.Popen([blossom5, "-e", tmpfile_in.name, "-w", tmpfile_out.name, "-V"], stdout=FNULL).wait()
        # open solution:
        # skip first line
        v, e = tmpfile_out.readline().strip().split()
        for l in tmpfile_out:
            solution.append(map(int,l.strip().split()))
    except ValueError: # empty solution if blossom5 crashes
        print >> sys.stderr, "Blossom5 didn't return a solution, file:", tmpfile_in.name
    finally:
        os.remove(tmpfile_in.name)
        os.remove(tmpfile_out.name)
    return solution
