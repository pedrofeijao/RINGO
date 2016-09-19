#!/usr/bin/env python2

import matplotlib

matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt

n = 32
dcj_data = []
ml_data = []
step_range = range(4, 44, 4)
for step in step_range:

    ml_rep_data = []
    dcj_rep_data = []
    for rep in range(1, 11):
        with open("rw.n%d.step%d.rep%d/genomes.txt.ml" % (n, step, rep)) as f:
            f.readline()
            dcj, ml = map(int, f.readline().strip().split())
            dcj_rep_data.append(dcj)
            ml_rep_data.append(ml)
    dcj_data.append(dcj_rep_data)
    ml_data.append(ml_rep_data)
# import ipdb;ipdb.set_trace()
plt.boxplot(dcj_data)
# plt.boxplot(ml_data)
#
plt.plot([x/4 for x in step_range], step_range)
plt.xticks(range(1,len(step_range)+1), step_range)
# plt.yticks(step_range)

plt.savefig('dcj_vs_real.pdf', bbox_inches='tight')
# plt.savefig('ml_vs_real.pdf', bbox_inches='tight')


