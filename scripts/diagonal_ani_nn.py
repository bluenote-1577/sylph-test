import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from collections import defaultdict
import seaborn as sns
from dataclasses import dataclass
from natsort import natsorted
cmap = sns.color_palette("muted")

np.random.seed(0)
def rand_jitter(arr):
    stdev = .10 
    return arr + np.random.randn(len(arr)) * stdev

@dataclass
class result:
    mean_cov: float
    adj_ani: float
    naive_ani: float
    median_cov: float
    ref_file: str
    query_file: str
    low_ani: float
    high_ani: float
    lam: float
    true_eff_cov: float
    final_ani: float
    low: bool


cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


results = []
true_results = []
prita_files = [
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-ill-c100.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-ill-c1000.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-nano-c100.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-nano-c1000.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-pac-c100.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-pac-c1000.tsv",
        ]

truth_files = [
        "/home/jshaw/projects/prita_test/gtdb-on-reads/true-gtdb-on-on-mock-c100.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/true-gtdb-on-on-mock.tsv",
        ]

mash_files = [
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-ill-mash-s1000.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-nano-mash-s1000.tsv",
        "/home/jshaw/projects/prita_test/gtdb-on-reads/gtdb-on-pac-mash-s1000.tsv",
        ]

mash_results = []
for file in mash_files:
    mash_results.append([])
    for line in open(file,'r'):
        ani = line.split()[0]
        file = line.split()[4].split('/')[-1]
        mash_results[-1].append((file,float(ani) * 100))

for file in prita_files:
    results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split()
        ref_file = spl[1].split('/')[-1]
        query_file = spl[0]
        final_ani = float(spl[2])
        naive_ani = float(spl[3])
        adj_ani = None
        low = False
        if "LOW" in spl[6]:
            low = True
        if spl[5] == "LOW" or spl[5] == "HIGH" or "NA" in spl[4]:
            adj_ani = None
            lam = None
            cis = (None, None)
        else:
            ci = spl[4].split('-')
            adj_ani = float(spl[2])
            lam = float(spl[5])
        if "NA" in spl[4]:
            cis = (None, None)
        else:
            ci = spl[4].split('-')
            cis = float(ci[0]), float(ci[1])
        mean_cov = float(spl[9])
        median_cov = float(spl[8])
        res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam, 0, final_ani, low)
        results[-1].append(res)


for file in truth_files:
    true_results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split()
        ref_file = spl[1].split('/')[-1]
        query_file = spl[0]
        final_ani = float(spl[2])
        naive_ani = float(spl[3])
        adj_ani = None
        low = False
        if "LOW" in spl[6]:
            low = True
        if spl[5] == "LOW" or spl[5] == "HIGH" or "NA" in spl[4]:
            adj_ani = None
            lam = None
            cis = (None, None)
        else:
            ci = spl[4].split('-')
            adj_ani = float(spl[2])
            lam = float(spl[5])
        if "NA" in spl[4]:
            cis = (None, None)
        else:
            ci = spl[4].split('-')
            cis = float(ci[0]), float(ci[1])
        mean_cov = float(spl[9])
        median_cov = float(spl[8])
        res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam, 0, final_ani, low)
        true_results[-1].append(res)

query_to_ani = defaultdict(list)
query_to_err = defaultdict(list)
query_to_anis_nn = defaultdict(list)

for res in true_results[1]:
    if res.final_ani > 90:
        query_to_anis_nn[res.ref_file].append(res.final_ani)
for (key,res) in query_to_anis_nn.items():
    query_to_ani[key].append(np.max(res))
for res in results[3]:
    if not res.low and res.lam != None:
        query_to_ani[res.ref_file].append(res.final_ani)
        query_to_err[res.ref_file].append(res.high_ani)
        query_to_err[res.ref_file].append(res.low_ani)
for res in mash_results[1]:
    query_to_ani[res[0]].append(res[1])

x = []
y = []
z = []
ymax = []
ymin = []

for (key,val) in query_to_ani.items():
    if len(val) >= 3:
        print(key)
        print(query_to_err[key])
        x.append(val[0])
        y.append(val[1])
        z.append(val[2])
        ymax.append(query_to_err[key][0] - val[1])
        ymin.append(np.abs(query_to_err[key][1] - val[1]))
plt.errorbar(x,y,yerr = [ymin,ymax], fmt = 'o', c = 'orange')
plt.scatter(x,z)
plt.plot([90,100],[90,100],'--', c = 'black')
plt.show()


