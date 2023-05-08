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
mash_results = []
sour_results = []


prita_files = [
        "/home/jshaw/projects/prita_test/results/mock2_ill_c100",
        "/home/jshaw/projects/prita_test/results/mock2_ill_c1000",
        "/home/jshaw/projects/prita_test/results/mock2_nano_c100",
        "/home/jshaw/projects/prita_test/results/mock2_nano_c1000",
        "/home/jshaw/projects/prita_test/results/mock2_pac_c100",
        "/home/jshaw/projects/prita_test/results/mock2_pac_c1000",]

mash_files = [
        "/home/jshaw/projects/prita_test/results/mash_ill",
        "/home/jshaw/projects/prita_test/results/mash_nano",
        "/home/jshaw/projects/prita_test/results/mash_pac",]

sourmash_files = [
        "/home/jshaw/projects/prita_test/results/sour_ill",
        "/home/jshaw/projects/prita_test/results/sour_nano",
        "/home/jshaw/projects/prita_test/results/sour_pac",]

for file in mash_files:
    mash_results.append([])
    for line in open(file,'r'):
        ani = line.split()[0]
        mash_results[-1].append(float(ani) * 100)

for file in sourmash_files:
    sour_results.append([])
    for line in open(file,'r'):
        if 'ani' in line:
            continue
        ani = line.split(',')[-1]
        sour_results[-1].append(float(ani) * 100)


for file in prita_files:
    results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split()
        ref_file = spl[1]
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

fig, ax = plt.subplots(ncols = 3, figsize = (16* cm , 10 * cm), sharey = True)
res_anis = [[x.final_ani for x in res] for res in results]
res_anis_low = [[x.final_ani for x in res if x.low] for res in results]
res_anis_pass = [[x.adj_ani for x in res if not x.low] for res in results]
boxes = []
all_res = []
s = 7

offset = 0.0
width = 0.6
labels = []


for (i,x) in enumerate(['Illumina', 'nanopore-old', 'PacBio']):
    ax[i].set_title(x)
    positions = []
    positions.append(i*4 - 1)
    positions.append(i*4 + 0 - offset)
    positions.append(i*4 + 1 - 2 * offset)
    positions.append(i*4 + 2 - 3 * offset)

    boxes = []
    box_c100 = res_anis[2*i]
    box_c1000 = res_anis[2*i+1]
    box_mash = mash_results[i]
    box_sour = sour_results[i]

    boxes.append(box_c100)
    boxes.append(box_c1000)
    boxes.append(box_mash)
    boxes.append(box_sour)
    
    pos_c100_low = rand_jitter([positions[0]  for x in range(len(res_anis_low[2*i]))])
    pos_c100_pass = rand_jitter([positions[0]  for x in range(len(res_anis_pass[2*i]))])
    pos_c1000_low = rand_jitter([positions[1]  for x in range(len(res_anis_low[2*i + 1]))])
    pos_c1000_pass = rand_jitter([positions[1]  for x in range(len(res_anis_pass[2*i+ 1]))])
    pos_mash = rand_jitter([positions[2]  for x in range(len(box_mash))])
    pos_sour = rand_jitter([positions[3]  for x in range(len(box_sour))])

    dot_label = []
    for box in boxes:
        geq99 = 0
        geq95 = 0
        for val in box:
            if val >= 99:
                geq99 += 1;
                geq95 += 1
            elif val >= 95:
                geq95 += 1
        print(geq99, geq95)
        dot_label.append(str(geq99) + "/87, " + str(geq95) + "/87")


    ax[i].scatter(pos_c100_low,res_anis_low[2*i], s = s, c = 'black', marker = "s")
    ax[i].scatter(pos_c100_pass,res_anis_pass[2*i],  s = s, c = cmap[0], label = dot_label[0])

    ax[i].scatter( pos_c1000_low,res_anis_low[2*i+1], s = s, c = 'black', marker = "s")
    ax[i].scatter( pos_c1000_pass, res_anis_pass[2*i+1],s = s, c = cmap[4], label=dot_label[1])

    ax[i].scatter( pos_mash,box_mash,s = s, c = cmap[2], label=dot_label[2])
    ax[i].scatter( pos_sour,box_sour, s = s, c = cmap[3], label=dot_label[3])
    ax[i].legend(frameon = False,title="# ANI > 99, 95")

    labels = []
    labels.append("prita -c 100")
    labels.append("prita -c 1000")
    labels.append("mash screen")
    labels.append("sourmash")
    ax[i].boxplot(boxes, showfliers=False, positions = positions, widths = width, labels=labels)
    ax[i].tick_params(axis='x', labelrotation=60)



for a in ax:
    a.spines[['right', 'top']].set_visible(False)

    #plt.scatter(rand_jitter([positions[4*i]  for x in range(len(res_anis_low[2*i]))]), res_anis_low[2*i], s = s, c = 'black')
    #plt.scatter(rand_jitter([positions[4*i]  for x in range(len(res_anis_pass[2*i]))]), res_anis_pass[2*i], s = s)

    #plt.scatter(rand_jitter([positions[4*i+1]  for x in range(len(res_anis_low[2*i+1]))]), res_anis_low[2*i+1], s = s, c = 'black')
    #plt.scatter(rand_jitter([positions[4*i+1]  for x in range(len(res_anis_pass[2*i+1]))]), res_anis_pass[2*i+1], s = s)

    #plt.scatter(rand_jitter([positions[4*i+2]  for x in range(len(box_mash))]), box_mash, s = s)
    #plt.scatter(rand_jitter([positions[4*i+3]  for x in range(len(box_sour))]), box_sour, s = s)

#ax.boxplot(boxes, showfliers=False, positions = positions, widths = width, vert=False)
plt.show()

 
