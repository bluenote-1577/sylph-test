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
cmap = sns.color_palette("deep")
plt.rcParams.update({'font.size': 6.5})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
results = []
mash_results = []
sour_results = []

sylph_files = [
        "results_meslier-oct15/mock2_ill_c200",
        "results_meslier-oct15/mock2_nano_c200",
        "results_meslier-oct15/mock2_pac_c200",
        ]


#sylph_files = [
#        "/home/jshaw/projects/sylph_test/results/mock2_ill_c100",
#        "/home/jshaw/projects/sylph_test/results/mock2_ill_c1000",
#        "/home/jshaw/projects/sylph_test/results/mock2_nano_c100",
#        "/home/jshaw/projects/sylph_test/results/mock2_nano_c1000",
#        "/home/jshaw/projects/sylph_test/results/mock2_pac_c100",
#        "/home/jshaw/projects/sylph_test/results/mock2_pac_c1000",]

mash_files = [
        "results_oct15/mash_ill",
        "results_oct15/mash_nano",
        "results_oct15/mash_pac",]

sourmash_files = [
        "results_oct15/sourmash_ill",
        "results_oct15/sourmash_nano",
        "results_oct15/sourmash_pac",]

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


for file in sylph_files:
    results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split()
        ref_file = spl[1]
        query_file = spl[0]
        final_ani = float(spl[2])
        naive_ani = float(spl[10])
        adj_ani = None
        low = False
        if "LOW" in spl[5]:
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
        mean_cov = float(spl[8])
        median_cov = float(spl[7])
        res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam, 0, final_ani, low)
        results[-1].append(res)

fig, ax = plt.subplots(ncols = 3, figsize = (16* cm, 11 * cm * 0.5), sharey = True)
res_anis = [[x.final_ani for x in res] for res in results]
res_anis_low = [[x.final_ani for x in res if x.low] for res in results]
res_anis_pass = [[x.adj_ani for x in res if not x.low] for res in results]
boxes = []
all_res = []
s = 7

offset = 0.0
width = 0.6
labels = []
num_methods = 3


for (i,x) in enumerate(['Illumina', 'Nanopore-old', 'PacBio']):
    ax[i].set_title(x)
    positions = []
    positions.append(i*num_methods - 1)
    positions.append(i*num_methods + 0 - offset)
    positions.append(i*num_methods + 1 - 2 * offset)

    boxes = []
    box_c100 = res_anis[i]
    box_mash = mash_results[i]
    box_sour = sour_results[i]

    boxes.append(box_c100)
    boxes.append(box_mash)
    boxes.append(box_sour)
    
    pos_c100_low = rand_jitter([positions[0]  for x in range(len(res_anis_low[i]))])
    pos_c100_pass = rand_jitter([positions[0]  for x in range(len(res_anis_pass[i]))])
    pos_mash = rand_jitter([positions[1]  for x in range(len(box_mash))])
    pos_sour = rand_jitter([positions[2]  for x in range(len(box_sour))])

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
        #dot_label.append(str(geq99) + "/87\n " + str(geq95) + "/87")
        dot_label.append(str(geq99) + "\n" + str(geq95))


    print(res_anis_low[i])
    ax[i].scatter(pos_c100_low,res_anis_low[i], s = s, color = 'black', marker = "s")
    #ax[i].scatter(pos_c100_pass,res_anis_pass[i],  s = s, color= cmap[0], label = dot_label[0] + ' sylph query')
    ax[i].scatter(pos_c100_pass,res_anis_pass[i],  s = s, color= cmap[0])
   # ax[i].set_ylim([70,100])

    #ax[i].scatter( pos_c1000_low,res_anis_low[2*i+1], s = s, color = 'black', marker = "s")
    #ax[i].scatter( pos_c1000_pass, res_anis_pass[2*i+1],s = s, color = cmap[4], label=dot_label[1] + ' sylph -c1000')
    #ax[i].scatter( pos_c1000_pass, res_anis_pass[2*i+1],s = s, color = cmap[4])

    #ax[i].scatter( pos_mash,box_mash,s = s, color = cmap[2], label=dot_label[2] + ' mash screen')
    ax[i].scatter( pos_mash,box_mash,s = s, color = cmap[1])

    #ax[i].scatter( pos_sour,box_sour, s = s, color = cmap[3], label=dot_label[3] + ' sourmash')
    ax[i].scatter( pos_sour,box_sour, s = s, color = cmap[2])
    #ax[i].legend(frameon = False,title="# ANI > 99, 95")

    labels = []
    labels.append("sylph query\n\n" + dot_label[0])
    #labels.append("sylph\n-c 1000\n\n" + dot_label[1])
    labels.append("mash screen\n\n" + dot_label[1])
    labels.append("sourmash\n\n" + dot_label[2])
    bp = ax[i].boxplot(boxes, showfliers=False, positions = positions, widths = width, labels=labels)
    for median in bp['medians']:
        median.set_color('black')
        print(median.get_ydata())
    ax[i].tick_params(axis='x', labelrotation=0)

    if i == 0:
        ax[i].set_ylabel('Query ANI')



for a in ax:
    a.spines[['right', 'top']].set_visible(False)

    #plt.scatter(rand_jitter([positions[4*i]  for x in range(len(res_anis_low[2*i]))]), res_anis_low[2*i], s = s, c = 'black')
    #plt.scatter(rand_jitter([positions[4*i]  for x in range(len(res_anis_pass[2*i]))]), res_anis_pass[2*i], s = s)

    #plt.scatter(rand_jitter([positions[4*i+1]  for x in range(len(res_anis_low[2*i+1]))]), res_anis_low[2*i+1], s = s, c = 'black')
    #plt.scatter(rand_jitter([positions[4*i+1]  for x in range(len(res_anis_pass[2*i+1]))]), res_anis_pass[2*i+1], s = s)

    #plt.scatter(rand_jitter([positions[4*i+2]  for x in range(len(box_mash))]), box_mash, s = s)
    #plt.scatter(rand_jitter([positions[4*i+3]  for x in range(len(box_sour))]), box_sour, s = s)

#ax.boxplot(boxes, showfliers=False, positions = positions, widths = width, vert=False)
plt.savefig("figures/mock_community_box.pdf")
plt.show()

 
