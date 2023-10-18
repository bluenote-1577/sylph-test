import sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from collections import defaultdict
import seaborn as sns
from dataclasses import dataclass
from natsort import natsorted
cmap = sns.color_palette("muted")
plt_diag = True

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
plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


results = []
true_results = []

sylph_files_gap = [
        "gtdb-on-reads-oct15/gtdb-on-ill-c200.tsv",
        "gtdb-on-reads-oct15/gtdb-on-nano-c200.tsv",
        "gtdb-on-reads-oct15/gtdb-on-pac-c200.tsv",
        ]


truth_files = [
        "gtdb-on-reads/true-gtdb-on-on-mock-c100.tsv",
        #"gtdb-on-reads/true-gtdb-on-on-mock.tsv",
        ]

mash_files = [
        "gtdb-gap/gtdb-on-reads/gtdb-on-ill-mash-s1000.tsv",
        "gtdb-gap/gtdb-on-reads/gtdb-on-nano-mash-s1000.tsv",
        "gtdb-gap/gtdb-on-reads/gtdb-on-pac-mash-s1000.tsv",
        ]

sylph_files = sylph_files_gap

mash_results = []
for file in mash_files:
    mash_results.append([])
    for line in open(file,'r'):
        ani = line.split()[0]
        file = line.split()[4].split('/')[-1]
        mash_results[-1].append((file,float(ani) * 100))

for file in sylph_files:
    results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split()
        ref_file = spl[1].split('/')[-1]
        query_file = spl[0]
        final_ani = float(spl[2])
        naive_ani = float(spl[10])
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
        mean_cov = float(spl[8])
        median_cov = float(spl[7])
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


fig, ax = plt.subplots(1, 3, figsize = (16* cm , 7 * cm), sharey = True, sharex = True)
s = 8

for (i,name) in enumerate(['Illumina', 'Nanopore-old', 'PacBio']):
    for j in range(1):
        query_to_ani = defaultdict(list)
        query_to_err = defaultdict(list)
        query_to_anis_nn = defaultdict(list)
        query_to_diff = defaultdict(list)
        query_to_eff_cov = dict()
        for res in true_results[j]:
            if res.final_ani > 90:
                query_to_anis_nn[res.ref_file].append(res.final_ani)
        for (key,res) in query_to_anis_nn.items():
            query_to_ani[key].append(np.max(res))
        for res in results[i]:
            #if not res.low and res.lam != None:
            if not res.low:
                query_to_eff_cov[res.ref_file] = res.lam
                if res.lam != None:
                    query_to_ani[res.ref_file].append(res.final_ani)
                else:
                    query_to_ani[res.ref_file].append(res.naive_ani)
                query_to_err[res.ref_file].append(res.high_ani)
                query_to_err[res.ref_file].append(res.low_ani)
                query_to_diff[res.ref_file].append(res.final_ani - query_to_ani[res.ref_file][0])
        for res in mash_results[i]:
            query_to_ani[res[0]].append(res[1])
            query_to_diff[res[0]].append(res[1] - query_to_ani[res[0]][0])

        x = []
        x_lam = []
        y = []
        z = []
        ymax = []
        ymin = []

        y_diff = []
        z_diff = []

        covered = 0
        total = 0
        for (key,val) in query_to_ani.items():
            if len(val) >= 3:
                x.append(val[0])
                y.append(val[1])
                z.append(val[2])
                if query_to_err[key][0] == None:
                    continue
                x_lam.append(query_to_eff_cov[key])
                y_diff.append(query_to_diff[key][0])
                z_diff.append(query_to_diff[key][1])
                if z_diff[-1] > y_diff[-1]:
                    print(key,val)
                ymax.append(query_to_err[key][0] - val[1])
                ymin.append(np.abs(query_to_err[key][1] - val[1]))
                if x_lam[-1] < 10:
                    if val[0] <= query_to_err[key][0] and val[0] >= query_to_err[key][1]:
                        covered += 1
                    total += 1
        #ax[j][i].errorbar(x,y,yerr = [ymin,ymax], fmt = 'o', c = 'orange', ms = 2)

        ax[i].spines[['right', 'top']].set_visible(False)
        if plt_diag:
            good_lr = stats.linregress(x,y)
            naive_lr = stats.linregress(x,z)

            ax[i].set_xlabel("True containment ANI")
            if j == 1:
                ax[i].scatter(x,y,s = s, color = cmap[2], alpha = 0.5, label = rf"c = 1000, $R^2$ = {round(good_lr.rvalue**2,3)}")
            else:
                if i == 0:
                    ax[i].set_title("Illumina", fontsize = plt.rcParams['font.size'])
                elif i == 1:
                    ax[i].set_title("Nanopore-old", fontsize = plt.rcParams['font.size'])
                elif i == 2:
                    ax[i].set_title("PacBio", fontsize = plt.rcParams['font.size'])

                ax[i].scatter(x,y,s = s, color = cmap[0], alpha = 0.5, label = rf"sylph, $R^2$ = {round(good_lr.rvalue**2,3)}")
            ax[i].scatter(x,z, s= s, color = cmap[3], alpha = 0.5, label = rf"Naive, $R^2$ = {round(naive_lr.rvalue**2,3)}")
            ax[i].plot([90,100],[90,100],'--', c = 'black')
            print('covered: ' + str(covered/total) + 'total: ' + str(total))
            ax[i].set_ylim([85,100])
            if i == 0:
                ax[i].set_ylabel("Estimated ANI")


        else:
            ax[i].set_xlabel("Estimated lambda")
            if j == 1:
                ax[i].scatter(x_lam,y_diff,s = s, color = cmap[2], alpha = 0.5, label = 'c = 1000')
            else:
                if i == 0:
                    ax[i].set_title("Illumina", fontsize = plt.rcParams['font.size'])
                elif i == 1:
                    ax[i].set_title("Nanopore-old", fontsize = plt.rcParams['font.size'])
                elif i == 2:
                    ax[i].set_title("PacBio", fontsize = plt.rcParams['font.size'])

                ax[i].scatter(x_lam,y_diff,s = s, color = cmap[0], alpha = 0.5, label = 'sylph')
            ax[i].scatter(x_lam,z_diff, s= s, color = cmap[3], alpha = 0.5, label = 'Naive')
            ax[i].set_xscale('log')
            ax[i].axhline(0, c = 'black', ls = '--')
            #ax[j][i].plot([90,100],[90,100],'--', c = 'black')
            print('covered: ' + str(covered/total) + 'total: ' + str(total))
            if i == 0:
                ax[i].set_ylabel("ANI deviation from truth")

        if plt_diag:
            lg = ax[i].legend(frameon = False, loc='upper left' )
        else:
            lg = ax[i].legend(frameon = False, loc='lower right' )
        #change the marker size manually for both lines
        lg.legendHandles[0]._sizes = [20]
        lg.legendHandles[1]._sizes = [20] 

            #ax[j][i].set_ylim([85,100])
if plt_diag:
    plt.savefig("figures/mock_diag.svg")
else:
    plt.savefig("figures/mock_deviation.svg")
plt.show()

    #plt.scatter(x_lam, y_diff)
    #plt.scatter(x_lam, z_diff)
    ##plt.xscale('log')
    #plt.show()


