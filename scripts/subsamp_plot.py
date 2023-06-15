import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import defaultdict
from dataclasses import dataclass

#fig, axs = plt.subplots(1,3,sharey = True, sharex = True)  # create a figure with 3 subplots

#cmap = sns.color_palette("hls", 4)
cm = 1/2.54  # centimeters in inches\n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
fig, axs = plt.subplots(1, 3, figsize = (16* cm , 8 * cm), sharey = True, sharex = True)

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

files = ['./real_prita_results/13712885_downsample_results_85.txt','./real_prita_results/5131958_downsample_results_85.txt','./real_prita_results/214966_downsample_results_85.txt']
for (i,file) in enumerate(files):
    ref_file_to_true_ani = dict()
    ref_file_to_true_mean_cov = dict()
    ref_file_to_res = defaultdict(list)
    re_str = "(\d+\.?\d+)-.+\.fastq"
    cutoff = 5
    lam_cutoff = 0.00

    start = True
    for line in open(file,'r'):
        if start:
            start = False
            continue
        x = re.findall(re_str, line)
        if len(x) == 0:
            spl = line.split()
            ref_file = spl[1]
            true_ani = float(spl[2])
            est_cov = float(spl[3])
            median_cov = float(spl[8])
            mean_cov = float(spl[9])
            if median_cov > cutoff and true_ani > 97:# and est_cov - median_cov < median_cov :
                #print(line)
                ref_file_to_true_ani[ref_file] = true_ani
                ref_file_to_true_mean_cov[ref_file] = mean_cov
        else:
            spl = line.split()
            ref_file = spl[1]
            query_file = spl[0]
            naive_ani = float(spl[3])
            adj_ani = None
            if spl[5] == "NA" or "NA" in spl[4]:
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
            res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam)
            ref_file_to_res[ref_file].append(res)

    mc_x = []
    diff_y = []
    naive_y = []
    max_diffs = []
    num_in = 0
    num_out = 0

    #print(ref_file_to_true_ani)

    for (ref,reses) in ref_file_to_res.items():
        if ref in ref_file_to_true_ani:
            for r in reses:
                if r.adj_ani is not None:
                    if r.lam < lam_cutoff:
                        continue
                    true_ani = ref_file_to_true_ani[ref]
                    true_mc = ref_file_to_true_mean_cov[ref] 
                    mc_x.append(r.mean_cov)
                    diff_y.append(r.adj_ani - true_ani)
                    max_diffs.append([np.abs(r.adj_ani- true_ani), ref, r.query_file])
                    naive_y.append(r.naive_ani - true_ani)
                    if r.high_ani is not None:
                        if true_ani <= r.high_ani and true_ani >= r.low_ani:
                            num_in += 1
                        else:
                            #print(r)
                            num_out += 1

    max_diffs.sort(reverse = True)
    print(max_diffs[0:20])
    print(num_in, num_out, num_in / (num_in + num_out+1))
    axs[i].plot(mc_x, diff_y, 'o', ms = 5, alpha = 0.3, c =  cmap[0])
    axs[i].plot(mc_x, naive_y, 'o', ms = 5, alpha = 0.3, c = cmap[3])
    from sklearn import datasets, linear_model
    axs[i].axhline(0, ls= '--', c = 'black' )
    if i == 0:
        axs[i].set_title("Human gut illumina")
        axs[i].set_ylabel("ANI deviation after subsampling")
    elif i == 1:
        axs[i].set_title("Ocean illumina")
        axs[i].set_xlabel("sylph estimated lambda")
    else:
        axs[i].set_title("Human saliva nanopore")

    regr = linear_model.LinearRegression()
    regr.fit(np.array(np.log(mc_x)).reshape(-1,1), np.array(diff_y))
    xvals = np.linspace(np.min(mc_x),3,100)
    y_pred = regr.predict(np.array(np.log(xvals)).reshape(-1,1))
    axs[i].plot(xvals, y_pred)

    regr = linear_model.LinearRegression()
    regr.fit(np.array(np.log(mc_x)).reshape(-1,1), np.array(naive_y))
    xvals = np.linspace(np.min(mc_x),3,100)
    y_pred = regr.predict(np.array(np.log(xvals)).reshape(-1,1))
    axs[i].plot(xvals, y_pred)
    axs[i].set_xscale('log')
plt.show()







