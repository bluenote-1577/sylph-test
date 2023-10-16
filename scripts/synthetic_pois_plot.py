import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from collections import defaultdict
import seaborn as sns
from dataclasses import dataclass
import glob
from natsort import natsorted

np.random.seed(0)
def rand_jitter(arr):
    stdev = 0.00 
    return arr + np.random.randn(len(arr)) * stdev


cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
#cmap = sns.color_palette()
#cmap = sns.color_palette("Set2")
cmap = sns.color_palette("husl", 4)


plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

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


#21-500-10-0.48621887049186563.tsv
#files_g = [glob.glob('synthetic_simulated_results/31-*-*'), glob.glob('synthetic_simulated_results_95/31-*-*'), glob.glob('synthetic_simulated_results_85/31-*-*')]
#files_g = [glob.glob('synthetic_simulated_results_jul3/31-*-*'), glob.glob('synthetic_simulated_results_95_jul3/31-*-*'), glob.glob('synthetic_simulated_results_85_jul3/31-*-*')]
files_g = [glob.glob('synthetic_simulated_results_oct15/31-*-*'), glob.glob('synthetic_simulated_results_95_oct15/31-*-*'), glob.glob('synthetic_simulated_results_85_oct15/31-*-*')]

#files = sys.argv[2:]
true_anis = [100, 96.2, 90.5]
#true_ani = float(sys.argv[1])
#re_str = "(\d+)-(\d+\.?\d+)-.+\.fastq.gz"
re_str = "(\d+)-(\d+)-(\d+)-(\d+\.?\d+).tsv"
cov_plot = False

if not cov_plot: 
    fig, ax = plt.subplots(3,1, figsize = (14* cm , 12*cm))
else:
    fig, ax = plt.subplots(3,3, figsize = (16* cm , 16 * cm))

for index in range(3):
    cs = set()
    ks = set()
    ck_res = defaultdict(lambda: defaultdict(list))

    files = files_g[index]
    true_ani = true_anis[index]
    eff_covs = set()
    for file in natsorted(files):
        x = re.findall(re_str, file)

        k = int(x[0][0])
        c = int(x[0][1])
        cs.add(c)
        ks.add(k)
        it = int(x[0][2])
        abund = round(float(x[0][3]),4)
        eff_cov = round((150 - k + 1) / 150 * abund,4)
        eff_covs.add(eff_cov)

        for line in open(file,'r'):
            if 'Naive' in line:
                continue
            spl = line.split('\t')
            ref_file = spl[1]
            query_file = spl[0]
            naive_ani = float(spl[-2])
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
            mean_cov = float(spl[8])
            median_cov = float(spl[7])
            res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam, eff_cov)
            ck_res[c][k].append(res)

    cs = sorted(list(cs))
    ks = sorted(list(ks))
    print(cs,ks)
    
    k = 31
    naive_ani = []
    xs = []
    naive_xs = []
    ys = []
    ls = []
    cov_probs = []

    for c in cs:
        print(c)
        results = ck_res[c][k]
        x = [] 
        y = []
        l = []
        cov_prob = []
        for res in results:
            x.append(res.true_eff_cov)
            if res.adj_ani != None:
                y.append(res.adj_ani)
                l.append((res.true_eff_cov, res.lam))
                if true_ani <= res.high_ani and true_ani >= res.low_ani:
                    cov_prob.append(1)
                else:
                    cov_prob.append(0)
            else:
                y.append(res.naive_ani)
            if c == 100:
                naive_ani.append(res.naive_ani)
                naive_xs.append(res.true_eff_cov)

        print(np.mean(cov_prob))
        xs.append(x)
        ys.append(y)
        ls.append(l)
        cov_probs.append(cov_prob)
    

    
    xs.append(naive_xs)
    ys.append(naive_ani)
    #box_colors = [cmap[2], cmap[4], cmap[0]]
    boxes = []
    labels = []
    positions = []
    offset = 0.15
    width = 0.9
    s = 12

    

    if not cov_plot:

        xticks = []
        xticklabs = []
        pos_to_index = dict()
        for i in range(len(cs) + 1):
            xvals = sorted(list(set(xs[i])))
            for j,x in enumerate(xvals):
                pos_to_index[x] = j

            ypos = ys[i]
            center = len(cs)/2 * offset - i*offset
            xpos = np.array([pos_to_index[x] for x in xs[i]]) + center
            if i != len(cs):
                label = f"-c {cs[i]}"
            else:
                label = "Naive ANI"
            ax[index].scatter(xpos, ypos, s = s, color = cmap[i], label = label)
            print(xs[i])

        true_ani = true_anis[index]
        ax[index].axhline(y=100, linestyle='--', color='black')
        ax[index].axhline(y=true_ani, linestyle='--', color='red')
        #ax[index].legend(lines, labels, frameon=False, loc = 'lower right')
        ax[index].spines[['right', 'top']].set_visible(False)
        ax[index].set_xlabel("True effective coverage")
        if index == 0:
            ax[index].set_ylabel("Estimated ANI\n(Truth 100)")
        elif index == 1:
            ax[index].set_ylabel("Estimated ANI\n(Truth 96.2)")
        else:
            ax[index].set_ylabel("Estimated ANI\n(Truth 90.5)")

        ax[index].set_xticks(sorted(list(pos_to_index.values())))
        ax[index].set_xticklabels(sorted(list(pos_to_index.keys())))
        if index == 0:
            ax[index].legend(frameon=False)
        



        #xticks = [i*3+0.5 for i in range(len(grouped_data))]
        #ax[index].set_xticks(xticks)
        ##ax[index].set_ylim([72,101])
        #ax[index].set_xticklabels(sorted(grouped_data.keys()))
        ##ax.axhline(y=100, linestyle='--', color='gray')
        #ax[index].axhline(y=100, linestyle='--', color='black')
        #ax[index].axhline(y=true_ani, linestyle='--', color='red')

        ## Create dummy Line2D objects for the legend
        #lines = [plt.Line2D([0], [0], color=color, linewidth=3, linestyle='-') for color in box_colors]
        #labels = ['Naive containment ANI', 'sylph c = 1000', 'sylph c = 100']

        ## Add the legend
        #ax[index].legend(lines, labels, frameon=False, loc = 'lower right')
        #ax[index].spines[['right', 'top']].set_visible(False)
        #plt.xlabel("True effective coverage")
        #if index == 0:
        #    ax[index].set_ylabel("Estimated ANI\n(Truth 100)")
        #elif index == 1:
        #    ax[index].set_ylabel("Estimated ANI\n(Truth 96.2)")
        #else:
        #    ax[index].set_ylabel("Estimated ANI\n(Truth 90.5)")
    
    if cov_plot:
        r = [0.008,2.4]
        if index == 0:
            ani_label = ', ANI = 100'
        if index == 1:
            ani_label = ', ANI = 96.2'
        if index == 2:
            ani_label = ', ANI = 90.5'
        ax[0][index].scatter(np.array(l)[:,0], np.array(l)[:,1], color = cmap[2], s = 10, alpha = 1.0, facecolors='none',label='c = 1000' + ani_label)
        ax[0][index].plot(r,r, '--', c = 'black')

        ax[0][0].set_ylabel("Predicted effective coverage (lambda)")
        ax[1][0].set_ylabel("Predicted effective coverage (lambda)")
        ax[2][0].set_ylabel("Predicted effective coverage (lambda)")
        ax[2][0].set_xlabel("True effective coverage")
        ax[2][1].set_xlabel("True effective coverage")
        ax[2][2].set_xlabel("True effective coverage")

        ax[1][index].scatter(np.array(l500)[:,0], np.array(l500)[:,1], color = cmap[1], s = 10, alpha = 1.0, facecolors='none',label = 'c = 500' + ani_label)
        ax[1][index].plot(r,r, '--', c = 'black')
        #ax[1][index].yaxis.set_visible(False)
        ax[2][index].scatter(np.array(l100)[:,0], np.array(l100)[:,1], color = cmap[0], s = 10, alpha = 1.0, facecolors='none',label = 'c = 100' + ani_label)
        ax[2][index].plot(r,r,'--', c = 'black')
        #ax[2][index].yaxis.set_visible(False)
        for b in ax:
            for a in b:
                a.set_yscale('log')
                a.set_xscale('log')
                a.spines[['right', 'top']].set_visible(False)
        #plt.scatter(np.array(l100), 'o', c = cmap[2])
        ax[0][index].legend(frameon=False)
        ax[1][index].legend(frameon=False)
        ax[2][index].legend(frameon=False)
        #ax[0].set_xticks(sorted(grouped_data.keys()))
        #ax[0].set_xticklabels(sorted(grouped_data.keys()))
        #ax[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.tight_layout()
if cov_plot:
    plt.savefig("figures/coverage_est_plot.pdf")
else:
    plt.savefig("figures/synthetic_kleb_plot.pdf")

#plt.legend()
plt.show()
