import sys


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from collections import defaultdict
import seaborn as sns
from dataclasses import dataclass
from natsort import natsorted
from scipy import stats

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
    eff_cov: float
    final_ani: float
    low: bool
    contig_name: str

read_file_to_pair_dict = dict()
read_file_to_status= dict()

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

results = []
#metadata = '/home/jshaw/projects/sylph_test/atopic_controls_table.txt'
metadata = 'atopic_controls_table.txt'
#metadata = '/home/jshaw/projects/sylph_test/all_case_control_table_nofail.txt'
#metadata = '/home/jshaw/projects/sylph_test/spt_case.txt'

boxplot_covs = True

prev_file = None
for line in open(metadata,'r'):
    file = line.split()[0]
    case_control = line.split()[1].rstrip()
    if prev_file == None:
        prev_file = file
    else:
        read_file_to_pair_dict[prev_file] = file
        read_file_to_pair_dict[file] = prev_file
        prev_file = None
    read_file_to_status[file] = case_control

prita_files = [
        #"/home/jshaw/projects/sylph_test/results/c100-fungi+aureus.tsv",
        #"results/c100-fungi-jul2.txt",
        "results_oct15/chng-fungi.tsv",
        ]

for file in prita_files:
    results.append([])
    for line in open(file,'r'):
        if 'Naive' in line:
            continue
        spl = line.split('\t')
        ref_file = spl[1].split('/')[-1]
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
        contig_nam = spl[-1]
        res = result(mean_cov, adj_ani, naive_ani, median_cov, ref_file, query_file, cis[0], cis[1], lam, float(spl[3]), final_ani, low, contig_nam)
        results[-1].append(res)


pair_to_res = defaultdict(list)
pair_to_res_naive = defaultdict(list)

pair_to_res_globo = defaultdict(list)
pair_to_res_naive_globo = defaultdict(list)

case_res = []
control_res = []


case_res_n = []
control_res_n = []
res_used_pairs = set()

case_globo = []
control_globo = []

case_globo_n = []
control_globo_n = []

globo_used_pairs = set()


fig, ax = plt.subplots(2,2,figsize = (8* cm , 8 * cm))
for res in results[0]:
    #srr = res.query_file.split('/')[-1].split('_')[0]
    srr = res.query_file.split('/')[-1].split('.')[0].split('_')[0]
    if srr not in read_file_to_pair_dict:
        continue
    pair = natsorted([srr,read_file_to_pair_dict[srr]])
    if 'restrict' in res.contig_name:
        pair_to_res[tuple(pair)].append((res.final_ani))
        pair_to_res_naive[tuple(pair)].append(res.naive_ani)
        if tuple(pair) in res_used_pairs:
            continue
        if read_file_to_status[srr] == 'Case':
            if not boxplot_covs:
                case_res.append(res.final_ani)
                case_res_n.append(res.naive_ani)

            else:
                case_res.append(res.mean_cov)
                case_res_n.append(res.mean_cov)
        else:
            if not boxplot_covs:
                control_res.append(res.final_ani)
                control_res_n.append(res.naive_ani)

            else:
                control_res.append(res.mean_cov)
                control_res_n.append(res.mean_cov)

        res_used_pairs.add(tuple(pair))
    if 'globo' in res.contig_name:
        pair_to_res_globo[tuple(pair)].append((res.final_ani))
        pair_to_res_naive_globo[tuple(pair)].append(res.naive_ani)

        if tuple(pair) in globo_used_pairs:
            continue
        if read_file_to_status[srr] == 'Case':
            if not boxplot_covs:
                case_globo.append(res.final_ani)
                case_globo_n.append(res.naive_ani)

            else:
                case_globo.append(res.mean_cov)
                case_globo_n.append(res.mean_cov)

        else:
            if not boxplot_covs:
                control_globo.append(res.final_ani)
                control_globo_n.append(res.naive_ani)
            else:
                control_globo.append(res.mean_cov)
                control_globo_n.append(res.mean_cov)

        globo_used_pairs.add(tuple(pair))

#print(len(globo_used_pairs))
print(list(pair_to_res.values()))

sc = np.array(list(pair_to_res.values()))
sc2 = np.array(list(pair_to_res_naive.values()))

sc_g = np.array(list(pair_to_res_globo.values()))
sc2_g = np.array(list(pair_to_res_naive_globo.values()))

s = 7
ax[0][0].scatter(sc[:,0], sc[:,1], alpha = 1.00, s = s, color = cmap[0], label = 'sylph adjusted')
ax[0][1].scatter(sc2[:,0], sc2[:,1], alpha = 1.00, s = s, color = cmap[3], label =  'naive containment')
ax[1][0].scatter(sc_g[:,0], sc_g[:,1], alpha = 1.00, s = s, color = cmap[0], label = 'sylph adjusted')
ax[1][1].scatter(sc2_g[:,0], sc2_g[:,1], alpha = 1.00, s = s, color = cmap[3], label =  'naive containment')
ax[0][0].plot([82,100],[82,100],'--',c='black')
ax[1][0].plot([82,100],[82,100],'--',c='black')
ax[0][1].plot([82,100],[82,100],'--',c='black')
ax[1][1].plot([82,100],[82,100],'--',c='black')
ax[0][0].set_ylabel('Right sample ANI')
ax[1][0].set_ylabel('Right sample ANI')
ax[1][1].set_xlabel('Left sample ANI')
ax[1][0].set_xlabel('Left sample ANI')


lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes[0:2]]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines, labels, loc = 'upper center', frameon=False,
           bbox_transform = plt.gcf().transFigure, ncol = 2)

sr = stats.linregress(sc[:,0], sc[:,1])
nr = stats.linregress(sc2[:,0], sc2[:,1])
sg = stats.linregress(sc_g[:,0], sc_g[:,1])
ng = stats.linregress(sc2_g[:,0], sc2_g[:,1])

ax[0][0].text(.05, .99, rf"M. restricta ", ha='left', va='top', transform=ax[0][0].transAxes)
ax[0][1].text(.05, .99, 'M. restricta', ha='left', va='top', transform=ax[0][1].transAxes)
ax[1][1].text(.05, .99, 'M. globosa', ha='left', va='top', transform=ax[1][1].transAxes)
ax[1][0].text(.05, .99, 'M. globosa', ha='left', va='top', transform=ax[1][0].transAxes)


ax[0][0].text(.55, .25, rf"$R^2$ = {round(sr.rvalue**2,3)}", ha='left', va='top', transform=ax[0][0].transAxes)
ax[0][1].text(.55, .25, rf"$R^2$ = {round(nr.rvalue**2,3)}", ha='left', va='top', transform=ax[0][1].transAxes)
ax[1][1].text(.55, .25, rf"$R^2$ = {round(ng.rvalue**2,3)}", ha='left', va='top', transform=ax[1][1].transAxes)
ax[1][0].text(.55, .25, rf"$R^2$ = {round(sg.rvalue**2,3)}", ha='left', va='top', transform=ax[1][0].transAxes)

for a in ax:
    for b in a:
        b.spines[['right', 'top']].set_visible(False)
        #b.legend(frameon=False, loc='lower right')
#fig.tight_layout()
fig.tight_layout(rect=(0,0,1,0.95))
plt.savefig("figures/chng_left_right_concordance.pdf")
plt.show()


from scipy.stats import ranksums
#fig, ax = plt.subplots(2,2, figsize = (8* cm , 8 * cm))
#ax[0][0].boxplot([case_res, control_res])
#ax[0][1].boxplot([case_globo, control_globo])
#ax[1][0].boxplot([case_res_n, control_res_n])
#ax[1][1].boxplot([case_globo_n, control_globo_n])
#
#plt.show()
#
#print(len(case_res), len(control_res))
#print(ranksums(case_res,control_res))
#print(ranksums(case_globo,control_globo))
#
#print(ranksums(case_res_n,control_res_n))
#print(ranksums(case_globo_n,control_globo_n))

plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


def add_stat_annotation(ax, bp, pval):
    # Determine the height of the annotation based on the y-axis limits
    ylim = ax.get_ylim()
    y_range = ylim[1] - ylim[0]
    y_annotation = ylim[1] - 0.00 * y_range

    # Determine the x positions of the two boxes being compared
    #x1 = bp["boxes"][0].get_xdata().mean()
    #x2 = bp["boxes"][1].get_xdata().mean()
    x1 = 1
    x2 = 2

    # Add the asterisk if the p-value is significant
    scale = 0.085
    if pval < 0.001:
        ax.plot([x1, x1, x2, x2], [y_annotation, y_annotation + scale * y_range, y_annotation + scale * y_range, y_annotation], lw=1.5, c='k')
        ax.text((x1+x2)/2, y_annotation + 0.01 * y_range, "***", ha='center', va='bottom', fontsize = 11)
    elif pval < 0.01:
        ax.plot([x1, x1, x2, x2], [y_annotation, y_annotation + scale * y_range, y_annotation + scale * y_range, y_annotation], lw=1.5, c='k')
        ax.text((x1+x2)/2, y_annotation + 0.01 * y_range, "**", ha='center', va='bottom', fontsize = 11)
    elif pval < 0.05:
        ax.plot([x1, x1, x2, x2], [y_annotation, y_annotation + scale * y_range, y_annotation + scale * y_range, y_annotation], lw=1.5, c='k')
        ax.text((x1+x2)/2, y_annotation + 0.01 * y_range, "*", ha='center', va='bottom', fontsize = 11)



    # Add the p-value as text
    ax.text((x1+x2)/2, y_annotation - 0.50, "p={:.4f}".format(pval), ha='center', va='bottom', fontsize=7)

# Create the figure and axes
fig, ax = plt.subplots(2,2, figsize=(8*cm, 8*cm))

# Add the boxplots to the axes
bp1 = ax[0][0].boxplot([case_res, control_res], patch_artist=True, boxprops=dict(facecolor=cmap[0]), labels = ['case', 'control'])
bp2 = ax[0][1].boxplot([case_res_n, control_res_n], patch_artist=True, boxprops=dict(facecolor=cmap[3]),labels = ['case', 'control'])
bp3 = ax[1][0].boxplot([case_globo, control_globo], patch_artist=True, boxprops=dict(facecolor=cmap[0]), labels = ['case','control'])
bp4 = ax[1][1].boxplot([case_globo_n, control_globo_n], patch_artist=True, boxprops=dict(facecolor=cmap[3]), labels = ['case', 'control'])
boxplots = [bp1, bp2, bp3, bp4]

print(f"n = {len(case_res) + len(case_globo)}")

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes[0:2]]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines, labels, loc = 'upper center', frameon=False,
           bbox_transform = plt.gcf().transFigure, ncol = 2)



# Set the x and y axis labels for the boxplots
ax[0][0].set_ylabel('ANI')
ax[0][1].set_ylabel('ANI')
ax[1][1].set_ylabel('ANI')
ax[1][0].set_ylabel('ANI')

#for l in range(2):
    #for k in range(2):
        #ax[l][k].set_ylim([85,103])

ax[0][0].text(.05, 1.15, rf"M. restricta ", ha='left', va='top', transform=ax[0][0].transAxes, fontsize = 8)
ax[0][1].text(.05, 1.15, 'M. restricta', ha='left', va='top', transform=ax[0][1].transAxes, fontsize = 8)
ax[1][1].text(.05, 1.15, 'M. globosa', ha='left', va='top', transform=ax[1][1].transAxes, fontsize = 8)
ax[1][0].text(.05, 1.15, 'M. globosa', ha='left', va='top', transform=ax[1][0].transAxes, fontsize = 8)



# Set the titles for the subplots
#ax[0][0].set_title('M. restricta ANI sylph')
#ax[0][1].set_title('M. restricta ANI naive')
#ax[1][0].set_title('M. globosal ANI sylph')
#ax[1][1].set_title('M. globosa ANI naive')

p_values = [ranksums(case_res, control_res)[1],  ranksums(case_res_n, control_res_n)[1], ranksums(case_globo, control_globo)[1], ranksums(case_globo_n, control_globo_n)[1]]
print(p_values)

add_stat_annotation(ax[0][0], bp1, p_values[0])
add_stat_annotation(ax[0][1], bp2, p_values[1])
add_stat_annotation(ax[1][0], bp3, p_values[2])
add_stat_annotation(ax[1][1], bp4, p_values[3])

bps = [bp1,bp2,bp3,bp4]
for bp in bps:
    for median in bp['medians']:
        median.set_color('black')

for a in ax:
    for b in a:
        b.spines[['right', 'top']].set_visible(False)

fig.legend([bp1["boxes"][0], bp2["boxes"][0]], ['sylph adjusted', 'naive containment'], loc='upper center', frameon=False, bbox_transform = plt.gcf().transFigure, ncol = 2)
fig.tight_layout(rect=(0,0,1,0.95))
plt.savefig("figures/chng_adjusted_ani_pvalues.pdf")
plt.show()
