import matplotlib.pyplot as plt
import scipy.stats as stats

import statsmodels.api as sm

from scipy.stats.distributions import norm,uniform

import matplotlib.cm as cm
import numpy as np
cmap = plt.get_cmap('tab20')
plt.set_cmap(cmap)

def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
fig, ax = plt.subplots(figsize = (8* cm , 8 * cm))

pvals = []
c = []
s = []
ind = 0
derep = False
seen_cs = set()
mag_to_pval = dict()
#pval_file = './no-adj_ordered_pvals.txt'
#pval_file = './ordered_pvals.txt'
#pval_file = './98-ordered_pvals.txt'
pval_file = './98-fixed-covs-ordered_pvals.txt'
#pval_file = './95-fixed-covs-ordered_pvals.txt'
order_file = './ordered_mags_by_cluster.txt'
mag_to_c = dict()

for line in open(pval_file,'r'):
    spl = line.split()
    #if spl[ind] == "nan":
    #    continue
    #print(spl)
    if derep:
        if spl[4] not in seen_cs:
            c.append(int(spl[4])%20)
            pvals.append(np.log10(float(spl[ind])))
            s.append(float(spl[ind]))
        seen_cs.add(spl[4])
    else:
        c.append(int(spl[4])%20)
        mag_to_c[spl[-1].rstrip()] = int(spl[4])
        #pvals.append(np.log10(float(spl[ind])))
        mag_to_pval[spl[-1].rstrip()] = np.log10(float(spl[ind]))
        if spl[4] not in seen_cs:
            s.append(float(spl[ind]))
            seen_cs.add(spl[4])


qq = stats.probplot(s,dist=uniform)

ax.scatter(-np.log10(qq[0][0]), -np.log10(qq[0][1]), s = 4, label = 'Species representative')
#ax.plot([np.min(-np.log10(qq[0][0])), np.max(-np.log10(qq[0][0]))], [np.min(-np.log10(qq[0][1])), np.max(-np.log10(qq[0][1]))], 'r-')
ax.plot([0,5],[0,5], 'r-')
#plt.title('Q-Q plot of -log10(p-values) against uniform distribution')
plt.ylabel('Expected P-value (-log10 scale)')
plt.xlabel('Observed P-value (-log10 scale)')
#plt.grid(True)
ax.spines[['right', 'top']].set_visible(False)
plt.legend( frameon=False)
plt.show()
#exit()

s = sorted(pvals, reverse=True) 
count = 0 
val = 0
q = 0.001
for p in s:
    if  (len(pvals) - count) / len(pvals) * q > 10**p:
        val = p
        print(10**val, q)
        break

cs = []
seen_cs = set()
for line in open(order_file,'r'):
    mag = line.rstrip()
    if mag in mag_to_pval:
        pvals.append(mag_to_pval[mag])
        seen_cs.add(mag_to_c[mag])
        cs.append(len(seen_cs)%20)
fig, ax = plt.subplots(figsize = (16* cm , 10 * cm))
plt.scatter(range(len(pvals)), -np.array(pvals),c = cs,s = 1)
ax.spines[['right', 'top']].set_visible(False)

#plt.scatter(range(len(pvals)), fd,c = c,s = 1)
plt.axhline(-np.log10(0.01))
plt.ylabel("-log10(p-val)")
plt.xlabel("UHGG MAGs coloured by species")
plt.xticks([])
plt.show()
