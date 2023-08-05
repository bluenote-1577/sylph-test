import matplotlib.pyplot as plt
import scipy.stats as stats

import statsmodels.api as sm
import statsmodels.stats as ss

from scipy.stats.distributions import norm,uniform

import matplotlib.cm as cm
import numpy as np
cmap = plt.get_cmap('tab20')
plt.set_cmap(cmap)

def fdr(p_vals, alpha ):
    not_used, pvals = ss.multitest.fdrcorrection(p_vals, alpha = alpha)
    return pvals


cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
fig, ax = plt.subplots(figsize = (6.5* cm , 6.5 * cm))

pvals = []
c = []
s = []
ind = 0
derep = False
seen_cs = set()
mag_to_pval = dict()
mag_to_effect = dict()
#pval_file = './no-adj_ordered_pvals.txt'
#pval_file = './ordered_pvals.txt'
#pval_file = './98-ordered_pvals.txt'
pval_file = './98-fixed-covs-ordered_pvals.txt'
#pval_file = './95-fixed-covs-ordered_pvals.txt'
order_file = './121ordered_mags_by_cluster.txt'
metadata = './genomes-all_metadata.tsv'
mag_to_c = dict()
first = True
gn_to_rep = dict()
gn_to_mag_rep = dict()
for line in open(metadata,'r'):
    if first:
        first = False
        continue
    spl = line.split()
    gn = spl[0]
    rep = spl[14] +" " + spl[15]
    gn_to_rep[gn] = rep
    gn_to_mag_rep[gn] = spl[13]
    #print(rep)


s_and_mag_pair = []
for line in open(pval_file,'r'):
    spl = line.split()
    if spl[ind] == "nan":
        continue
    #print(spl)
    if derep:
        if spl[4] not in seen_cs:
            c.append(int(spl[4])%20)
            #pvals.append(np.log10(float(spl[ind])))
            s.append(float(spl[ind]))
            mag_to_c[spl[-1].rstrip()] = int(spl[4])
            mag_to_pval[spl[-1].rstrip()] = np.log10(float(spl[ind]))
            seen_cs.add(spl[4])
    else:
        c.append(int(spl[4])%20)
        mag_to_c[spl[-1].rstrip()] = int(spl[4])
        #pvals.append(np.log10(float(spl[ind])))
        mag_to_pval[spl[-1].rstrip()] = np.log10(float(spl[ind]))
        mag_to_effect[spl[-1].rstrip()] = float(spl[2])
        if spl[4] not in seen_cs:
            s.append(float(spl[ind]))
            s_and_mag_pair.append([float(spl[ind]),spl[-1][0:-3]])
            seen_cs.add(spl[4])


qq = stats.probplot(s,dist=uniform)
s_and_mag_pair = sorted(s_and_mag_pair,reverse=False)
for i in range(10):
    print(gn_to_rep[s_and_mag_pair[i][1]])
    #exit()

ax.scatter(-np.log10(qq[0][0]), -np.log10(qq[0][1]), s = 4, label = 'Arbitrary species representative')
#ax.plot([np.min(-np.log10(qq[0][0])), np.max(-np.log10(qq[0][0]))], [np.min(-np.log10(qq[0][1])), np.max(-np.log10(qq[0][1]))], 'r-')
ax.plot([0,5],[0,5], 'r-')
#plt.title('Q-Q plot of -log10(p-values) against uniform distribution')
plt.ylabel('Expected P-value (-log10 scale)')
plt.xlabel('Observed P-value (-log10 scale)')
#plt.grid(True)
ax.spines[['right', 'top']].set_visible(False)
plt.legend( frameon=False)

plt.savefig("figures/qqplot.png", dpi = 300)
plt.show()
#exit()

s = sorted(pvals, reverse=True) 
count = 0 
val = 0
q = 0.05
for p in s:
    if  (len(pvals) - count) / len(pvals) * q > 10**p:
        val = p
        #print(10**val, q)
        break

cs = []
mag_pval_list = []
seen_cs = set()
for line in open(order_file,'r'):
    mag = line.rstrip()
    if mag in mag_to_pval:
        pvals.append(mag_to_pval[mag])
        seen_cs.add(mag_to_c[mag])
        cs.append(len(seen_cs)%20)
        mag_pval_list.append(mag)

fig, ax = plt.subplots(figsize = (16* cm , 6 * cm))
agatho_qvals = []
agatho_rep = 'MGYG000002492'
#print(pvals)
qvals = fdr(np.power(10,pvals), q)
qvals = np.log10(qvals)
seen_reps = set()
for (i,pval) in enumerate(qvals):
    m = mag_pval_list[i][:-3]
    if pval < np.log10(q):
        rep = gn_to_rep[m]
        if rep not in seen_reps:
            s = f"{10**pvals[i]},{10**pval},{rep},{gn_to_mag_rep[m]},{m},{mag_to_effect[mag_pval_list[i]]})"
            print(s)
        seen_reps.add(rep)
    if gn_to_mag_rep[m] == agatho_rep:
        agatho_qvals.append(pval)
#print(qvals)
if derep:
    size = 3
else:
    size = 0.3
plt.scatter(range(len(qvals)), -np.array(qvals),c = cs,s = size)
ax.spines[['right', 'top']].set_visible(False)
plt.xticks([])
plt.axhline(-np.log10(q))

#plt.scatter(range(len(qvals)), fd,c = c,s = 1)
plt.ylabel("-log10(q-val)")
plt.xlabel("UHGG MAGs coloured and clustered by species")
plt.savefig('figures/manhat.png', dpi = 300)
plt.show()

cmap = plt.get_cmap('tab20')
plt.set_cmap(cmap)

fig, ax = plt.subplots(figsize = (6.5* cm , 6.5 * cm))
plt.ylabel("-log10(q-val)")
plt.xlabel("Agathobacter rectalis")
mag_c = mag_to_c[agatho_rep + '.fa']
#print(mag_c)
it_tab20 = [plt.cm.tab20(i) for i in range(20)]
plt.scatter(range(len(agatho_qvals)), -np.array(agatho_qvals), c = [it_tab20[mag_c+8] for x in agatho_qvals], s = size, cmap = cmap)
plt.xticks([])
ax.spines[['right', 'top']].set_visible(False)
plt.axhline(-np.log10(q))
plt.savefig("figures/agatho.png", dpi = 300)
plt.show()
