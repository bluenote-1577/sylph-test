import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

ab_t = 0.00


def read_metaphlan(file, cut, genus=True):
    species_abs = dict()
    with open(file,'r') as f:
        for line in f:
            if not genus:
                if 's__' in line:
                    ab = float(line.split('\t')[1])
                    spec = line.split('\t')[0].split('|')[-1]
                    if cut:
                        if ab > ab_t:
                            species_abs[spec] = ab
                    else:
                        if ab < ab_t:
                            species_abs[spec] = ab
            else:
                if 's__' in line and 't__' not in line:
                    ab = float(line.split('\t')[1])
                    if 'motus' in file:
                        spec = line.split('\t')[0].split(';')[-2]
                    else:
                        spec = line.split('\t')[0].split('|')[-2]
                    print(spec)
                    if cut:
                        if ab > ab_t:
                            species_abs[spec] = ab
                    else:
                        if ab < ab_t:
                            species_abs[spec] = ab
    return species_abs


def mp_plot(file1, file2, sample, cut, motus = False, first = False):
    c = 'black'
    if cut == False:
        c = 'red'

    # Read the files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = read_metaphlan(file2, cut)
    #print(df1,df2)

    # Filter the first file for the specific sample
    df1 = df1[df1['Sample_file'].str.contains(sample) & ~df1['Sample_file'].str.contains('subsamp')]
    if cut:
        df1 = df1[df1['Taxonomic_abundance'] > ab_t]
    else:
        df1 = df1[df1['Taxonomic_abundance'] < ab_t]
    #print(df1)

    # Create a mapping from Genome_file to Genome
    df1['Genome'] = df1['Genome_file'].apply(lambda x: x.split('/')[-1].replace('.fna.gz', ''))

    # Merge the two dataframes on the Genome column
    num_gn_df1 = len(df1)
    num_gn_df2 = len(df2)
    if num_gn_df1 == 0 or num_gn_df2 == 0:
        return

    if motus:
        print(sample, 'sylph', num_gn_df1, 'mOTUs', num_gn_df2, 'diff', num_gn_df1 - num_gn_df2)
    else:
        print(sample, 'sylph', num_gn_df1, 'metaphlan', num_gn_df2, 'diff', num_gn_df1 - num_gn_df2)

    if cut:
        a = ax
    else:
        a = ax
    a.set_ylabel('sylph species detected')
    a.set_xlabel('MetaPhlAn4/mOTUs3 species detected')
    if motus:
        c = 'blue'
        marker = 'x'
        l = 'mOTUs3'
    else:
        c = 'red'
        marker = 'o'
        l = 'MetaPhlAn4'
    if first:
        a.plot(num_gn_df2, num_gn_df1, marker , c = c, ms = 3, label = l)
    else:
        a.plot(num_gn_df2, num_gn_df1, marker , c = c, ms = 3)
    #plt.title(f'Log-Log Plot of True Coverage vs Mean for {sample}')
    #plt.xscale('log')
    #plt.yscale('log')
    #a.grid(True)

def get_sample(file):
    return file.split('/')[-1].split('_')[0]

fig, ax = plt.subplots(1,1, figsize=(4*cm, 4*cm))
# Usage
file1 = sys.argv[1]  # Replace with the path to your first file
file2 = sys.argv[2:]  # Replace with the path to your second file
motfirst = True
metfirst = True
for d in [True]:
    for f in file2:
        if 'subsamp' in f:
            continue
        if 'motus' in f:
            motus = True
        else:
            motus = False
        sample = get_sample(f)
        print(sample)
        if motfirst and motus:
            mp_plot(file1, f, sample, d, motus, first = True)
            motfirst = False
            continue
        if metfirst and not motus:
            mp_plot(file1, f, sample, d, motus, first = True)
            metfirst = False
            continue
        mp_plot(file1, f, sample, d, motus)
        # set the legend for the first plot only
ax.plot([100,800],[100,800],'--', c = 'black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(frameon=False)
#ax[1].plot([0,400],[0,400],'-')
#ax[1].spines['top'].set_visible(False)
#ax[1].spines['right'].set_visible(False)
#ax[1].legend(['Species with <= 0.01% abundance'],frameon=False,)

#ax.set_title("# detected species (50 real gut metagenomes)")
plt.tight_layout()
plt.savefig("supp_figs/metaphlan-vs-sylph_50_gut.svg")
plt.show()

