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

ab_t = 0.01


def read_metaphlan(file, cut):
    species_abs = dict()
    with open(file,'r') as f:
        for line in f:
            if 's__' in line:
                ab = float(line.split('\t')[1])
                spec = line.split('\t')[0].split('|')[-1]
                if cut:
                    if ab > ab_t:
                        species_abs[spec] = ab
                else:
                    if ab < ab_t:
                        species_abs[spec] = ab
    return species_abs


def mp_plot(file1, file2, sample, cut):
    c = 'black'
    if cut == False:
        c = 'red'

    # Read the files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = read_metaphlan(file2, cut)
    #print(df1,df2)

    # Filter the first file for the specific sample
    df1 = df1[df1['Sample_file'] == 'real_gut_reads/' + sample + '_1.fastq.gz']
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

    print(sample, 'sylph', num_gn_df1, 'metaphlan', num_gn_df2, 'diff', num_gn_df1 - num_gn_df2)

    if cut:
        a = ax[0]
    else:
        a = ax[1]
    a.set_ylabel('Sylph species detected')
    a.set_xlabel('MetaPhlAn4 species detected')
    a.plot(num_gn_df2, num_gn_df1, 'o', c = c, ms = 3)
    #plt.title(f'Log-Log Plot of True Coverage vs Mean for {sample}')
    #plt.xscale('log')
    #plt.yscale('log')
    a.grid(True)

def get_sample(file):
    return file.split('/')[-1].split('_')[0]

fig, ax = plt.subplots(2,1, figsize=(14*cm, 12*cm))
# Usage
file1 = sys.argv[1]  # Replace with the path to your first file
file2 = sys.argv[2:]  # Replace with the path to your second file
for d in [True, False]:
    for f in file2:
        sample = get_sample(f)
        #print(sample)
        mp_plot(file1, f, sample, d)
ax[0].plot([0,300],[0,300],'-')
ax[1].plot([0,300],[0,300],'-')
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)

ax[0].legend(['Species with > 0.01% abundance'],frameon=False)
ax[1].legend(['Species with <= 0.01% abundance'],frameon=False,)
ax[0].set_title("Sylph vs MetaPhlAn4 number of species detected on 50 real gut metagenomes")
plt.savefig("supp_figs/metaphlan-vs-sylph_50_gut.svg")
plt.show()

