import sys
import seaborn as sns
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from collections import defaultdict
cmap = sns.color_palette("muted")
cm = 1/2.54  # centimeters in inches\n",
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

cm = 1/2.54

coverm = sys.argv[-1]
sylph = sys.argv[1]

genome_to_cov  = defaultdict(list)
df1 = pd.read_csv(coverm, sep='\t')
df2 = pd.read_csv(sylph, sep='\t')

# Define function to process Genome_file to Genome
def process_genome_file(genome_file):
    return genome_file.split('/')[-1].split('.f')[0]
def process_sample(x):
    return x.split('/')[-1].split('_')[0]
def process_low(x):
    if x == 0 or x == None:
        return 0
    spl = x.split('-')
    if len(spl) > 1:
        return float(x.split('-')[0]) * 1.45
    else:
        return 0
def process_high(x):
    if x == 0 or x == None:
        return 0
    spl = x.split('-')
    if len(spl) > 1:
        return float(x.split('-')[1]) * 1.45
    else:
        return 0

def process_ci(x):
    if x == None or x == 0:
        return 0
    else:
        return [float(y)*1.45 for y in x.split('-')]


cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
fig, axes = plt.subplots(1,1,figsize=(7 * cm, 7* cm))


# Process the Genome_file in df2 to make it comparable with df1
df2['Genome'] = df2['Genome_file'].apply(process_genome_file)
df2['Sample_file'] = df2['Sample_file'].apply(process_sample)

# Merge DataFrames
merged_df = pd.merge(df1, df2, on='Genome',how='left').fillna({"True_cov":0.001}).fillna(0)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(merged_df)
print(merged_df['True_cov'])

# Iterate over each Sample_file
sample_files = merged_df['Sample_file'].unique()
maxx = 0
labels = ['sylph with deduplication', 'sylph without deduplication']
print(sample_files)
i = 0
for sample in sample_files:
    if sample == 0:
        continue

    # Filter merged_df for each Sample_file
    sample_df = merged_df[merged_df['Sample_file'] == sample]

    # Extract corresponding columns
    eff_cov_column = 'Sequence_abundance'  # Assuming Eff_cov corresponds to 'Taxonomic_abundance'

    sample_df['Lambda_5-95_percentile'].mask(sample_df['Lambda_5-95_percentile'] == 'NA-NA', None, inplace=True)
    ci1 = sample_df['Lambda_5-95_percentile'].apply(process_high)
    ci2 = sample_df['Lambda_5-95_percentile'].apply(process_low)

    eff_cov = sample_df[eff_cov_column] + 0.000
    eff_cov = eff_cov / np.sum(eff_cov) * 100
    sample_mean = sample_df.iloc[:,1] + 0.000
    sample_mean = sample_mean / np.sum(sample_mean) * 100
    #sample_mean = sample_df.iloc[:,1]/2
    if max(sample_mean) > maxx:
        maxx = max(sample_mean)

    l = labels[i] + f'\n{len(eff_cov)} species detected'
    plt.scatter(sample_mean, eff_cov, label=l, s=8, c = cmap[i])
    #    plt.scatter(sample_mean, ci1, s=12, label = '95% upper conf')
    #    plt.scatter(sample_mean, ci2, s=12, label = '5% low conf')

    plt.ylabel('Sylph abundance (%)')
    plt.xlabel('Read alignment abundance (%)')
    #plt.title(f'')
    plt.xscale('log')
    plt.yscale('log')
    i+=1

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.plot([0,maxx],[0,maxx], '--', c = 'black')
plt.legend(frameon=False)
plt.savefig('../supp_figs/' + sys.argv[1] + '.svg')
plt.show()

