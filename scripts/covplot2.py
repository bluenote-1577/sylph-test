import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

def stem(x):
    return x.split('/')[-1]

def log_log_plot(file1, file2, sample):
    # Read the files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t', names=['Genome', 'Mean'])
    print(df1,df2)

    # Filter the first file for the specific sample
    df1 = df1[df1['Sample_file'].apply(stem) == sample]

    # Create a mapping from Genome_file to Genome
    df1['Genome'] = df1['Genome_file'].apply(lambda x: x.split('/')[-1].replace('.fna.gz', ''))

    # Merge the two dataframes on the Genome column
    merged_df = pd.merge(df1, df2, on='Genome')

    # Plotting
    plt.plot(pd.to_numeric(merged_df['Mean']), pd.to_numeric(merged_df['True_cov']), 'o', ms = 2, c = 'black', alpha = 0.2)
    #plt.plot([np.min(pd.to_numeric(merged_df['True_cov']))], [np.max(pd.to_numeric(merged_df['True_cov']))])
    plt.ylabel('Sylph coverage (log scale)')
    plt.xlabel('Alignment coverage (log scale)')
    plt.title(f'Log-Log Plot of True Coverage vs Mean for {sample}')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)

def get_sample(file):
    with open(file,'r') as f:
        first = f.readline()
        return first.split()[1]

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


plt.figure(figsize=(12*cm, 12*cm))
# Usage
file1 = sys.argv[1]  # Replace with the path to your first file
file2 = sys.argv[2:]  # Replace with the path to your second file
for f in file2:
    sample = get_sample(f)
    print(sample)
    log_log_plot(file1, f, sample)
plt.loglog([0.01,500],[0.01,500],'-')
plt.show()

plt.savefig('supp_figs/real_gut_cov_v0.5.svg')
