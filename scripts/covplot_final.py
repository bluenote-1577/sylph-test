import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats
import seaborn as sns
def stem(x):
    return x.split('/')[-1]


    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
cm = 1/2.54  # centimeters in inches\n",
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


ratio = False

vec_x = []
vec_y = []

def log_log_plot(file1, file2, sample, vec_x, vec_y, first):
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
    merged_df_H = merged_df[merged_df['Eff_lambda'] == 'HIGH']
    merged_df_L = merged_df[merged_df['Eff_lambda'] != 'HIGH']

    rel_ab = pd.to_numeric(merged_df['Mean']) 
    rel_ab_sylph = pd.to_numeric(merged_df['True_cov'])

    rel_ab_H = pd.to_numeric(merged_df_H['Mean'])
    rel_ab_sylph_H = pd.to_numeric(merged_df_H['True_cov'])

    rel_ab_L = pd.to_numeric(merged_df_L['Mean']) 
    rel_ab_sylph_L = pd.to_numeric(merged_df_L['True_cov']) 



    vec_x += list(rel_ab)
    vec_y += list(rel_ab_sylph)
    al = 0.3

    # Plotting
    if ratio:
        plt.plot(rel_ab, rel_ab_sylph / rel_ab, 'o', ms = 2, c = 'black', alpha = al)
    else:
        if first:
            plt.plot(rel_ab_L, rel_ab_sylph_L, 'o', ms = 2, c = 'blue', alpha = al, label = 'Estimated from statistical model')
            plt.plot(rel_ab_H, rel_ab_sylph_H, 'o', ms = 2, c = 'red', alpha = al, label = 'Not estimated from statistical model')
        else:
            plt.plot(rel_ab_L, rel_ab_sylph_L, 'o', ms = 2, c = 'blue', alpha = al)
            plt.plot(rel_ab_H, rel_ab_sylph_H, 'o', ms = 2, c = 'red', alpha = al)

    #plt.plot([np.min(pd.to_numeric(merged_df['True_cov']))], [np.max(pd.to_numeric(merged_df['True_cov']))])
    plt.ylabel('Sylph coverage (log scale)')
    plt.xlabel('Alignment coverage (log scale)')
    #plt.title(f'Log-Log Plot of True Coverage vs Mean for {sample}')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)

def get_sample(file):
    with open(file,'r') as f:
        first = f.readline()
        return first.split()[1]

plt.figure(figsize=(12*cm, 12*cm))
# Usage
file1 = sys.argv[1]  # Replace with the path to your first file
file2 = sys.argv[2:]  # Replace with the path to your second file
first = True
for f in file2:
    sample = get_sample(f)
    print(sample)
    log_log_plot(file1, f, sample, vec_x, vec_y, first)
    first = False

vec_x = np.array(vec_x)
vec_y = np.array(vec_y)
x = np.logspace(-2,2.6,100)
# Calculating the linear regression without an intercept
slope_no_intercept = np.sum(vec_x * vec_y) / np.sum(vec_x * vec_x)

# Creating the line of best fit without intercept
line_no_intercept = slope_no_intercept * x

r = stats.pearsonr(vec_x,vec_y)[0]
r = stats.pearsonr(vec_x,vec_y)[0]
spearman = stats.spearmanr(vec_x,vec_y)
print(spearman)
plt.loglog(x, x, '--', color='gray', label=f'Y = X')
#plt.loglog(x, line_no_intercept, '--', color='black', label=f'Y = {slope_no_intercept:.2f}X, R = {r:.2f}')
#    plt.plot([0.001,100],[0.001,100],'-', label=f'Pearson R = {r:.2f}')
plt.title("Read alignment versus sylph coverage on 50 real gut metagenomes")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend(frameon=False)
plt.savefig("supp_figs/read_align_50_gut_cov.svg")
plt.show()

