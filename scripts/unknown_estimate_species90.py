import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 8})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})



def sum_second_column(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, nrows=50)
    return df.iloc[:, 1].sum()

def sum_sequence_abundance(unknown_df, identifier):
    filtered_df = unknown_df[unknown_df['Sample_file'].str.split('_').str[1] == identifier]
    return filtered_df['Sequence_abundance'].sum()

# Read the species90-unknown-est-sylph.tsv file
unknown_df = pd.read_csv('profile_results/species90-unknown-est-sylph.tsv', sep='\t')

# Get all species90_x.tsv files
tsv_files = glob.glob('profile_results/species90_*.tsv')

sums_tsv = []
sums_fqgz = []

for file in tsv_files:
    # Extract identifier (assuming it's the part after the last underscore and before .tsv)
    identifier = os.path.basename(file).split('_')[-1].split('.')[0]
    print(identifier)

    sum_tsv = sum_second_column(file) * 100
    sum_fqgz = sum_sequence_abundance(unknown_df, identifier)
    print(sum_fqgz)

    sums_tsv.append(sum_tsv)
    sums_fqgz.append(sum_fqgz)

# Scatter plot
import statsmodels.api as sm
figure = plt.figure(figsize=(10*cm, 10*cm))
plt.scatter(sums_tsv, sums_fqgz)
plt.xlabel('True sequence abundance of species in database')
plt.ylabel('Sylph sequence abundance estimate (with option -u)')
plt.title('% of reads detected in database\nfor the 50 known species')
#get pearson correlation and p value
r = pd.Series(sums_tsv).corr(pd.Series(sums_fqgz))

# Add a constant to the input data to include an intercept in the regression model
x_values_with_constant = sm.add_constant(sums_tsv)

# Create a model and fit it
model = sm.OLS(sums_fqgz, x_values_with_constant)
results = model.fit()

# Get the p-value
p_value = results.pvalues[1]  # p-value for the slope coefficient

print(f"P-value for the regression: {p_value}")

plt.plot([0, 50], [0, 50], color='red', linestyle='--', label = 'y = x')
# plot the linaer regression 
#show p-value exp format
plt.plot(sums_tsv, results.predict(), color='green', label = f"Pearson R: {r:.2f},  p-value: {p_value:.2e}")
plt.legend(frameon=False)
sns.despine()
plt.tight_layout()
plt.savefig('./supp_figs/unknown_estimate_species90.pdf')
plt.show()
