import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})



# Read the files
skani_df = pd.read_csv('skani-test_D2-skani.tsv', sep='\t')
contain_df = pd.read_csv('skani-test_D2_contain.tsv', sep='\t')
# Process the 'Ref_file' and 'Query_file' to match the format in contain_df
skani_df['Ref_file_processed'] = skani_df['Ref_file'].apply(lambda x: x.split('/')[-1])
skani_df['Query_file_processed'] = skani_df['Query_file'].apply(lambda x: x.split('/')[-1])
skani_df = skani_df[skani_df['ANI'] > 90]  # Filter out ANI < 95
# Process the 'Sample_file' and 'Genome_file' to match the format in skani_df
contain_df['Sample_file'] = contain_df['Sample_file'].apply(lambda x: x.split('/')[-1])
contain_df['Genome_file'] = contain_df['Genome_file'].apply(lambda x: x.split('/')[-1])

# Find corresponding entries
corresponding_entries = pd.merge(skani_df, contain_df, left_on=['Ref_file_processed', 'Query_file_processed'],
                                 right_on=['Sample_file', 'Genome_file'])

# Calculate average Aligned Fraction for the first file and ANI difference
corresponding_entries['Average_Align_Fraction'] = (corresponding_entries['Align_fraction_ref'] + corresponding_entries['Align_fraction_query']) / 2
corresponding_entries['ANI_Difference'] = corresponding_entries['Naive_ANI'] - corresponding_entries['ANI']

# Plot
plt.figure(figsize=(8*cm, 5*cm))
#plt.scatter(corresponding_entries['Average_Align_Fraction'], corresponding_entries['ANI_Difference'],s=3, alpha=0.5, color='black')
plt.scatter(corresponding_entries['ANI'], corresponding_entries['Naive_ANI'],s=3, alpha=0.5, color='black')
plt.xlabel('Aligned fraction')
plt.ylabel('containment ANI\nminus skani ANI')
#dashed
#plt.axhline(0, color='red', linestyle='--', linewidth=1.0)
plt.plot([90, 100], [90, 100], color='red', linestyle='--', linewidth=1.0)
plt.grid()
sns.despine()
plt.title('Average aligned fraction versus ANI (k=31)\nB. anthracis genomes (skani ANI > 90)')
plt.tight_layout()
plt.savefig('figures/D2_ani_vs_af.svg')
plt.show()

