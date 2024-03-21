import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
cm = 1/2.54  # centimeters in inches\n",
cmap = sns.color_palette("hls")
plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


# Replace these with the paths to your files
#codes = ['ERR7739005', 'ERR7745346','ERR7803603','ERR7745335','ERR7745291']
#codes = ['ERR7745291']
#file_path1 = './ERR7745346_subsamp0.10_1.fastq.sylphmpa'
#file_path2 = './ERR7745346_1.fastq.gz.sylphmpa'


def read_and_filter(file_path):
    # Reading the file
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df.fillna('NA', inplace=True)
    df.rename(columns={df.columns[0]: 'Species', df.columns[1]: 'Relative Abundance'}, inplace=True)
    species_df = df[df[df.columns[0]].str.contains('s__') & ~df[df.columns[0]].str.contains('t__')]
    if 'sylphmpa' in file_path:
        species_df['Species'] = species_df['Species'].str.replace('|', ';', regex=False)
    if 'motus' in file_path:
        species_df['Relative Abundance'] *= 100

    sorted_df = species_df.sort_values(by='Relative Abundance', ascending=False)  # Sort the DataFrame in descending order
    top_quantile_percent_df = sorted_df
    return top_quantile_percent_df

# find codes in the folders metaphlan_results, sylph_results, motus_results and intersect them
import os
metaphlan_codes = [f.split('_')[0] for f in os.listdir('metaphlan_results')]
sylph_codes = [f.split('_')[0] for f in os.listdir('sylph_results_round1')]
motus_codes = [f.split('_')[0] for f in os.listdir('motus_results')]
codes = list(set(metaphlan_codes) & set(sylph_codes) & set(motus_codes))
print(codes)

plt.figure(figsize=(10, 6))
corrs = []
for code in codes:
    #file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
    #file_paths2 = [f'metaphlan_results/{code}_subsamp0.10_metaphlan.txt.gtdb', f'sylph_results/{code}_subsamp0.10_1.fastq.sylphmpa', f'motus_results/{code}_subsamp0.10_motus.txt.motus']

    file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results_round1/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
    file_paths2 = [f'sylph_results_round1/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus', f'metaphlan_results/{code}_metaphlan.txt.gtdb']


    for file_path1, file_path2 in zip(file_paths1, file_paths2):
        # Reading and filtering data
        species_abundance1 = read_and_filter(file_path1).set_index('Species')
        species_abundance2 = read_and_filter(file_path2).set_index('Species')

        # Merging dataframes by index (Species ID)
        merged_df = species_abundance1.join(species_abundance2, lsuffix='_File1', rsuffix='_File2', how='outer')

        # Replacing NaN values with 0 (assuming absence of species)
        merged_df.fillna(0, inplace=True)

        # Log-log plot
        plt.scatter(merged_df['Relative Abundance_File1'], merged_df['Relative Abundance_File2'], label=file_path1, s=15,alpha=1.0)
        # print pearson correlation between relative abundance
        a = merged_df['Relative Abundance_File1']
        print(len(a),file_path1)
        b = merged_df['Relative Abundance_File2']
        #spearman corr
        #spearman = a.corr(b, method='spearman')
        spearman = a.corr(b)
        #print(a.corr(b), file_path1, file_path2)
        #l1 dist
        #spearman = abs(a-b).mean()
        #corrs.append(l1)
        corrs.append(spearman)

    #plt.plot([1, 1], [0, 1], color='black', linestyle='--')
#    plt.xlim(1e-4, 100)
#    plt.ylim(1e-4, 100)
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlabel('Relative Abundance in File 1')
#    plt.ylabel('Relative Abundance in File 2')
#    plt.title(f'Correlation between {file_path1} and {file_path2} (Spearman: {spearman:.2f})')
#    plt.grid(True)
#    plt.legend()
#    plt.show()
#
#every third
#plt.figure(figsize=(10, 6))
plt.figure(figsize=(6*cm, 5*cm))
sylph_metaphlan = corrs[0::3]
sylph_motus = corrs[1::3]
motus_metaphlan = corrs[2::3]

#print median
print('sylph_metaphlan', pd.Series(sylph_metaphlan).median())
print('sylph_motus', pd.Series(sylph_motus).median())
print('motus_metaphlan', pd.Series(motus_metaphlan).median())
#annotate median values
sns.stripplot(data=[sylph_metaphlan, sylph_motus, motus_metaphlan], palette='hls', alpha=0.8, marker='o', size=3)

sns.boxplot(data=[sylph_metaphlan, sylph_motus, motus_metaphlan],boxprops=dict(alpha=.5), showfliers=False, palette='hls')
plt.text(0+.125, 1.05, f"{pd.Series(sylph_metaphlan).median():.2f}", ha='right', va='center', fontsize=7)
plt.text(1+.125, 1.05, f"{pd.Series(sylph_motus).median():.2f}", ha='right', va='center', fontsize=7)
plt.text(2+.125, 1.05, f"{pd.Series(motus_metaphlan).median():.2f}", ha='right', va='center', fontsize=7)


plt.xticks([0, 1, 2], ['sylph\nMetaPhlAn4', 'sylph\nmOTUs3', 'mOTUs3\nMetaPhlAn4'])
plt.title("Abundance correlation between methods\n(50 gut metagenomes)", fontsize=7)
plt.ylabel("Pearson correlation")
sns.despine()
plt.savefig('supp_figs/correlation_methods.svg')
plt.show()
