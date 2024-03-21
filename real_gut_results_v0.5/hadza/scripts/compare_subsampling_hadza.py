import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

only_motus = True

# Replace these with the paths to your files
codes = ['ERR7739005', 'ERR7745346','ERR7803603','ERR7745335','ERR7745291']
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


plt.figure(figsize=(5.5*cm, 5.5*cm))
corrs = []
l1s = []
for code in codes:
    file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
    file_paths2 = [f'metaphlan_results/{code}_subsamp0.10_metaphlan.txt.gtdb', f'sylph_results/{code}_subsamp0.10_1.fastq.gz.sylphmpa', f'motus_results/{code}_subsamp0.10_motus.txt.motus']

    #file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
    #file_paths2 = [f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus', f'metaphlan_results/{code}_metaphlan.txt.gtdb']




    for file_path1, file_path2 in zip(file_paths1, file_paths2):
        if only_motus and 'motus' not in file_path1:
            continue
        if not only_motus and 'motus' in file_path1:
            continue

        # Reading and filtering data
        species_abundance1 = read_and_filter(file_path1).set_index('Species')
        species_abundance2 = read_and_filter(file_path2).set_index('Species')
        print(species_abundance1)

        # Merging dataframes by index (Species ID)
        merged_df = species_abundance1.join(species_abundance2, lsuffix='_File1', rsuffix='_File2', how='outer')

        # Replacing NaN values with 0 (assuming absence of species)
        merged_df.fillna(0, inplace=True)

        # Log-log plot
        ty = file_path1.split('_')[0]
        spearman = merged_df['Relative Abundance_File1'].corr(merged_df['Relative Abundance_File2'], method='spearman')
        pearson = merged_df['Relative Abundance_File1'].corr(merged_df['Relative Abundance_File2'], method='pearson')
        #print l1 distance between 
        l1 = abs(merged_df['Relative Abundance_File1'] - merged_df['Relative Abundance_File2']).sum()


        #plt.scatter(merged_df['Relative Abundance_File1'], merged_df['Relative Abundance_File2'], label=ty + f'\nL1:{l1:.2f}',s=5,alpha=0.5, color =cmap[6])
        if 'sylph' in file_path1:
            c = cmap[0]
            marker = 'o'
        if 'metaphlan' in file_path1:
            c = cmap[4]
            marker = 'x'
        if 'motus' in file_path1:
            c = cmap[6]
            marker = 'o'
        plt.scatter(merged_df['Relative Abundance_File1'], merged_df['Relative Abundance_File2'], s=5,alpha=0.5, color = c, marker=marker)
        # print pearson correlation between relative abundance
        a = merged_df['Relative Abundance_File1']
        b = merged_df['Relative Abundance_File2']
        #spearman corr
        #spearman = a.corr(b, method='spearman')
        spearman = a.corr(b)
        #print(a.corr(b), file_path1, file_path2)
        corrs.append(spearman)
        l1s.append(l1)

        plt.plot([0, 1], [0, 1], color='black', linestyle='--')
        plt.xlim(1e-4, 10)
        plt.ylim(1e-4, 10)
        plt.xlabel('Abundance (no subsampling)')
        plt.ylabel('Abundance (subsampled reads)')
        #plt.grid(True)
        #plt.show()
if only_motus:
    mean_spearman = sum(corrs)/len(corrs)
    #make legend marker large
    plt.scatter(1,1, s=5,alpha=0.5, color = cmap[6], label = 'mOTUs3\nSpear. ' + f'{mean_spearman:.2f}')
else:
    mean_spearman_sy = sum(corrs[0::3])/len(corrs[0::3])
    mean_spearman_sm = sum(corrs[1::3])/len(corrs[1::3])
    plt.scatter(1,1, marker='o', s=5,alpha=0.5, color = cmap[0], label = 'sylph\nSpear. ' + f'{mean_spearman_sy:.2f}')
    plt.scatter(1,1, marker='x', s=5,alpha=0.5, color = cmap[4], label = 'MetaPhlAn4\nSpear. ' + f'{mean_spearman_sm:.2f}')
    #plt.annotate('sylph\nSpearman: ' + f'{mean_spearman_sy:.2f}', xy=(1,1), xytext=(1,1), textcoords='offset points', color=cmap[0])

sns.despine()
if only_motus:
    plt.legend(frameon=False,  fontsize=7, markerscale=2)
else:
    plt.legend(frameon=False,  fontsize=7, markerscale=2)

plt.xscale('log')
plt.yscale('log')
plt.tight_layout()

if only_motus:
    #plt.title('mOTUs')
    plt.savefig('hadza_figs/correlation_subsampling_motus.png', dpi=300)
else:
    #plt.title("sylph and MetaPhlAn4")
    plt.savefig('hadza_figs/correlation_subsampling.png', dpi=300)
plt.show()


#every third
#plt.figure(figsize=(10, 6))
sylph_metaphlan = corrs[0::3]
sylph_motus = corrs[1::3]
motus_metaphlan = corrs[2::3]

l1_sylph_metaphlan = l1s[0::3]
l1_sylph_motus = l1s[1::3]
l1_motus_metaphlan = l1s[2::3]

##mean
#print(sum(sylph_metaphlan)/len(sylph_metaphlan))
#print(sum(sylph_motus)/len(sylph_motus))
#print(sum(motus_metaphlan)/len(motus_metaphlan))
#
#print(sum(l1_sylph_metaphlan)/len(l1_sylph_metaphlan))
#print(sum(l1_sylph_motus)/len(l1_sylph_motus))
#print(sum(l1_motus_metaphlan)/len(l1_motus_metaphlan))
#
#
#sns.stripplot(data=[sylph_metaphlan, sylph_motus, motus_metaphlan])
#sns.boxplot(data=[sylph_metaphlan, sylph_motus, motus_metaphlan],boxprops=dict(alpha=.5))
#plt.show()
