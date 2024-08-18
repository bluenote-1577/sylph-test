import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("hls", 3)
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

genus = True
species = ~genus

# Replace these with the paths to your files
codes = ['ERR7739005', 'ERR7745346','ERR7803603','ERR7745335','ERR7745291']
sec = ['ERR7747615', 'ERR7746688', 'ERR7738937', 'ERR7738959', 'ERR7738938']
codes += sec
#codes = ['ERR7745291']
#file_path1 = './ERR7745346_subsamp0.10_1.fastq.sylphmpa'
#file_path2 = './ERR7745346_1.fastq.gz.sylphmpa'


def read_and_filter_genus(file_path):
    # Reading the file
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df.fillna('NA', inplace=True)
    df.rename(columns={df.columns[0]: 'Species', df.columns[1]: 'Relative Abundance'}, inplace=True)
    species_df = df[df[df.columns[0]].str.contains('g__') & ~df[df.columns[0]].str.contains('s__')]
    if 'motus' in file_path:
        species_df = df[df[df.columns[0]].str.contains('s__')]
        species_df['Species'] = species_df['Species'].str.split(';').str[:-1].apply(';'.join)
        species_df['Relative Abundance'] *= 100
        species_df = species_df.groupby('Species').sum().reset_index()
    if 'sylphmpa' in file_path:
        species_df['Species'] = species_df['Species'].str.replace('|', ';', regex=False)

    sorted_df = species_df.sort_values(by='Relative Abundance', ascending=False)  # Sort the DataFrame in descending order
    top_quantile_percent_df = sorted_df
    return top_quantile_percent_df




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


#plt.figure(3,1,figsize=(5.5*cm, 5.5*cm))
#plot with 3 subplots 3 columns


l1s = []
fig, axs = plt.subplots(1, 3, figsize=(12.1*cm, 4.5*cm))
sns.despine()
stats = {'l1': [], 'spearman': [], 'pearson': []}
filtered_stats = {'l1': [], 'spearman': [], 'pearson': []}
for code in codes:
#    file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
#    file_paths2 = [f'metaphlan_results/{code}_subsamp0.10_metaphlan.txt.gtdb', f'sylph_results/{code}_subsamp0.10_1.fastq.gz.sylphmpa', f'motus_results/{code}_subsamp0.10_motus.txt.motus']

    file_paths1 = [f'metaphlan_results//{code}_metaphlan.txt.gtdb', f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus']
    file_paths2 = [f'sylph_results/{code}_1.fastq.gz.sylphmpa', f'motus_results/{code}_motus.txt.motus', f'metaphlan_results/{code}_metaphlan.txt.gtdb']

    i = 0


    for file_path1, file_path2 in zip(file_paths1, file_paths2):

        # Reading and filtering data
        if genus:
            species_abundance1 = read_and_filter_genus(file_path1).set_index('Species')
            species_abundance2 = read_and_filter_genus(file_path2).set_index('Species')
        else:
            species_abundance1 = read_and_filter(file_path1).set_index('Species')
            species_abundance2 = read_and_filter(file_path2).set_index('Species')


        # Merging dataframes by index (Species ID)
        merged_df = species_abundance1.join(species_abundance2, lsuffix='_File1', rsuffix='_File2', how='outer')

        # Replacing NaN values with 0 (assuming absence of species)
        merged_df.fillna(0.00001, inplace=True)
        print(merged_df.head())

        #merged_df = merged_df[~((merged_df['Relative Abundance_File1'] < 0.0001) & (merged_df['Relative Abundance_File2'] < 0.0001))]
        merged_df_f = merged_df[(merged_df['Relative Abundance_File2'] > 0.00002) & (merged_df['Relative Abundance_File1'] > 0.00002)]

        #get 

        # Log-log plot
        ty = file_path1.split('_')[0]
        spearman = merged_df['Relative Abundance_File1'].corr(merged_df['Relative Abundance_File2'], method='spearman')
        pearson = merged_df['Relative Abundance_File1'].corr(merged_df['Relative Abundance_File2'], method='pearson')
        #print l1 distance between 
        l1 = abs(merged_df['Relative Abundance_File1'] - merged_df['Relative Abundance_File2']).sum()

        filtered_spearman = merged_df_f['Relative Abundance_File1'].corr(merged_df_f['Relative Abundance_File2'], method='spearman')
        filtered_pearson = merged_df_f['Relative Abundance_File1'].corr(merged_df_f['Relative Abundance_File2'], method='pearson')
        #print l1 distance between
        filtered_l1 = abs(merged_df_f['Relative Abundance_File1'] - merged_df_f['Relative Abundance_File2']).sum()


        #plt.scatter(merged_df['Relative Abundance_File1'], merged_df['Relative Abundance_File2'], label=ty + f'\nL1:{l1:.2f}',s=5,alpha=0.5, color =cmap[6])
        if i == 0:
            xl = 'MetaPhlAn4 abund. (%)'
            yl = 'sylph abund. (%)'
            if not genus:
                ty='MetaPhlAn4 vs sylph (s)'
            else:
                ty='MetaPhlAn4 vs sylph (g)'
            c = cmap[0]
            marker = 'o'
        if i == 1:
            xl = 'sylph abund. (%)'
            yl = 'mOTUs3 abund. (%)'
            c = cmap[1]
            marker = 'o'
            if not genus:
                ty='sylph vs mOTUs3 (s)'
            else:
                ty='sylph vs mOTUs3 (g)'
        if i == 2:
            xl = 'mOTUs3 abund. (%)'
            yl = 'MetaPhlAn4 abund. (%)'
            c = cmap[2]
            marker = 'o'
            if not genus:
                ty='mOTUs3 vs MetaPhlAn4 (s)'
            else:
                ty='mOTUs3 vs MetaPhlAn4 (g)'
        sns.scatterplot(x=merged_df['Relative Abundance_File1'], y=merged_df['Relative Abundance_File2'], s=7,alpha=0.5, color = c, marker=marker,ax=axs[i])
        # print pearson correlation between relative abundance
        a = merged_df['Relative Abundance_File1']
        b = merged_df['Relative Abundance_File2']
        #spearman = a.corr(b)
        #print(a.corr(b), file_path1, file_path2)
        #corrs.append(pearson)
        stats['l1'].append(l1)
        stats['spearman'].append(spearman)
        stats['pearson'].append(pearson)
        filtered_stats['l1'].append(filtered_l1)
        filtered_stats['spearman'].append(filtered_spearman)
        filtered_stats['pearson'].append(filtered_pearson)

        axs[i].plot([0, 10], [0, 10], color='black', linestyle='--')
        axs[i].set_xlim(1e-5, 100)
        axs[i].set_ylim(1e-5, 100)
        axs[i].set_xlabel(xl)
        axs[i].set_ylabel(yl)
        axs[i].set_title(ty, fontsize=7)
        if i > 0:
            axs[i].set_yticklabels([])
        #axs[i].set_xscale('symlog',linthresh=0.1)
        #axs[i].set_yscale('symlog',linthresh=0.1)
        axs[i].set_xscale('log')
        axs[i].set_yscale('log')

        i += 1
        #plt.grid(True)
        #plt.show()
    
mean_spearman_sy = sum(stats['spearman'][0::3])/len(stats['spearman'][0::3])
mean_spearman_sm = sum(stats['spearman'][1::3])/len(stats['spearman'][1::3])
mean_spearman_mm = sum(stats['spearman'][2::3])/len(stats['spearman'][2::3])

mean_pearson_sy = sum(stats['pearson'][0::3])/len(stats['pearson'][0::3])
mean_pearson_sm = sum(stats['pearson'][1::3])/len(stats['pearson'][1::3])
mean_pearson_mm = sum(stats['pearson'][2::3])/len(stats['pearson'][2::3])

mean_l1_sy = sum(stats['l1'][0::3])/len(stats['l1'][0::3])
mean_l1_sm = sum(stats['l1'][1::3])/len(stats['l1'][1::3])
mean_l1_mm = sum(stats['l1'][2::3])/len(stats['l1'][2::3])

filtered_mean_spearman_sy = sum(filtered_stats['spearman'][0::3])/len(filtered_stats['spearman'][0::3])
filtered_mean_spearman_sm = sum(filtered_stats['spearman'][1::3])/len(filtered_stats['spearman'][1::3])
filtered_mean_spearman_mm = sum(filtered_stats['spearman'][2::3])/len(filtered_stats['spearman'][2::3])

filtered_mean_pearson_sy = sum(filtered_stats['pearson'][0::3])/len(filtered_stats['pearson'][0::3])
filtered_mean_pearson_sm = sum(filtered_stats['pearson'][1::3])/len(filtered_stats['pearson'][1::3])
filtered_mean_pearson_mm = sum(filtered_stats['pearson'][2::3])/len(filtered_stats['pearson'][2::3]) 

filtered_mean_l1_sy = sum(filtered_stats['l1'][0::3])/len(filtered_stats['l1'][0::3])
filtered_mean_l1_sm = sum(filtered_stats['l1'][1::3])/len(filtered_stats['l1'][1::3])
filtered_mean_l1_mm = sum(filtered_stats['l1'][2::3])/len(filtered_stats['l1'][2::3])
axs[0].scatter(1,1, s=0,alpha=0.5, color = cmap[0], label = f'L1, Spear., Pear.\n{mean_l1_sy/100:.2f}, {mean_spearman_sy:.2f}, {mean_pearson_sy:.2f}')
axs[1].scatter(1,1, s=0,alpha=0.5, color = cmap[1], label = f'L1, Spear., Pear.\n{mean_l1_sm/100:.2f}, {mean_spearman_sm:.2f}, {mean_pearson_sm:.2f}')
axs[2].scatter(1,1, s=0,alpha=0.5, color = cmap[2], label = f'L1, Spear., Pear.\n{mean_l1_mm/100:.2f}, {mean_spearman_mm:.2f}, {mean_pearson_mm:.2f}')
#axs[1].scatter(1,1, marker='o', s=5,alpha=0.5, color = cmap[4], label = 'sylph-mot\n
#axs[2].scatter(1,1, marker='x', s=5,alpha=0.5, color = cmap[6], label = 'met-mot\nSpear. ' + f'{mean_spearman_mm:.2f}')
#plt.annotate('sylph\nSpearman: ' + f'{mean_spearman_sy:.2f}', xy=(1,1), xytext=(1,1), textcoords='offset points', color=cmap[0])

for ax in axs:
    ax.legend(frameon=False,  fontsize=6.5, markerscale=0, loc='lower right')

#if genus:
    #plt.suptitle('Genus',fontsize=7)
#else:
#    plt.suptitle('Species',fontsize=7)

plt.tight_layout()
if genus:
    plt.savefig('./hadza_figs/genus_comparison.png', dpi=300)
else:
    plt.savefig('./hadza_figs/species_comparison.png', dpi=300)
plt.show()


#every third
#plt.figure(figsize=(10, 6))

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
