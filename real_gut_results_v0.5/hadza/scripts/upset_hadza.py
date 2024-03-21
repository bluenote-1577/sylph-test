import pandas as pd
from glob import glob
import seaborn as sns
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import UpSet
from upsetplot import plot as upplot

quantile_share = True
venn = False
upset = False
quantiles = [0.10, 0.2, 0.30, 0.4, 0.5,0.6,0.7,0.8,0.9, 1]
#quantiles = [1,0.25]
cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


# Replace these with the paths to your files
#file_path2 = './ERR7745346_metaphlan_gtdb.txt'
#file_path2 = './ERR7745346_1.fastq.gz.sylphmpa'
#file_path2 = './ERR7745346_1.fastq.gz.sylphmpa'

## subsampled
#file_path1 = './motus.csv'
#file_path1 = './motus_noambig.csv'
#file_path2 = './ERR7745346_subsamp0.10_1.fastq.sylphmpa'
#file_path2 = './ERR7745346_subsamp_metaphlan_gtdb.txt'

metaphlan_files = glob('metaphlan_results/*.gtdb')
metaphlan_files = [x for x in metaphlan_files if 'subsamp' not in x]
sylph_files = glob('sylph_results/*.sylphmpa')
sylph_files = [x for x in sylph_files if 'subsamp' not in x]
motus_files = glob('motus_results/*.motus')
motus_files = [x for x in motus_files if 'subsamp' not in x]
#sort
metaphlan_files.sort()
sylph_files.sort()
motus_files.sort()
#quantile = 0.25


#file_path1 = './ERR7803603_motus.txt.motus.noambig'
#file_path2 = './ERR7803603_1.fastq.gz.sylphmpa'
#file_paths = [file_path1, file_path2, file_path3]


def read_and_filter(file_path, species, quantile=1):
    # Reading the file
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df.fillna('NA', inplace=True)
    df.rename(columns={df.columns[0]: 'Species', df.columns[1]: 'Relative Abundance'}, inplace=True)
    #if 'sylphmpa' in file_path:
    #    df.rename(columns={df.columns[0]: 'Species', df.columns[2]: 'Relative Abundance'}, inplace=True)
    #else:
    #    df.rename(columns={df.columns[0]: 'Species', df.columns[1]: 'Relative Abundance'}, inplace=True)
    if species:
        species_df = df[df[df.columns[0]].str.contains('s__') & ~df[df.columns[0]].str.contains('t__')]
    else:
        if 'motus' in file_path:
            species_df = df[df[df.columns[0]].str.contains('g__')]
            #split the column by delim ;, and then do not take the last part
            #DO NOT take hte last part, combine the rest
            species_df['Species'] = species_df['Species'].str.split(';').str[:-1].apply(';'.join)
        else:
            species_df = df[df[df.columns[0]].str.contains('g__') & ~df[df.columns[0]].str.contains('s__')]
    if 'sylphmpa' in file_path:
        species_df['Species'] = species_df['Species'].str.replace('|', ';', regex=False)
    if 'motus' in file_path:
        species_df['Relative Abundance'] *= 100

    sorted_df = species_df.sort_values(by='Relative Abundance', ascending=False)  # Sort the DataFrame in descending order
    print(file_path, sorted_df.head())
    # Calculate the number of rows to select for the top quantile%
    num_rows = int(len(sorted_df) * quantile)
    # Select the top 95% of entries
    top_quantile_percent_df = sorted_df.head(num_rows)
    return top_quantile_percent_df



for spec in [True, False]:
    unique_sylph_q = []
    unique_metaphlan_q = []
    unique_motus_q = []
    shared_all = []
    for quantile in quantiles:
        df = pd.DataFrame(index=pd.MultiIndex.from_tuples([], names=['source1', 'source2', 'source3']))
        for file_path1, file_path2, file_path3 in zip(metaphlan_files, sylph_files, motus_files):
            file_paths = [file_path1, file_path2, file_path3]
            # Reading and filtering data
            species_abundances = [read_and_filter(file_path, spec, quantile) for file_path in file_paths]

            # Define log bins

            # Initialize a dictionary to store Jaccard indices
            jaccard_indices = []
            contain1 = []
            contain2 = []

            species_in_files = [set(species_abundance['Species']) for species_abundance in species_abundances]
            #print(species_in_files)
            #for bin in bins:

                #species_in_bin_files = [set(species_abundance[species_abundance['Relative Abundance'] > bin].index) for species_abundance in species_abundances]

                # Calculate Jaccard index
        #        intersection1 = species_in_bin_file1.intersection(species_in_2)
        #        intersection2 = species_in_bin_file2.intersection(species_in_1)
        #        diff1 = species_in_bin_file1.difference(species_in_bin_file2)
        #        diff2 = species_in_bin_file2.difference(species_in_bin_file1)
        #        print('metaphlan unique')
        #        for x in diff1:
        #            print(x)
        #        print('sylph unique')
        #        for x in diff2:
        #            print(x)
        #
                #contain1.append(len(intersection1)/len(species_in_bin_file1))
                #contain2.append(len(intersection2)/len(species_in_bin_file2))

            # Create a DataFrame to hold the membership information
            set1 = species_in_files[0]
            set2 = species_in_files[1]
            set3 = species_in_files[2]
            all_elements = set().union(set1, set2, set3)
            index = pd.MultiIndex.from_tuples(
                [(element in set1, element in set2, element in set3) for element in all_elements],
                names=[file_path1, file_path2, file_path3]
            )

            # Create a DataFrame with counts for each combination
            data = pd.Series(1, index=index).groupby(level=list(range(3))).sum()/len(sylph_files)
            data = data.round()
            #concatenate data to an existing dataframe
            data.name = "count"
            df = pd.concat([df, data], axis=0, ignore_index=False)
            # Assuming 's' is your Series
            df = df.groupby(level=[0, 1, 2]).sum()
            print(df)

                # Check current level names
        # Reset the index to turn the MultiIndex into columns
        df = df.reset_index()

        # Rename the columns
        df.columns = ['MetaPhlAn4', 'sylph', 'mOTUs3', 'value']

        # Convert the DataFrame back to a Series with a name
        df = df.set_index(['MetaPhlAn4', 'sylph', 'mOTUs3'])['value']
        print(df)


        if quantile_share:
            motus_unique = df[(df.index.get_level_values('MetaPhlAn4') == False) & (df.index.get_level_values('sylph') == False) & (df.index.get_level_values('mOTUs3') == True)]
            sylph_unique = df[(df.index.get_level_values('MetaPhlAn4') == False) & (df.index.get_level_values('sylph') == True) & (df.index.get_level_values('mOTUs3') == False)]
            metaphlan_unique = df[(df.index.get_level_values('MetaPhlAn4') == True) & (df.index.get_level_values('sylph') == False) & (df.index.get_level_values('mOTUs3') == False)]

            if len(sylph_unique) == 0:
                unique_sylph_q.append(0)
            else:
                unique_sylph_q.append(float(sylph_unique))
            if len(motus_unique) == 0:
                unique_motus_q.append(0)
            else:
                unique_motus_q.append(motus_unique)
            if len(metaphlan_unique) == 0:
                unique_metaphlan_q.append(0)
            else:
                unique_metaphlan_q.append(metaphlan_unique)
            shared_df = df[(df.index.get_level_values('MetaPhlAn4') == True) & (df.index.get_level_values('sylph') == True) & (df.index.get_level_values('mOTUs3') == True)]
            shared_all.append(shared_df)
            
        if venn:
            # Convert to dictionary for venn3
            venn_dict = {}
            for idx, count in df.items():
                # Convert index tuple to a string of '1's and '0's
                key = ''.join(['1' if x else '0' for x in idx])
                venn_dict[key] = count

            # Plotting Venn diagram
            plt.figure(figsize=(5*cm, 5*cm))
            venn3(subsets=venn_dict, set_labels=('MetaPhlAn4', 'sylph', 'mOTUs3'))
            if spec:
                if quantile == 1:
                    plt.title("Average number of species detected\n(all species)", fontsize=7)
                else:
                    plt.title("Average number of species detected\n(top {}% abundance quantile)".format(quantile*100), fontsize=7)
            else:
                if quantile == 1:
                    plt.title("Average number of genera detected\n(all genera)", fontsize=7)
                else:
                    plt.title("Average number of genera detected\n(top {}% abundance quantile)".format(quantile*100), fontsize=7)
            plt.tight_layout()
            plt.savefig(f'hadza_figs/venn3-quantile{quantile}{spec}.pdf')
            plt.show()


        elif upset:
            # Show the new Series
            fig = plt.figure(figsize=(8*cm, 8*cm))
            try:
                #upplot(df, sort_by='cardinality')
                up = UpSet(df, totals_plot_elements=0, element_size=None)
                up.plot(fig=fig)
                if spec:
                    plt.ylabel('Number of species')
                else:
                    plt.ylabel('Number of genera')
                #upplot(df,fig=fig, element_size=None,total_plot_elements=0)
            except ValueError as e:
                print("Error in plotting:", e)

            # rename the true/false columns to the file names
            df.columns = file_paths
            plt.show()
                # Create the UpSet plot
                # Check the type of index

            
    if quantile_share:
        plt.figure(figsize=(6*cm, 6*cm))
        print(unique_metaphlan_q)
        plt.plot(quantiles, unique_motus_q, 'o-', label='mOTUs3', color=cmap[6], ms=4, linestyle='--')
        plt.plot(quantiles, unique_metaphlan_q, 'o-', label='MetaPhlAn4', color=cmap[4],ms=4, linestyle='-.')
        plt.plot(quantiles, unique_sylph_q, 'o-', label='sylph', color=cmap[0],ms=4)
        plt.plot(quantiles, shared_all, 'o-', label='All shared', color='black', linestyle=':',ms=4)
        plt.xlabel('Top fraction of genomes\n(sorted by abundance)')
        if spec:
            plt.ylabel('Number of species unique to method')
            plt.title("Species")
        else:
            plt.ylabel('Number of genera unique to method')
            plt.title("Genus")
        plt.legend(frameon=False)
        sns.despine()
        plt.tight_layout()
        if spec:
            plt.savefig('hadza_figs/unique-quantile.pdf')
        else:
            plt.savefig('hadza_figs/unique-quantile-genera.pdf')
        plt.show()
