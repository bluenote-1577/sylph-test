import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
mute = sns.color_palette("muted")
#cmap = sns.color_palette("hls",3)
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


def read_and_filter(file_path, species=True, quantile=1.00, subsamp = False):
    if not subsamp:
        if 'subsamp' in file_path:
            return pd.DataFrame()
    # Reading the file
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df.fillna('NA', inplace=True)
    rl_name = 'Relative Abundance'
    if 'motus' in file_path:
        rl_name = 'Relative Abundance Motus'
        df['Algorithm'] = 'mOTUs'
    elif 'sylph' in file_path:
        rl_name = 'Relative Abundance Sylph'
        df['Algorithm'] = 'sylph'
        print(df)
    elif 'metaphlan' in file_path:
        rl_name = 'Relative Abundance MetaPhlAn4'
        df['Algorithm'] = 'MetaPhlAn4'
    df.rename(columns={df.columns[0]: 'Species', df.columns[1]: rl_name}, inplace=True)
    if 'sylph' in file_path and species:
        df = df[['Species', rl_name, 'Effective coverage', 'Algorithm']]
    else:
        df = df[['Species', rl_name, 'Algorithm']]
    if species:
        species_df = df[df[df.columns[0]].str.contains('s__')]
    else:
        if 'motus' in file_path:
            species_df = df[df[df.columns[0]].str.contains('g__')]
            #split the column by delim ;, and then do not take the last part
            #DO NOT take hte last part, combine the rest
            species_df['Species'] = species_df['Species'].str.split(';').str[:-1].apply(';'.join)
        else:
            species_df = df[df[df.columns[0]].str.contains('g__') & ~df[df.columns[0]].str.contains('s__')]
    if 'sylphmpa' in file_path:
        species_df = species_df[species_df['Species'].str.contains('t__')]
        species_df['Species'] = species_df['Species'].str.split('|').str[:-1].apply(';'.join)
        species_df['Species'] = species_df['Species'].str.replace('|', ';', regex=False)
        print(species_df)
        if species:
            species_df['Effective coverage'] = species_df['Effective coverage'].astype(float)
    if 'motus' in file_path:
        species_df[rl_name] *= 100

    sorted_df = species_df.sort_values(by=rl_name, ascending=False)  # Sort the DataFrame in descending order
    # Calculate the number of rows to select for the top quantile%
    num_rows = int(len(sorted_df) * quantile)
    # Select the top 95% of entries
    top_quantile_percent_df = sorted_df.head(num_rows)
    top_quantile_percent_df.set_index('Species', inplace=True)
    top_quantile_percent_df['Dataset'] = file_path.split('/')[-1].split('_')[0]
    return top_quantile_percent_df

metaphlan_files = glob('./metaphlan_results/*gtdb')
motus_files = glob('./motus_results/*.motus')
sylph_files = glob('./sylph_results/*.sylphmpa')
metaphlan_species = pd.concat([read_and_filter(file) for file in metaphlan_files]).reset_index()
motus_species = pd.concat([read_and_filter(file) for file in motus_files]).reset_index()
sylph_species = pd.concat([read_and_filter(file) for file in sylph_files]).reset_index()

# Merge the DataFrames on 'index' and 'Dataset' columns
intersection_df = metaphlan_species.merge(motus_species, on=['Species', 'Dataset'], how='inner')
intersection_df = intersection_df.merge(sylph_species, on=['Species', 'Dataset'], how='inner')
#remove duplicates
intersection_df = intersection_df.drop_duplicates(subset=['Species', 'Dataset'])
intersection_df = intersection_df.set_index('Species')
print(intersection_df)
intersection_df.to_csv('intersection_df.csv')
# Read and process the first two files
#file1_species = read_and_filter('ERR7745346_metaphlan_gtdb.txt')
#file2_species = read_and_filter('ERR7745346_1.fastq.gz.sylphmpa')

# Find species that are shared between the first two files
#shared_species = file1_species.index.intersection(file2_species.index)
#print(shared_species)

#Find shared species for each 'Dataset' column between all three methods
# only find interesection grouped by datasets
#shared_species = metaphlan_species

# Read the third file (subsampled dataset)
subsampled_sylph_files = glob('./sylph_results/*subsamp*.sylphmpa')
subsampled_motus_files = glob('./motus_results/*subsamp*.motus')
subsampled_metaphlan_files = glob('./metaphlan_results/*subsamp*.gtdb')

subsampled_sylph_df = pd.concat([read_and_filter(file, subsamp=True) for file in subsampled_sylph_files])
subsampled_motus_df = pd.concat([read_and_filter(file, subsamp=True) for file in subsampled_motus_files])
subsampled_metaphlan_df = pd.concat([read_and_filter(file, subsamp=True) for file in subsampled_metaphlan_files])
#ss_all_df = pd.concat([subsampled_sylph_df, subsampled_motus_df, subsampled_metaphlan_df])
subsampled_sylph_df.to_csv('subsampled_sylph_df.csv')
all_ss = [subsampled_sylph_df, subsampled_motus_df, subsampled_metaphlan_df]
for x in all_ss:
    x.reset_index(inplace=True)
# Create a list to store the results
results_data = []

# Populate the list with results
intersection_df = intersection_df.reset_index()
for _, row in intersection_df.iterrows():
    species = row['Species']
    is_present_sylph = ((subsampled_sylph_df['Species'] == row['Species']) &
                  (subsampled_sylph_df['Dataset'] == row['Dataset'])).any()
    is_present_motus = ((subsampled_motus_df['Species'] == row['Species']) &
                  (subsampled_motus_df['Dataset'] == row['Dataset'])).any()
    is_present_metaphlan = ((subsampled_metaphlan_df['Species'] == row['Species']) &
                    (subsampled_metaphlan_df['Dataset'] == row['Dataset'])).any()

    results_data.append({
        'Species': row['Species'], 
        'detect': is_present_sylph,
        'Estimated coverage': row['Effective coverage'],
        'Dataset': row['Dataset'],
        'Algorithm': 'sylph'
    })

    results_data.append({
        'Species': row['Species'], 
        'detect': is_present_motus,
        'Estimated coverage': row['Effective coverage'],
        'Dataset': row['Dataset'],
        'Algorithm': 'mOTUs'
    })

    results_data.append({
        'Species': row['Species'], 
        'detect': is_present_metaphlan,
        'Estimated coverage': row['Effective coverage'],
        'Dataset': row['Dataset'],
        'Algorithm': 'MetaPhlAn4'
        })


# Convert the list of dictionaries to a DataFrame
detect_df = pd.DataFrame(results_data)
detect_df.to_csv('results.tsv',sep='\t')

all_subsamp_df = pd.concat([subsampled_sylph_df, subsampled_motus_df, subsampled_metaphlan_df])
print(all_subsamp_df.head())
results_fp = []

for _, row in all_subsamp_df.iterrows():
    species = row['Species']
    if 'sylph' in row['Algorithm']:
        inter_df = sylph_species
        ab = 'Relative Abundance Sylph'
    elif 'mOTUs' in row['Algorithm']:
        inter_df = motus_species
        ab = 'Relative Abundance Motus'
    elif 'MetaPhlAn4' in row['Algorithm']:
        inter_df = metaphlan_species
        ab = 'Relative Abundance MetaPhlAn4'
    else:
        print(row)
        exit()

    is_present_inter = ((inter_df['Species'] == row['Species']) &
                    (inter_df['Dataset'] == row['Dataset'])).any()

    results_fp.append({
        'Species': row['Species'], 
        'False positive': not is_present_inter,
        'Estimated coverage': row[ab],
        'Dataset': row['Dataset'],
        'Algorithm': row['Algorithm']
    })

fp_df = pd.DataFrame(results_fp)
#get false posities for each algorithm

import numpy as np

for i in [0]:
    if i == 0:
        results = detect_df
    else:
        results = fp_df

    # Convert 'Original Abundance' to numeric and drop rows where it's not available
    results['Estimated coverage'] = pd.to_numeric(results['Estimated coverage'], errors='coerce')
    if i == 0:
        results['Estimated coverage'] = pd.to_numeric(results['Estimated coverage'], errors='coerce')/10
    #results = results.dropna(subset=['Estimated coverage'])

    # Create log-spaced bins
    min_abundance = results['Estimated coverage'].min()
    max_abundance = 1.0
    bins = np.logspace(np.log10(min_abundance), np.log10(max_abundance), num=15)

    # Bin data
    results['Abundance Bin'] = pd.cut(results['Estimated coverage'], bins, include_lowest=True).apply(lambda x: x.mid)

    #turn true to 1 and false to 0 in the above dataframe
    if i == 0:
        results['detect'] = results['detect'].astype(int)
    else:
        results['False positive'] = results['False positive'].astype(int)
    sylph_res = results[results['Algorithm'] == 'sylph']
    motus_res = results[results['Algorithm'] == 'mOTUs']
    metaphlan_res = results[results['Algorithm'] == 'MetaPhlAn4']
    cmap = [mute[0], mute[4], mute[6]]

    # Plotting
    plt.figure(figsize=(6.5/2.54, 6/2.54))
    if i == 0:
        y = 'detect'
    else:
        y = 'False positive'
    sns.lineplot(data=sylph_res, x="Abundance Bin", y=y, marker='o', label='sylph',color=cmap[0])
    sns.lineplot(data=motus_res, x="Abundance Bin", y=y, marker='o', label = 'mOTUs3', color=cmap[2], linestyle='--')
    sns.lineplot(data=metaphlan_res, x="Abundance Bin", y=y, marker='o', label='MetaPhlAn4', color=cmap[1], linestyle='-.')
    plt.xscale('log')
    plt.xlabel('Estimated effective coverage')
    plt.ylabel('Detection probability')
    plt.title('Species detection probability under subsampling', fontsize=7)
    #plt.xticks(rotation=45)
    plt.gca().spines[['top','right']].set_visible(False)
    plt.tight_layout()


    plt.legend(frameon=False)
    if i == 0:
        plt.savefig('./hadza_figs/detection_probability.svg')
    else:
        plt.savefig('./hadza_figs/false_positive.svg')
    plt.show()


