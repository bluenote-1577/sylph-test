import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
cm = 1/2.54  # centimeters in inches\n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

# Directories for Sylph, Motus, and MetaPhlAn results
sylph_dir = "hadza/sylph_results"
motus_dir = "hadza/motus_results"
metaphlan_dir = "hadza/metaphlan_results"


#sylph_dir = "sylph_results_round1"
#motus_dir = "motus_results"
#metaphlan_dir = "metaphlan_results"

# Function to count unique genera/species
def count_unique(file_path, skiprows=1, genera=True):
    df = pd.read_csv(file_path, sep='\t', header=None, skiprows=skiprows, comment = '#')
    col = df.columns[0]
    if genera:
        #first col as string
        df[col] = df[col].astype(str)
        df = df[df[col].str.contains('s__')]
        df = df[~(df[col].str.contains('t__'))]
        df[col] = df[col].str.split('s__').str[0]
        return len(df[col].unique())
    else:
        df[col] = df[col].astype(str)
        df = df[df[col].str.contains('s__')]
        df = df[~(df[col].str.contains('t__'))]
        return len(df[col].unique())

# Prepare data for plotting
data = {'SampleID': [], 'Sylph': [], 'Motus': [], 'MetaPhlAn': [], 'Sylph gen': [], 'Motus gen': [], 'MetaPhlAn gen': []}

# Loop through Sylph result files
for filename in os.listdir(sylph_dir):
    if filename.endswith(".sylphmpa") and "_subsamp" not in filename:
        sample_id = filename.split("_")[0]
        print(filename, sample_id)
        motus_file = os.path.join(motus_dir, f"{sample_id}_motus.txt.motus")
        metaphlan_file = os.path.join(metaphlan_dir, f"{sample_id}_metaphlan.txt.gtdb")
        sylph_file = os.path.join(sylph_dir, filename)

        if not os.path.exists(motus_file) or not os.path.exists(metaphlan_file) or not os.path.exists(sylph_file):
            print(f"Skipping {motus_file}")
            continue

        sylph_count = count_unique(os.path.join(sylph_dir, filename))
        motus_count = count_unique(os.path.join(motus_dir, f"{sample_id}_motus.txt.motus"),skiprows=0)
        metaphlan_count = count_unique(os.path.join(metaphlan_dir, f"{sample_id}_metaphlan.txt.gtdb"), skiprows=2)

        sylph_count_gen = count_unique(os.path.join(sylph_dir, filename), genera=False)
        motus_count_gen = count_unique(os.path.join(motus_dir, f"{sample_id}_motus.txt.motus"),skiprows=0, genera=False)
        metaphlan_count_gen = count_unique(os.path.join(metaphlan_dir, f"{sample_id}_metaphlan.txt.gtdb"), skiprows=2, genera=False)

        data['SampleID'].append(sample_id)
        data['Sylph'].append(sylph_count)
        data['Motus'].append(motus_count)
        data['MetaPhlAn'].append(metaphlan_count)


        data['Sylph gen'].append(sylph_count_gen)
        data['Motus gen'].append(motus_count_gen)
        data['MetaPhlAn gen'].append(metaphlan_count_gen)



# Create DataFrame
df = pd.DataFrame(data)
#remove rows with 0 in any of the columns
df = df[(df.T != 0).all()]

# Plotting
fig,axs = plt.subplots(2,1, figsize=(4.5*cm, 8*cm))

mot_c = cmap[6]
meta_c = cmap[4]

# Sylph vs Motus
ms = 20
sns.scatterplot(ax=axs[1],y=df['Motus'], x=df['Sylph'], label = 'mOTUs3', color=mot_c, s=ms)
sns.scatterplot(ax=axs[1],y=df['MetaPhlAn'],x=df['Sylph'], label='MP4', color=meta_c, marker='s', s=ms)

sns.scatterplot(ax=axs[0],y=df['Motus gen'], x=df['Sylph gen'], label = 'mOTUs3', color=mot_c, s=ms)
sns.scatterplot(ax=axs[0],y=df['MetaPhlAn gen'],x=df['Sylph gen'], label='MP4', color=meta_c, marker='s', s=ms)
min_gen = min(df['Sylph'].min(), df['Motus'].min(), df['MetaPhlAn'].min())
min_spec = min(df['Sylph gen'].min(), df['Motus gen'].min(), df['MetaPhlAn gen'].min())
max_gen = max(df['Sylph'].max(), df['Motus'].max(), df['MetaPhlAn'].max())
max_spec = max(df['Sylph gen'].max(), df['Motus gen'].max(), df['MetaPhlAn gen'].max())
axs[1].plot([min_gen, max_gen], [min_gen, max_gen], color='black', linestyle='--',linewidth=1.0)
axs[0].plot([min_spec, max_spec], [min_spec, max_spec], color='black', linestyle='--',linewidth=1.0)
#axs[0].plot([300, 750], [300, 750], color='black', linestyle='--')
#axs[1].plot([180, 350], [180, 350], color='black', linestyle='--')

print(df['Motus gen'].mean(),df['Sylph gen'].mean(), df['MetaPhlAn gen'].mean())
print(df['Motus'].mean(), df['Sylph'].mean(), df['MetaPhlAn'].mean())

sns.despine()
axs[0].set_xlabel('Sylph # species')
axs[0].set_ylabel('# species')

axs[1].set_xlabel('Sylph # genera')
axs[1].set_ylabel('# genera')

axs[0].legend(frameon=False,markerscale=1.0,loc='upper left')
axs[1].legend(frameon=False,markerscale=1.0,loc='upper left')
#grid
#axs[0].grid(True)
#axs[1].grid(True)
plt.tight_layout()
plt.savefig('hadza/hadza_figs/sylph_motus_metaphlan_number.svg')
plt.show()
