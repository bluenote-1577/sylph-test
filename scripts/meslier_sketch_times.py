import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt


cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("deep")
plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


def parse_time_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        user_time = float(lines[1].split(":")[1].strip().split(" ")[0])
        system_time = float(lines[2].split(":")[1].strip().split(" ")[0])
        wall_time = lines[4].split(" ")[-1].strip()
        #split by :, reverse the iteration and add by 60^i
        wall_time = sum(float(a) * 60 ** i for i, a in enumerate(reversed(wall_time.split(':'))))
        memory_usage = float(lines[9].split(":")[1].strip().split(" ")[0])/1_000_000
        cpu_time = user_time + system_time
        return cpu_time, wall_time, memory_usage

# List of files
files = [['./meslier_times/mash_screen_ill.time','./meslier_times/sourmash_search_ill.time', './meslier_times/sourmash_sketch_ill.time', './meslier_times/sylph_sketch_ill.time', './meslier_times/sylph_profile_ill.time'],['./meslier_times/mash_screen_nano.time','./meslier_times/sourmash_search_nano.time', './meslier_times/sourmash_sketch_nano.time', './meslier_times/sylph_sketch_nano.time', './meslier_times/sylph_profile_nano.time'],['./meslier_times/mash_screen_pac.time','./meslier_times/sourmash_search_pac.time', './meslier_times/sourmash_sketch_pac.time', './meslier_times/sylph_sketch_pac.time', './meslier_times/sylph_profile_pac.time']]
titles = ['Illumina', 'Nanopore (old)', 'PacBio']
labels = ['sylph query', 'Mash screen', 'Sourmash search']
re = [1,0,2]

# Parsing data from files
dfs = []
for j in range(3):
    data = []
    for i,file in enumerate(files[j]):
        user_time, wall_time, memory_usage = parse_time_file(file)
        data.append({
            'File': 'Mash screen',
            'CPU Time (s)': user_time,
            'Wall Time': wall_time,
            'Memory Usage (GB)': memory_usage,
            'Technology': titles[j]
        })

    # Creating DataFrame
    data.append({'File': 'sylph query', 'CPU Time (s)': data[4]['CPU Time (s)'] + data[3]['CPU Time (s)'], 'Wall Time': data[4]['Wall Time'] + data[3]['Wall Time'], 'Memory Usage (GB)': np.max([data[4]['Memory Usage (GB)'], data[3]['Memory Usage (GB)']]), 'Technology': titles[j]})
    #delete sylph profile and sylph sketch

    # Creating DataFrame
    data.append({'File': 'Sourmash search', 'CPU Time (s)': data[2]['CPU Time (s)'] + data[1]['CPU Time (s)'], 'Wall Time': data[2]['Wall Time'] + data[1]['Wall Time'], 'Memory Usage (GB)': np.max([data[2]['Memory Usage (GB)'], data[1]['Memory Usage (GB)']]), 'Technology': titles[j]})
    #delete sylph profile and sylph sketch
    df = pd.DataFrame(data)
    df = df.drop([1,2,3,4])
    print(df)
    dfs.append(df)

df = pd.concat(dfs)
print(df)

# Convert wall time to minutes for easier interpretation
df['Wall Time (s)'] = df['Wall Time']

#df['Wall Time (s)'] = df['Wall Time (s)']/50
# Set up the matplotlib figure
#9 plots, 3 rows, 3 columns, wall, cpu, memory
fig, ax = plt.subplots(3, 3, figsize=(14*cm, 14*cm))

cmap = [cmap[0],cmap[1],cmap[2]]
for i , title in enumerate(titles):
    # Wall time
    sns.barplot(x='File', y='Wall Time (s)', data=df[df['Technology'] == title], palette=cmap, ax=ax[0,i], order=labels)
    ax[0,i].set_title(title)
    ax[0,i].set_xticklabels(ax[0,i].get_xticklabels(), rotation=45)
    ax[0,i].set_xlabel("")
    ax[0,i].set_ylabel("Average Wall Time (s)")
    ax[0,i].spines['top'].set_visible(False)
    ax[0,i].spines['right'].set_visible(False)
    # annotate bar with actual value
    for index, value in enumerate(df[df['Technology'] == title]['Wall Time (s)']):
        ax[0,i].text(re[index], value, str(round(value)), ha='center', va='bottom')
    # CPU time
    sns.barplot(x='File', y='CPU Time (s)', data=df[df['Technology'] == title], palette=cmap, ax=ax[1,i], order=labels)
    #ax[1,i].set_title(title)
    ax[1,i].set_xticklabels(ax[1,i].get_xticklabels(), rotation=45)
    ax[1,i].set_xlabel("")
    ax[1,i].set_ylabel("Average CPU Time (s)")
    ax[1,i].spines['top'].set_visible(False)
    ax[1,i].spines['right'].set_visible(False)
    for index, value in enumerate(df[df['Technology'] == title]['CPU Time (s)']):
        ax[1,i].text(re[index], value, str(round(value)), ha='center', va='bottom')

    # annotate bar with actual value
    #for index, value in enumerate(df['CPU Time (s)']):
    #    ax[i,1].text(index, value, str(round(value)), ha='center', va='bottom')
    # Memory usage
    sns.barplot(x='File', y='Memory Usage (GB)', data=df[df['Technology'] == title], palette=cmap, ax=ax[2,i], order=labels)
    #ax[2,i].set_title(title)
    ax[2,i].set_xticklabels(ax[2,i].get_xticklabels(), rotation=45)
    ax[2,i].set_xlabel("")
    ax[2,i].set_ylabel("Average Memory Usage (GB)")
    ax[2,i].spines['top'].set_visible(False)
    ax[2,i].spines['right'].set_visible(False)
    # annotate bar with actual value
    for index, value in enumerate(df[df['Technology'] == title]['Memory Usage (GB)']):
        ax[2,i].text(re[index], value, str(round(value,1)), ha='center', va='bottom')

    #for index, value in enumerate(df['Memory Usage (GB)']):
    #    ax[i,2].text(index, value, str(round(value)), ha='center', va='bottom')


# Show plot
plt.tight_layout()
plt.savefig('supp_figs/meslier_sketch_times.svg')
plt.show()
