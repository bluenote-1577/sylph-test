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
files = ['times/metaphlan_all.time', 'times/motus_all.time', 'times/sylph_profile.time', 'times/sylph_sketch.time']
labels = ['MetaPhlAn4', 'mOTUs3', 'sylph profile', 'sylph sketch' ]

# Parsing data from files
data = []
for i,file in enumerate(files):
    user_time, wall_time, memory_usage = parse_time_file(file)
    data.append({
        'File': labels[i],
        'CPU Time (s)': user_time/50,
        'Wall Time': wall_time/50,
        'Memory Usage (GB)': memory_usage
    })

# Creating DataFrame
data.append({'File': 'sylph', 'CPU Time (s)': data[2]['CPU Time (s)'] + data[3]['CPU Time (s)'], 'Wall Time': data[2]['Wall Time'] + data[3]['Wall Time'], 'Memory Usage (GB)': np.max([data[2]['Memory Usage (GB)'], data[3]['Memory Usage (GB)']])})
df = pd.DataFrame(data)
#delete sylph profile and sylph sketch
df = df.drop([2,3])


# Convert wall time to minutes for easier interpretation
df['Wall Time (s)'] = df['Wall Time']

#df['Wall Time (s)'] = df['Wall Time (s)']/50
# Set up the matplotlib figure
a = plt.figure(figsize=(10*cm, 5*cm))

cmap = [cmap[4],cmap[6],cmap[0]]
# User Time Plot
plt.subplot(1, 3, 1)
sns.barplot(x='File', y='CPU Time (s)', data=df, palette=cmap)
#plt.title('CPU Time (secs)')
plt.xticks(rotation=45)
plt.xlabel("")
plt.ylabel("Average CPU Time (s)")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
# annotate bar with actual value
for index, value in enumerate(df['CPU Time (s)']):
    plt.text(index, value, str(round(value)), ha='center', va='bottom')

#set spines tonone

# Wall Time Plot
plt.subplot(1, 3, 2)
sns.barplot(x='File', y='Wall Time (s)', data=df,palette=cmap)
#plt.title('Wall Time (secs)')
plt.xticks(rotation=45)
plt.xlabel("")
plt.ylabel("Average Wall Time (s)")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
for index, value in enumerate(df['Wall Time (s)']):
    plt.text(index, value, str(round(value,1)), ha='center', va='bottom')



# Memory Usage Plot
plt.subplot(1, 3, 3)
sns.barplot(x='File', y='Memory Usage (GB)', data=df,palette=cmap)
#plt.title('Memory Usage (KB)')
plt.xticks(rotation=45)
plt.xlabel("")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
for index, value in enumerate(df['Memory Usage (GB)']):
    plt.text(index, value, str(round(value,1)), ha='center', va='bottom')

# Add one title for all subplots
plt.suptitle('Mean time and max memory (50 gut metagenomes, 50 threads)', fontsize=7)

# Show plot
plt.tight_layout()
plt.savefig('supp_figs/timing.svg')
plt.show()
