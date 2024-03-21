import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from scipy import stats
from sklearn.linear_model import HuberRegressor

cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
cmap = sns.color_palette("muted")
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})


# Replace 'your_file.tsv' with the path to your TSV file
file_path = sys.argv[1]

# Reading the TSV file
df = pd.read_csv(file_path, sep='\t')
df = df[df['ANI'] >= 95]
df = df[df['ANI'] != 100]

# Calculating the average of Align_fraction_ref and Align_fraction_query
df['Average_Align_Fraction'] = df[['Align_fraction_ref', 'Align_fraction_query']].mean(axis=1)

for i in range(1,6):
    #x = df[df['ANI'] < 95+i]['Average_Align_Fraction'].mean()
    #x = df[df['ANI'] < 95+i and df['ANI'] >= 95+i-1]['Average_Align_Fraction'].mean()
    #between 95 and 96, 96 and 97, etc. pandas dataframe manip
    x = df[df['ANI'] < 95+i][df['ANI'] >= 95+i-1]['Average_Align_Fraction'].mean()
    print(x)

from sklearn.linear_model import RANSACRegressor

huber = HuberRegressor()
huber.fit(df[['ANI']], df['Average_Align_Fraction'])

ransac = RANSACRegressor(random_state=0, min_samples=25)
ransac.fit(df[['ANI']].values, df[['Average_Align_Fraction']].values)

# Create the line of best fit
#line = huber.predict(df[['ANI']].sort_values(by='ANI'))
line = ransac.predict(df[['ANI']].sort_values(by='ANI'))

#fit a ransac resssor to a dataframe



# Plotting
plt.figure(figsize=(8*cm, 5*cm))
#Change transparency of the points
plt.scatter(df['ANI'], df['Average_Align_Fraction'],color='blue', s=5, alpha=0.5, label='ANI vs Average Alignment Fraction')
plt.plot(df[['ANI']].sort_values(by='ANI'), line, 'r', label='Huber regression line')

#Remove spines,tight layout
sns.despine()
plt.title('GTDB R207 genomes vs R89 genomes')
plt.xlabel('ANI (%)')
plt.ylabel('Average aligned fraction')
plt.xlim(95, 100)
plt.ylim(50, 100)
plt.grid(True)
plt.tight_layout()
plt.savefig('supp_figs/ani_vs_af_r207_vs_r89.svg')
plt.show()
