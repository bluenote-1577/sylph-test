import pandas as pd

# Read the file into a DataFrame
df = pd.read_csv('./output_file.tsv', delimiter='\t')

# Apply filters
# Keep rows with ANI > 95 AND Align_fraction_ref > 80 AND Align_fraction_query > 80
condition1 = (df['ANI'] > 95) & (df['Align_fraction_ref'] > 80) & (df['Align_fraction_query'] > 80)
# Keep rows with ANI < 90
condition2 = (df['ANI'] < 90)

# Filter the DataFrame based on each condition
df_filtered_condition1 = df[condition1]
df_filtered_condition2 = df[condition2]

# Write the filtered DataFrame to separate files
df_filtered_condition1.to_csv('output_file_condition1.tsv', sep='\t', index=False)
df_filtered_condition2.to_csv('output_file_condition2.tsv', sep='\t', index=False)

print(f"Filtered data based on condition 1 has been written to 'output_file_condition1.tsv'")
print(f"Filtered data based on condition 2 has been written to 'output_file_condition2.tsv'")
