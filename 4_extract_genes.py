import pandas as pd

# Read the CSV file
df = pd.read_csv('complete_peak_annotation.csv')

# Filter rows where binding_type is either 'both' or 'exo_only'
filtered_df = df[df['binding_type'].isin(['both', 'exo_only'])]

# Get unique gene symbols, remove NaN values, and sort them alphabetically
gene_symbols = sorted(filtered_df['SYMBOL'].dropna().unique())

# Create a new dataframe with just the gene symbols
genes_df = pd.DataFrame(gene_symbols, columns=['Gene'])

# Save to a new CSV file without index
genes_df.to_csv('all_mecp2_targets_2.csv', index=False, header=False) 