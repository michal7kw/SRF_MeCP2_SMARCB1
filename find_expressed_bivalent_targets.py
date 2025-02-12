import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Find intersection between bivalent and highly expressed target genes')
parser.add_argument('--target_file', type=str, required=True,
                    help='Path to the target genes file')

args = parser.parse_args()

data_dir = "Gene_lists/bivalent"

# Read the files
bivalent_genes = pd.read_csv(f"{data_dir}/bivalent_NPCs.csv", header=None, names=['Gene'])
high_expr_targets = pd.read_csv(args.target_file, header=None, names=['Gene'])

# Find genes that are in both lists (bivalent)
expressed_bivalent = bivalent_genes[bivalent_genes['Gene'].isin(high_expr_targets['Gene'])]

# Find genes that are only in target list (non-bivalent)
non_bivalent = high_expr_targets[~high_expr_targets['Gene'].isin(bivalent_genes['Gene'])]

# Create output filenames based on input filename
base_name = args.target_file.split('/')[-1].split('.')[0].split('_')[3]
bivalent_output = f"expressed_targeted_bivalent_NPCs_{base_name}.csv"
non_bivalent_output = f"expressed_targeted_non_bivalent_NPCs_{base_name}.csv"

# Save the results
expressed_bivalent.to_csv(f"{data_dir}/{bivalent_output}", index=False, header=False)
non_bivalent.to_csv(f"{data_dir}/{non_bivalent_output}", index=False, header=False)

# Print statistics
print(f"Using target file: {args.target_file}")
print(f"Number of bivalent genes: {len(bivalent_genes)}")
print(f"Number of highly expressed target genes: {len(high_expr_targets)}")
print(f"Number of bivalent genes in intersection: {len(expressed_bivalent)}")
print(f"Number of non-bivalent target genes: {len(non_bivalent)}") 