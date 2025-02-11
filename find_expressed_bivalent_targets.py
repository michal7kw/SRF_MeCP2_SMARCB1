import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Find intersection between bivalent and highly expressed target genes')
parser.add_argument('--target_file', type=str, required=True,
                    help='Path to the target genes file')

args = parser.parse_args()

# Read the files
bivalent_genes = pd.read_csv('Gene_lists/bivalent/bivalent_NPCs.csv', header=None, names=['Gene'])
high_expr_targets = pd.read_csv(args.target_file, header=None, names=['Gene'])

# Find genes that are in both lists
expressed_bivalent = bivalent_genes[bivalent_genes['Gene'].isin(high_expr_targets['Gene'])]

# Create output filename based on input filename
output_filename = f"expressed_targeted_bivalent_NPCs_{args.target_file.split('/')[-1].split('.')[0].split('_')[3]}.csv"

# Save the result
expressed_bivalent.to_csv(output_filename, index=False, header=False)

# Print statistics
print(f"Using target file: {args.target_file}")
print(f"Number of bivalent genes: {len(bivalent_genes)}")
print(f"Number of highly expressed target genes: {len(high_expr_targets)}")
print(f"Number of genes in intersection: {len(expressed_bivalent)}") 