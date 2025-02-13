import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Filter genes based on baseMean threshold')
parser.add_argument('--threshold', type=float, default=100.0,
                    help='baseMean threshold value (default: 100.0)')
args = parser.parse_args()

# Read the files
dea = pd.read_csv('./Gene_lists/DEA_NSC.csv')
targets1 = pd.read_csv('./Gene_lists/targets/all_mecp2_targets_1.csv', header=None, names=['Gene'])
targets2 = pd.read_csv('./Gene_lists/targets/all_mecp2_targets_2.csv', header=None, names=['Gene'])

# Filter DEA data for genes with baseMean > threshold
high_expression_genes = set(dea[dea['baseMean'] > args.threshold]['gene'])

# Create set of all targets
all_targets = set(targets1['Gene']).union(set(targets2['Gene']))

# Find highly expressed genes that are not targets
high_expression_no_targets = high_expression_genes - all_targets

# Filter both target lists to keep only genes with high expression
filtered_targets1 = targets1[targets1['Gene'].isin(high_expression_genes)]
filtered_targets2 = targets2[targets2['Gene'].isin(high_expression_genes)]

# Save filtered lists to new files
filtered_targets1.to_csv(f'high_expression_targets1_{args.threshold}.csv', index=False, header=False)
filtered_targets2.to_csv(f'high_expression_targets2_{args.threshold}.csv', index=False, header=False)

# Save highly expressed no-targets to a new file
pd.DataFrame(list(high_expression_no_targets), columns=['Gene']).to_csv(
    f'./Gene_lists/targets/high_expression_no_targets_{args.threshold}.csv', index=False, header=False)

# Print some statistics
print(f"Using baseMean threshold: {args.threshold}")
print(f"Original number of targets in list 1: {len(targets1)}")
print(f"Filtered number of targets in list 1: {len(filtered_targets1)}")
print(f"Original number of targets in list 2: {len(targets2)}")
print(f"Filtered number of targets in list 2: {len(filtered_targets2)}")
print(f"Number of highly expressed non-target genes: {len(high_expression_no_targets)}") 