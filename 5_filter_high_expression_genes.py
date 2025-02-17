import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Filter genes based on baseMean threshold')
parser.add_argument('--threshold', type=float, default=100.0,
                    help='baseMean threshold value (default: 100.0)')
args = parser.parse_args()

"""
This script filters gene lists based on expression level (baseMean) threshold.

Input files:
- ./Gene_lists/DEA_NSC.csv: Differential expression analysis results containing gene names and baseMean values
- ./Gene_lists/targets/all_mecp2_targets_1.csv: First list of target genes
- ./Gene_lists/targets/all_mecp2_targets_2.csv: Second list of target genes

Output files:
- ./Gene_lists/targets/high_expression_targets1_{threshold}.csv: Filtered first target list containing only highly expressed genes
- ./Gene_lists/targets/high_expression_targets2_{threshold}.csv: Filtered second target list containing only highly expressed genes
- ./Gene_lists/targets/high_expression_no_targets_{threshold}.csv: List of highly expressed genes that are not targets
- ./Gene_lists/targets/all_targets1.csv: Complete unfiltered first target list
- ./Gene_lists/targets/all_targets2.csv: Complete unfiltered second target list
- ./Gene_lists/targets/all_no_targets.csv: Complete list of all non-target genes

P.S
`all_mecp2_targets_1.csv` --> joined: 
- `enriched_down_regulated.csv`
- `enriched_not_disregulated.csv`
- `enriched_up_regulated.csv`
(extracted from .xlsx files from email)

`all_mecp2_targets_2.csv` --> subseted (for: both and exon_only): 
- `complete_peak_annotation.csv`
(originated from  `SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_2_rep_in_peaks/`)
"""

# Paths
DEA_PATH = './Gene_lists'
TARGETS_PATH = './Gene_lists/targets'
OUTPUT_PATH = './Gene_lists/targets'

# Read the files
dea = pd.read_csv(f'{DEA_PATH}/DEA_NSC.csv')
targets1 = pd.read_csv(f'{TARGETS_PATH}/all_mecp2_targets_1.csv', header=None, names=['Gene'])
targets2 = pd.read_csv(f'{TARGETS_PATH}/all_mecp2_targets_2.csv', header=None, names=['Gene'])

# Get all genes from DEA data
all_genes = set(dea['gene'])

# Create set of all targets
# all_targets = set(targets1['Gene']).union(set(targets2['Gene']))
all_targets = set(targets2['Gene'])

# Find all genes that are not targets (without expression filtering)
all_no_targets = all_genes - all_targets

# Save complete unfiltered lists
targets1.to_csv(f'{OUTPUT_PATH}/all_targets1.csv', index=False, header=False)
targets2.to_csv(f'{OUTPUT_PATH}/all_targets2.csv', index=False, header=False)
pd.DataFrame(list(all_no_targets), columns=['Gene']).to_csv(
    f'{OUTPUT_PATH}/all_no_targets_final.csv', index=False, header=False)

pd.DataFrame(list(all_targets), columns=['Gene']).to_csv(
    f'{OUTPUT_PATH}/all_targets_final.csv', index=False, header=False)

# Filter DEA data for genes with baseMean > threshold
high_expression_genes = set(dea[dea['baseMean'] > args.threshold]['gene'])

# Find highly expressed genes that are not targets
high_expression_no_targets = high_expression_genes - all_targets

# Filter both target lists to keep only genes with high expression
filtered_targets1 = targets1[targets1['Gene'].isin(high_expression_genes)]
filtered_targets2 = targets2[targets2['Gene'].isin(high_expression_genes)]

# Save filtered lists to new files
filtered_targets1.to_csv(f'{OUTPUT_PATH}/high_expression_targets1_{args.threshold}.csv', index=False, header=False)
filtered_targets2.to_csv(f'{OUTPUT_PATH}/high_expression_targets2_{args.threshold}.csv', index=False, header=False)
pd.DataFrame(list(high_expression_no_targets), columns=['Gene']).to_csv(
    f'{OUTPUT_PATH}/high_expression_no_targets_{args.threshold}.csv', index=False, header=False)

# Print some statistics
print(f"\nUnfiltered statistics:")
print(f"Total number of targets in list 1: {len(targets1)}")
print(f"Total number of targets in list 2: {len(targets2)}")
print(f"Total number of non-target genes: {len(all_no_targets)}")

print(f"\nFiltered statistics (baseMean threshold: {args.threshold}):")
print(f"Number of targets in list 1: {len(filtered_targets1)}")
print(f"Number of targets in list 2: {len(filtered_targets2)}")
print(f"Number of highly expressed non-target genes: {len(high_expression_no_targets)}")