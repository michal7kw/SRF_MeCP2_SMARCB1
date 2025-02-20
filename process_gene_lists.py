import pandas as pd
import numpy as np

# Define significance thresholds
LOG2FC_THRESHOLD = 1.0  # Log2 fold change threshold
PADJ_THRESHOLD = 0.05   # Adjusted p-value threshold

# Read the input files
def load_and_process_genes():
    # Read target genes (one column, no header)
    target_genes = pd.read_csv('Gene_lists/targets/all_targets_final.csv', header=None)
    target_genes = target_genes[0].tolist()  # Convert to list
    
    # Read DEA results
    dea_results = pd.read_csv('Gene_lists/DEA_NSC.csv')
    
    # Initialize lists for categorized genes
    up_regulated = []
    down_regulated = []
    not_regulated = []
    
    # Categorize each target gene
    for gene in target_genes:
        # Find the gene in DEA results
        gene_data = dea_results[dea_results['gene'] == gene]
        
        if len(gene_data) == 0:
            # Skip genes not found in DEA results
            continue
            
        log2fc = gene_data['log2FoldChange'].values[0]
        padj = gene_data['padj'].values[0]
        
        # Categorize based on thresholds
        if padj <= PADJ_THRESHOLD:
            if log2fc >= LOG2FC_THRESHOLD:
                up_regulated.append(gene)
            elif log2fc <= -LOG2FC_THRESHOLD:
                down_regulated.append(gene)
            else:
                not_regulated.append(gene)
        else:
            not_regulated.append(gene)
    
    return up_regulated, down_regulated, not_regulated

def main():
    up_regulated, down_regulated, not_regulated = load_and_process_genes()
    
    # Print results
    print("Upregulated genes:", up_regulated)
    print("Downregulated genes:", down_regulated)
    print("Not differentially expressed genes:", not_regulated)
    
    # Optionally save results to files
    pd.Series(up_regulated).to_csv('Gene_lists/targets/all_targets_final_up_regulated.csv', index=False)
    pd.Series(down_regulated).to_csv('Gene_lists/targets/all_targets_final_down_regulated.csv', index=False)
    pd.Series(not_regulated).to_csv('Gene_lists/targets/all_targets_final_not_regulated.csv', index=False)

if __name__ == "__main__":
    main() 