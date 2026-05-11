#!/usr/bin/env python3

#' Compare SMARCB1 binding profiles between up-regulated, down-regulated, and non-regulated genes
#' 
#' Input files:
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/targets/all_targets_final_up_regulated.csv: Up-regulated genes
#' - Gene_lists/targets/all_targets_final_down_regulated.csv: Down-regulated genes
#' - Gene_lists/targets/all_targets_final_not_regulated.csv: Non-regulated genes
#'
#' Output files:
#' - results/metaprofiles_comparison_R/regulated_genes_heatmaps_bm.pdf: 
#'   Heatmap comparing SMARCB1 binding at differently regulated genes

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from typing import List, Tuple, Dict

# Suppress warnings
warnings.filterwarnings('ignore')

def load_gene_list(file_path: str) -> List[str]:
    """Load gene list from CSV file."""
    return pd.read_csv(file_path, header=None)[0].tolist()

def get_tss_regions(gtf_file: str, gene_list: List[str], upstream: int = 2500, downstream: int = 2500) -> pd.DataFrame:
    """Extract TSS regions for given genes."""
    # Read GTF file using pandas directly
    df = pd.read_csv(gtf_file, sep='\t', comment='#',
                    names=['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'attributes'])
    
    # Filter for genes
    df = df[df['feature'] == 'gene']
    
    # Extract gene_name from attributes
    def extract_gene_name(attr):
        for field in attr.split(';'):
            if 'gene_name' in field:
                return field.strip().split('"')[1]
        return None
    
    df['gene_name'] = df['attributes'].apply(extract_gene_name)
    df = df[df['gene_name'].isin(gene_list)]
    
    # Create TSS regions
    tss_regions = []
    for _, row in df.iterrows():
        if row['strand'] == '+':
            start = row['start'] - upstream
            end = row['start'] + downstream
        else:
            start = row['end'] - downstream
            end = row['end'] + upstream
            
        tss_regions.append({
            'chrom': row['seqname'],
            'start': start,
            'end': end,
            'gene_name': row['gene_name'],
            'strand': row['strand']
        })
    
    return pd.DataFrame(tss_regions)

def extract_signal(bw_file: str, regions: pd.DataFrame, bins: int = 100) -> np.ndarray:
    """Extract signal from bigWig file for given regions."""
    bw = pyBigWig.open(bw_file)
    matrix = np.zeros((len(regions), bins))
    
    for i, row in regions.iterrows():
        try:
            # Get signal values
            values = bw.values(row['chrom'], row['start'], row['end'])
            values = [v if v is not None else 0 for v in values]
            
            # Bin the values
            if len(values) > 0:
                binned = np.array_split(values, bins)
                matrix[i] = [np.mean(bin) for bin in binned]
                
            # Reverse the signal if on negative strand
            if row['strand'] == '-':
                matrix[i] = matrix[i][::-1]
                
        except Exception as e:
            print(f"Warning: Error processing region {row['gene_name']}: {str(e)}")
            continue
            
    bw.close()
    return matrix

def plot_heatmaps(up_matrix: np.ndarray, down_matrix: np.ndarray, 
                  not_matrix: np.ndarray, output_path: str):
    """Create and save comparative heatmaps for three categories."""
    # Set up the figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('SMARCB1 Signal at TSS of Regulated Genes', fontsize=16)
    
    # Function to create heatmap
    def create_heatmap(data, ax, title, vmin=None, vmax=None, cmap='YlOrRd'):
        sns.heatmap(data, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
                   xticklabels=False, yticklabels=False)
        ax.set_title(f"{title}\n(n={len(data)})")
        ax.set_xlabel('Distance from TSS')
        
        # Add custom x-axis labels at proper positions
        num_bins = data.shape[1]
        ax.set_xticks([0, num_bins//2, num_bins])
        ax.set_xticklabels(['-2.5kb', 'TSS', '+2.5kb'])
    
    # Find global max for consistent scale
    vmax = 18.0  # Using same scale as original script
    
    # Plot signal for all three groups with same scale
    create_heatmap(up_matrix, ax1, 'Up-regulated Genes', vmax=vmax)
    create_heatmap(down_matrix, ax2, 'Down-regulated Genes', vmax=vmax)
    create_heatmap(not_matrix, ax3, 'Non-regulated Genes', vmax=vmax)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Set paths
    bm_file = "results/bigwig/BM3_CPM.bw"
    gtf_file = "data/gencode.vM10.annotation.gtf.gz"
    up_genes = "Gene_lists/targets/all_targets_final_up_regulated.csv"
    down_genes = "Gene_lists/targets/all_targets_final_down_regulated.csv"
    not_regulated_genes = "Gene_lists/targets/all_targets_final_not_regulated.csv"
    
    # Load gene lists
    print("Loading gene lists...")
    up_list = load_gene_list(up_genes)
    down_list = load_gene_list(down_genes)
    not_list = load_gene_list(not_regulated_genes)
    
    # Get TSS regions
    print("Extracting TSS regions...")
    up_regions = get_tss_regions(gtf_file, up_list)
    down_regions = get_tss_regions(gtf_file, down_list)
    not_regions = get_tss_regions(gtf_file, not_list)
    
    # Process bigWig file for each category
    print("Processing bigWig files...")
    up_matrix = extract_signal(bm_file, up_regions)
    down_matrix = extract_signal(bm_file, down_regions)
    not_matrix = extract_signal(bm_file, not_regions)
    
    # Create heatmaps
    print("Generating heatmaps...")
    output_path = "results/metaprofiles_comparison_R/regulated_genes_heatmaps_bm.pdf"
    plot_heatmaps(up_matrix, down_matrix, not_matrix, output_path)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main() 