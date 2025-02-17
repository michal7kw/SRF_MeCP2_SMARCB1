#!/usr/bin/env python3

#' Compare SMARCB1 binding profiles between bivalent and non-bivalent genes using heatmaps
#' 
#' Input files:
#' - results/bigwig/BG1_CPM.bw: Control replicate 1 bigWig file
#' - results/bigwig/BG2_CPM.bw: Control replicate 2 bigWig file 
#' - results/bigwig/BG3_CPM.bw: Control replicate 3 bigWig file
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv: List of bivalent genes
#' - Gene_lists/bivalent/expressed_targeted_non_bivalent_NPCs_1000.csv: List of non-bivalent genes
#'
#' Output files:
#' - results/heatmaps_comparison_R/bivalent_nonbivalent_heatmaps.pdf: 
#'   Heatmap comparing SMARCB1 binding at bivalent vs non-bivalent genes
#'
#' The script generates heatmaps comparing SMARCB1 binding profiles around TSS regions
#' between bivalent and non-bivalent genes.

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from typing import List, Tuple, Dict
from multiprocessing import Pool
from functools import partial

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

def process_bigwig_files(bg_files: List[str], bm_file: str, regions: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Process multiple bigWig files and calculate mean signal."""
    # Process BG files in parallel
    with Pool() as pool:
        bg_matrices = pool.map(partial(extract_signal, regions=regions), bg_files)
    bg_matrix = np.mean(bg_matrices, axis=0)
    
    # Process BM file
    bm_matrix = extract_signal(bm_file, regions)
    
    return bg_matrix, bm_matrix

def plot_heatmaps(bg_bivalent: np.ndarray, bm_bivalent: np.ndarray, 
                  bg_nonbivalent: np.ndarray, bm_nonbivalent: np.ndarray,
                  output_path: str):
    """Create and save comparative heatmaps."""
    # Set up the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('SMARCB1 Signal at TSS', fontsize=16)
    
    # Function to create heatmap
    def create_heatmap(data, ax, title, vmin=None, vmax=None, cmap='YlOrRd'):
        sns.heatmap(data, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
                   xticklabels=False, yticklabels=False)
        ax.set_title(title)
        ax.set_xlabel('Distance from TSS')
        
        # Add custom x-axis labels at proper positions
        num_bins = data.shape[1]
        ax.set_xticks([0, num_bins//2, num_bins])
        ax.set_xticklabels(['-2.5kb', 'TSS', '+2.5kb'])
    
    # Find global max for consistent scale
    vmax = max(np.percentile(bm_bivalent, 99), np.percentile(bm_nonbivalent, 99))
    
    # Plot BM signal for both groups with same scale
    create_heatmap(bm_bivalent, ax1, 'Bivalent Genes', vmax=vmax)
    create_heatmap(bm_nonbivalent, ax2, 'Non-bivalent Genes', vmax=vmax)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Set paths
    bg_files = [
        "results/bigwig/BG1_CPM.bw",
        "results/bigwig/BG2_CPM.bw",
        "results/bigwig/BG3_CPM.bw"
    ]
    bm_file = "results/bigwig/BM3_CPM.bw"
    gtf_file = "data/gencode.vM10.annotation.gtf.gz"
    bivalent_genes = "Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv"
    non_bivalent_genes = "Gene_lists/bivalent/expressed_targeted_non_bivalent_NPCs_1000.csv"
    
    # Load gene lists
    print("Loading gene lists...")
    bivalent_list = load_gene_list(bivalent_genes)
    nonbivalent_list = load_gene_list(non_bivalent_genes)
    
    # Get TSS regions
    print("Extracting TSS regions...")
    bivalent_regions = get_tss_regions(gtf_file, bivalent_list)
    nonbivalent_regions = get_tss_regions(gtf_file, nonbivalent_list)
    
    # Process bigWig files
    print("Processing bigWig files for bivalent genes...")
    bg_bivalent, bm_bivalent = process_bigwig_files(bg_files, bm_file, bivalent_regions)
    
    print("Processing bigWig files for non-bivalent genes...")
    bg_nonbivalent, bm_nonbivalent = process_bigwig_files(bg_files, bm_file, nonbivalent_regions)
    
    # Create heatmaps
    print("Generating heatmaps...")
    output_path = "results/heatmaps_comparison_R/bivalent_nonbivalent_heatmaps.pdf"
    plot_heatmaps(bg_bivalent, bm_bivalent, bg_nonbivalent, bm_nonbivalent, output_path)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main() 