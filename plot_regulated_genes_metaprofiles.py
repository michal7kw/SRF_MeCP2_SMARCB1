#!/usr/bin/env python3

#' Compare SMARCB1 binding profiles between up-regulated, down-regulated, and non-regulated genes
#' 
#' Input files:
#' - results/bigwig/BG1_CPM.bw: Control replicate 1 bigWig file
#' - results/bigwig/BG2_CPM.bw: Control replicate 2 bigWig file 
#' - results/bigwig/BG3_CPM.bw: Control replicate 3 bigWig file
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/targets/all_targets_final_up_regulated.csv: Up-regulated genes
#' - Gene_lists/targets/all_targets_final_down_regulated.csv: Down-regulated genes
#' - Gene_lists/targets/all_targets_final_not_regulated.csv: Non-regulated genes
#'
#' Output files:
#' - results/metaprofiles_comparison_R/combined_regulated_genes_profile.pdf: 
#'   Combined plot showing binding profiles and fold changes

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

def process_bigwig_files(files: List[str], regions: pd.DataFrame) -> np.ndarray:
    """Process multiple bigWig files and calculate mean signal."""
    matrices = []
    for file in files:
        matrix = extract_signal(file, regions)
        matrices.append(matrix)
    return np.mean(matrices, axis=0)

def calculate_profile_stats(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate mean and standard error for the profile."""
    mean_profile = np.mean(matrix, axis=0)
    se_profile = np.std(matrix, axis=0) / np.sqrt(matrix.shape[0])
    return mean_profile, se_profile

def plot_profiles(bg_profiles: Dict, bm_profiles: Dict, output_path: str):
    """Create and save combined profile plots."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), height_ratios=[2, 1])
    x_pos = np.linspace(-2500, 2500, bg_profiles['up'][0].shape[0])
    
    # Colors for different categories
    colors = {
        'up': '#ff7f0e',
        'down': '#2ca02c',
        'not': '#1f77b4'
    }
    
    # Plot signal profiles
    for category in ['up', 'down', 'not']:
        bg_mean, bg_se = bg_profiles[category]
        bm_mean, bm_se = bm_profiles[category]
        
        # Plot BG
        ax1.plot(x_pos, bg_mean, color=colors[category], linestyle='--', 
                label=f'{category.capitalize()} regulated - BG')
        ax1.fill_between(x_pos, bg_mean - bg_se, bg_mean + bg_se, 
                        color=colors[category], alpha=0.1)
        
        # Plot BM
        ax1.plot(x_pos, bm_mean, color=colors[category], linestyle='-',
                label=f'{category.capitalize()} regulated - BM')
        ax1.fill_between(x_pos, bm_mean - bm_se, bm_mean + bm_se,
                        color=colors[category], alpha=0.2)
        
        # Calculate and plot fold change
        fold_change = np.log2(bm_mean / bg_mean)
        ax2.plot(x_pos, fold_change, color=colors[category],
                label=f'{category.capitalize()} regulated')
    
    # Customize plots
    ax1.set_title('SMARCB1 binding around TSS')
    ax1.set_xlabel('Distance from TSS (bp)')
    ax1.set_ylabel('Average RPKM')
    ax1.legend()
    
    ax2.set_title('Log2 Fold Change (BM/BG)')
    ax2.set_xlabel('Distance from TSS (bp)')
    ax2.set_ylabel('Log2 Fold Change')
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax2.legend()
    
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
    
    # Gene lists
    gene_files = {
        'up': "Gene_lists/targets/all_targets_final_up_regulated.csv",
        'down': "Gene_lists/targets/all_targets_final_down_regulated.csv",
        'not': "Gene_lists/targets/all_targets_final_not_regulated.csv"
    }
    
    # Process each category
    bg_profiles = {}
    bm_profiles = {}
    
    for category, file_path in gene_files.items():
        print(f"Processing {category} regulated genes...")
        
        # Load genes and get regions
        genes = load_gene_list(file_path)
        regions = get_tss_regions(gtf_file, genes)
        
        # Process background signal (average of replicates)
        bg_matrix = process_bigwig_files(bg_files, regions)
        bg_profiles[category] = calculate_profile_stats(bg_matrix)
        
        # Process SMARCB1 signal
        bm_matrix = extract_signal(bm_file, regions)
        bm_profiles[category] = calculate_profile_stats(bm_matrix)
        
        # Print summary statistics
        print(f"\nSummary Statistics for {category} regulated genes:")
        print(f"Number of genes: {len(genes)}")
        print(f"Mean BG signal: {np.mean(bg_matrix):.3f}")
        print(f"Mean BM signal: {np.mean(bm_matrix):.3f}")
        print(f"Mean log2 fold change: {np.mean(np.log2(np.mean(bm_matrix)/np.mean(bg_matrix))):.3f}\n")
    
    # Create plots
    print("Generating plots...")
    output_path = "results/metaprofiles_comparison_R/combined_regulated_genes_profile.pdf"
    plot_profiles(bg_profiles, bm_profiles, output_path)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main() 