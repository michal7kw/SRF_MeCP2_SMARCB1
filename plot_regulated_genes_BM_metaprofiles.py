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
#' - results/metaprofiles_comparison_R/regulated_genes_BM_profile.pdf: 
#'   Plot showing SMARCB1 binding profiles for different gene categories

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

def calculate_profile_stats(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate mean and standard error for the profile."""
    mean_profile = np.mean(matrix, axis=0)
    se_profile = np.std(matrix, axis=0) / np.sqrt(matrix.shape[0])
    return mean_profile, se_profile

def plot_profiles(profiles: Dict, output_path: str):
    """Create and save profile plot."""
    fig, ax = plt.subplots(figsize=(10, 6))
    x_pos = np.linspace(-2500, 2500, profiles['up'][0].shape[0])
    
    # Colors for different categories
    colors = {
        'up': '#ff7f0e',
        'down': '#2ca02c',
        'not': '#1f77b4'
    }
    
    # Plot signal profiles
    for category in ['up', 'down', 'not']:
        mean_profile, se_profile = profiles[category]
        
        # Plot mean profile
        ax.plot(x_pos, mean_profile, color=colors[category],
                label=f'{category.capitalize()} regulated')
        
        # Plot standard error bands
        ax.fill_between(x_pos, mean_profile - se_profile, mean_profile + se_profile,
                       color=colors[category], alpha=0.2)
    
    # Customize plot
    ax.set_title('SMARCB1 binding around TSS')
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average RPKM')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Set paths
    bm_file = "results/bigwig/BM3_CPM.bw"
    gtf_file = "data/gencode.vM10.annotation.gtf.gz"
    
    # Gene lists
    gene_files = {
        'up': "Gene_lists/targets/all_targets_final_up_regulated.csv",
        'down': "Gene_lists/targets/all_targets_final_down_regulated.csv",
        'not': "Gene_lists/targets/all_targets_final_not_regulated.csv"
    }
    
    # Process each category
    profiles = {}
    
    for category, file_path in gene_files.items():
        print(f"Processing {category} regulated genes...")
        
        # Load genes and get regions
        genes = load_gene_list(file_path)
        regions = get_tss_regions(gtf_file, genes)
        
        # Process SMARCB1 signal
        bm_matrix = extract_signal(bm_file, regions)
        profiles[category] = calculate_profile_stats(bm_matrix)
        
        # Print summary statistics
        print(f"\nSummary Statistics for {category} regulated genes:")
        print(f"Number of genes: {len(genes)}")
        print(f"Mean SMARCB1 signal: {np.mean(bm_matrix):.3f}\n")
    
    # Create plot
    print("Generating plot...")
    output_path = "results/metaprofiles_comparison_R/regulated_genes_BM_profile.pdf"
    plot_profiles(profiles, output_path)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main() 