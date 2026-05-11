#!/usr/bin/env python3

#' Compare SMARCB1 binding profiles between targeted and non-targeted genes using heatmaps
#' 
#' Input files:
#' - results/bigwig/BG1_CPM.bw: Control replicate 1 bigWig file
#' - results/bigwig/BG2_CPM.bw: Control replicate 2 bigWig file 
#' - results/bigwig/BG3_CPM.bw: Control replicate 3 bigWig file
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/targeted/expressed_targeted_targeted_NPCs_1000.csv: List of targeted genes
#' - Gene_lists/targeted/expressed_targeted_non_targeted_NPCs_1000.csv: List of non-targeted genes
#'
#' Output files:
#' - results/metaprofiles_comparison_R/targeted_nontargeted_heatmaps_bg.pdf: Heatmap comparing BG signal at targeted vs non-targeted genes.
#' - results/metaprofiles_comparison_R/targeted_nontargeted_heatmaps_bm.pdf: Heatmap comparing BM signal at targeted vs non-targeted genes.
#' - results/metaprofiles_comparison_R/targeted_nontargeted_metaprofile_bg.pdf: Metaprofile comparing average BG signal at targeted vs non-targeted genes.
#' - results/metaprofiles_comparison_R/targeted_nontargeted_metaprofile_bm.pdf: Metaprofile comparing average BM signal at targeted vs non-targeted genes.
#'
#' The script generates heatmaps and metaprofiles comparing SMARCB1 (BM) and Control (BG)
#' binding profiles around TSS regions between targeted and non-targeted genes.

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use('Agg') # Use the Agg backend for non-interactive plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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

def process_bigwig_files_bg(bg_files: List[str], regions: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Process multiple bigWig files and calculate mean signal."""
    # Process BG files in parallel
    with Pool() as pool:
        bg_matrices = pool.map(partial(extract_signal, regions=regions), bg_files)
    bg_matrix = np.mean(bg_matrices, axis=0)
    
    return bg_matrix # Corrected return statement

def process_bigwig_files_bm(bm_file: str, regions: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Process multiple bigWig files and calculate mean signal."""
    # Process BM file
    bm_matrix = extract_signal(bm_file, regions)
    
    return bm_matrix

# --- Define and Register Custom Colormaps ---
# White to Dark Yellow for BG/GFP
colors_bg = ["#FFFFFF", "#B8860B"] # White to DarkGoldenrod
cmap_name_bg = 'white_to_dark_yellow'
cm_bg = LinearSegmentedColormap.from_list(cmap_name_bg, colors_bg)
matplotlib.colormaps.register(cmap=cm_bg) # Correct way to register colormap
# -----------------------------------------

def plot_heatmaps_bg(bg_targeted: np.ndarray, bg_nontargeted: np.ndarray,
                  output_path: str):
    """Create and save comparative heatmaps for BG (GFP/Control) signal."""
    # Set up the figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 6)) # Adjusted figsize slightly
    fig.suptitle('GFP (Control) Signal at TSS', fontsize=16)
    
    # Choose colormap (custom white to dark yellow)
    cmap = plt.get_cmap('white_to_dark_yellow')
    # Removed cmap.set_under()

    # Combine data to find global scale, excluding NaNs
    # Set fixed vmin and vmax
    vmin, vmax = 0, 4

    # Function to create heatmap (uses shared scale)
    def create_heatmap(data, ax, title, shared_vmin, shared_vmax):
        # Sort data by mean signal intensity for better visualization
        sorted_indices = np.argsort(np.mean(data, axis=1))[::-1]
        sorted_data = data[sorted_indices]
        
        im = ax.imshow(sorted_data, aspect='auto', cmap=cmap, vmin=shared_vmin, vmax=shared_vmax,
                       interpolation='nearest')
        ax.set_title(title)
        ax.set_xlabel('Distance from TSS')
        ax.set_ylabel('Genes')
        
        # Add custom x-axis labels at proper positions
        num_bins = data.shape[1]
        ax.set_xticks([0, num_bins//2, num_bins-1]) # Adjust last tick to be within bounds
        ax.set_xticklabels(['-2.5kb', 'TSS', '+2.5kb'])
        ax.set_yticks([]) # Hide y-axis ticks (gene indices)
        return im

    # Plot heatmaps for both groups using the shared scale
    im1 = create_heatmap(bg_targeted, axes[0], 'Targeted Genes', vmin, vmax)
    im2 = create_heatmap(bg_nontargeted, axes[1], 'Non-targeted Genes', vmin, vmax)
    
    # Adjust layout to make space for shared colorbar
    plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1) 
    
    # Add a shared colorbar in a dedicated axis
    cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7]) # [left, bottom, width, height]
    fig.colorbar(im1, cax=cbar_ax, label='Signal Intensity (CPM)')

    plt.savefig(output_path, dpi=300, bbox_inches='tight') # Save as PDF
    png_output_path = output_path.replace('.pdf', '.png')
    plt.savefig(png_output_path, dpi=300, bbox_inches='tight') # Save as PNG
    plt.close()

def plot_heatmaps_bm(bm_targeted: np.ndarray, bm_nontargeted: np.ndarray,
                  output_path: str):
    """Create and save comparative heatmaps for BM (MeCP2/SMARCB1) signal."""
    # Set up the figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 6)) # Adjusted figsize slightly
    fig.suptitle('MeCP2 (SMARCB1) Signal at TSS', fontsize=16)

    # Choose colormap (white to dark orange)
    cmap = plt.get_cmap('Oranges')
    # Removed cmap.set_under()

    # Combine data to find global scale, excluding NaNs
    # Set fixed vmin and vmax
    vmin, vmax = 0, 8

    # Function to create heatmap (uses shared scale)
    def create_heatmap(data, ax, title, shared_vmin, shared_vmax):
        # Sort data by mean signal intensity for better visualization
        sorted_indices = np.argsort(np.mean(data, axis=1))[::-1]
        sorted_data = data[sorted_indices]
        
        im = ax.imshow(sorted_data, aspect='auto', cmap=cmap, vmin=shared_vmin, vmax=shared_vmax,
                       interpolation='nearest')
        ax.set_title(title)
        ax.set_xlabel('Distance from TSS')
        ax.set_ylabel('Genes')
        
        # Add custom x-axis labels at proper positions
        num_bins = data.shape[1]
        ax.set_xticks([0, num_bins//2, num_bins-1]) # Adjust last tick to be within bounds
        ax.set_xticklabels(['-2.5kb', 'TSS', '+2.5kb'])
        ax.set_yticks([]) # Hide y-axis ticks (gene indices)
        return im

    # Plot heatmaps for both groups using the shared scale
    im1 = create_heatmap(bm_targeted, axes[0], 'Targeted Genes', vmin, vmax)
    im2 = create_heatmap(bm_nontargeted, axes[1], 'Non-targeted Genes', vmin, vmax)

    # Adjust layout to make space for shared colorbar
    plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)

    # Add a shared colorbar in a dedicated axis
    cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7]) # [left, bottom, width, height]
    fig.colorbar(im1, cax=cbar_ax, label='Signal Intensity (CPM)')

    plt.savefig(output_path, dpi=300, bbox_inches='tight') # Save as PDF
    png_output_path = output_path.replace('.pdf', '.png')
    plt.savefig(png_output_path, dpi=300, bbox_inches='tight') # Save as PNG
    plt.close()

def plot_metaprofiles(targeted_mean: np.ndarray, nontargeted_mean: np.ndarray,
                      label1: str, label2: str, title: str, ylabel: str,
                      output_path: str, color1: str = 'darkorange', color2: str = 'grey'):
    """Plot comparative metaprofiles."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    num_bins = len(targeted_mean)
    x_values = np.linspace(-2.5, 2.5, num_bins) # Assuming 5kb window centered at TSS

    ax.plot(x_values, targeted_mean, label=label1, color=color1, linewidth=2)
    ax.plot(x_values, nontargeted_mean, label=label2, color=color2, linewidth=2)
    
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Add vertical line at TSS
    ax.axvline(0, color='black', linestyle=':', linewidth=1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    png_output_path = output_path.replace('.pdf', '.png')
    plt.savefig(png_output_path, dpi=300)
    plt.close()


def main():
    # Set paths
    bg_files = [
        "results/bigwig/BG1_CPM.bw",
        "results/bigwig/BG2_CPM.bw",
        "results/bigwig/BG3_CPM.bw"
    ]
    bm_file = "results/bigwig/BM3_CPM.bw"
    gtf_file = "data/gencode.vM10.basic.annotation.gtf.gz"
    # targeted_genes = "Gene_lists/targeted/expressed_targeted_targeted_NPCs_1000.csv"
    # non_targeted_genes = "Gene_lists/targeted/expressed_targeted_non_targeted_NPCs_1000.csv"
    targeted_genes = "Gene_lists/targets/all_targets_final.csv"
    non_targeted_genes = "Gene_lists/targets/all_no_targets_mm10.csv"
    
    # Load gene lists
    print("Loading gene lists...")
    targeted_list = load_gene_list(targeted_genes)
    nontargeted_list = load_gene_list(non_targeted_genes)
    
    # Get TSS regions
    print("Extracting TSS regions...")
    targeted_regions = get_tss_regions(gtf_file, targeted_list)
    nontargeted_regions = get_tss_regions(gtf_file, nontargeted_list)
    
    # Process bigWig files
    print("Processing bigWig files for targeted genes...")
    bg_targeted = process_bigwig_files_bg(bg_files, targeted_regions)
    bm_targeted = process_bigwig_files_bm(bm_file, targeted_regions) # Calling the function
    
    print("Processing bigWig files for non-targeted genes...")
    bg_nontargeted = process_bigwig_files_bg(bg_files, nontargeted_regions)
    bm_nontargeted = process_bigwig_files_bm(bm_file, nontargeted_regions) # Calling the function
    
    # Create heatmaps
    print("Generating heatmaps...")
    output_path_bg = "results/metaprofiles_comparison_R/targeted_nontargeted_heatmaps_bg.pdf"
    output_path_bm = "results/metaprofiles_comparison_R/targeted_nontargeted_heatmaps_bm.pdf"
    plot_heatmaps_bg(bg_targeted, bg_nontargeted, output_path_bg)
    plot_heatmaps_bm(bm_targeted, bm_nontargeted, output_path_bm)

    # Calculate mean signals for metaprofiles
    print("Calculating mean signals for metaprofiles...")
    bg_targeted_mean = np.nanmean(bg_targeted, axis=0)
    bg_nontargeted_mean = np.nanmean(bg_nontargeted, axis=0)
    bm_targeted_mean = np.nanmean(bm_targeted, axis=0)
    bm_nontargeted_mean = np.nanmean(bm_nontargeted, axis=0)

    # Create metaprofiles
    print("Generating metaprofiles...")
    output_path_meta_bg = "results/metaprofiles_comparison_R/targeted_nontargeted_metaprofile_bg.pdf"
    output_path_meta_bm = "results/metaprofiles_comparison_R/targeted_nontargeted_metaprofile_bm.pdf"

    plot_metaprofiles(bg_targeted_mean, bg_nontargeted_mean,
                      'Targeted', 'Non-targeted',
                      'Average GFP (Control) Signal at TSS', 'Mean Signal (CPM)',
                      output_path_meta_bg, color1='yellow', color2='grey') # Changed color1 to yellow

    plot_metaprofiles(bm_targeted_mean, bm_nontargeted_mean,
                      'Targeted', 'Non-targeted',
                      'Average MeCP2 (SMARCB1) Signal at TSS', 'Mean Signal (CPM)',
                      output_path_meta_bm, color1='orange', color2='grey') # Changed color1 to orange

    print("Analysis completed successfully!")

if __name__ == "__main__":
    main()
