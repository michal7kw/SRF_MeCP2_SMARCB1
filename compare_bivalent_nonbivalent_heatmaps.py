#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import gtfparse
import warnings
warnings.filterwarnings('ignore')

def load_gene_list(file_path):
    """Load gene list from CSV file."""
    return pd.read_csv(file_path, header=None)[0].tolist()

def parse_gtf(gtf_file):
    """Parse GTF file and extract gene information."""
    df = gtfparse.read_gtf(gtf_file)
    gene_df = df[df['feature'] == 'gene'].copy()
    return gene_df

def get_tss_regions(genes, gene_df, upstream=2500, downstream=2500):
    """Get TSS regions for given genes."""
    gene_info = gene_df[gene_df['gene_name'].isin(genes)].copy()
    
    tss_regions = []
    for _, row in gene_info.iterrows():
        if row['strand'] == '+':
            tss_start = row['start'] - upstream
            tss_end = row['start'] + downstream
        else:
            tss_start = row['end'] - downstream
            tss_end = row['end'] + upstream
            
        tss_regions.append({
            'chr': row['seqname'],
            'start': tss_start,
            'end': tss_end,
            'gene_name': row['gene_name'],
            'strand': row['strand']
        })
    
    return pd.DataFrame(tss_regions)

def calculate_signal_matrix(bw_file, regions, bins=100):
    """Calculate signal matrix for given regions using bigWig file."""
    bw = pyBigWig.open(bw_file)
    matrix = np.zeros((len(regions), bins))
    
    for i, row in regions.iterrows():
        try:
            values = bw.stats(row['chr'], 
                            int(row['start']), 
                            int(row['end']), 
                            nBins=bins,
                            type="mean")
            values = [0 if v is None else v for v in values]
            matrix[i] = values
        except:
            continue
            
    bw.close()
    return matrix

def define_enrichment_methods(total_genome_peaks):
    """Define different methods to calculate enrichment between samples."""
    methods = {
        'signal_ratio': lambda bg, bm: (
            np.mean(bm) / np.mean(bg)
            if np.mean(bg) > 0
            else float('inf') if np.mean(bm) > 0
            else 0.0
        ),
        
        'peak_count': lambda bg, bm: (
            np.sum(bm > 0) / max(np.sum(bg > 0), 1)
        ),
        
        'combined_score': lambda bg, bm: (
            (np.mean(bm) * np.sum(bm > 0)) / 
            max((np.mean(bg) * np.sum(bg > 0)), 1)
        ),
        
        'statistical': lambda bg, bm: (
            stats.fisher_exact([
                [np.sum(bm > 0), np.sum(bg > 0)],
                [total_genome_peaks - np.sum(bm > 0), 
                 total_genome_peaks - np.sum(bg > 0)]
            ])[1]
        ),
        
        'coverage_score': lambda bg, bm: (
            (np.sum(bm) / len(bm)) / max((np.sum(bg) / len(bg)), 1)
        )
    }
    return methods

def calculate_enrichment(bg_matrix, bm_matrix, method_func):
    """Calculate enrichment scores for each region."""
    n_regions = bg_matrix.shape[0]
    enrichment_scores = np.zeros(n_regions)
    
    for i in range(n_regions):
        enrichment_scores[i] = method_func(bg_matrix[i], bm_matrix[i])
    
    return enrichment_scores

def plot_heatmap(enrichment_scores, gene_list_name, method_name, output_dir):
    """Create and save heatmap for enrichment scores."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(enrichment_scores.reshape(-1, 1),
                cmap='YlOrRd',
                xticklabels=False,
                yticklabels=False)
    plt.title(f'{method_name} Enrichment\n{gene_list_name}')
    plt.xlabel('Position relative to TSS')
    plt.ylabel('Genes')
    
    output_path = Path(output_dir) / f'{gene_list_name}_{method_name}_heatmap.pdf'
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def main():
    # Set paths
    bg_files = [
        "results/bigwig/BG1_RPKM.bw",
        "results/bigwig/BG2_RPKM.bw",
        "results/bigwig/BG3_RPKM.bw"
    ]
    bm_file = "results/bigwig/BM3_RPKM.bw"
    gtf_file = "data/gencode.vM10.annotation.gtf.gz"
    
    # Gene lists
    gene_lists = {
        'bivalent': "Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv",
        'non_bivalent': "Gene_lists/bivalent/expressed_targeted_non_bivalent_NPCs_1000.csv"
    }
    
    # Create output directory
    output_dir = Path("results/heatmaps_comparison")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse GTF file
    print("Parsing GTF file...")
    gene_df = parse_gtf(gtf_file)
    
    # Process each gene list
    for list_name, gene_list_file in gene_lists.items():
        print(f"\nProcessing {list_name} genes...")
        genes = load_gene_list(gene_list_file)
        regions = get_tss_regions(genes, gene_df)
        
        # Calculate signal matrices
        print("Calculating BG signal matrix...")
        bg_matrices = [calculate_signal_matrix(f, regions) for f in bg_files]
        bg_matrix = np.mean(bg_matrices, axis=0)
        
        print("Calculating BM signal matrix...")
        bm_matrix = calculate_signal_matrix(bm_file, regions)
        
        # Calculate enrichment using different methods
        methods = define_enrichment_methods(len(gene_df))
        
        for method_name, method_func in methods.items():
            print(f"Calculating {method_name} enrichment...")
            enrichment_scores = calculate_enrichment(bg_matrix, bm_matrix, method_func)
            
            # Create heatmap
            print(f"Creating heatmap for {method_name}...")
            plot_heatmap(enrichment_scores, list_name, method_name, output_dir)

if __name__ == "__main__":
    main() 