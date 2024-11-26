import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

def compare_peaks(peak_count_files, output_file):
    # Read all peak count files
    peak_counts = {}
    first_df = None
    for count_file in peak_count_files:
        sample = os.path.basename(count_file).replace('_peak_counts.txt', '')
        df = pd.read_csv(count_file, sep='\t', 
                        names=['chr', 'start', 'end', 'count'])
        if first_df is None:
            first_df = df
        peak_counts[sample] = df['count']

    # Create a DataFrame with all counts
    counts_df = pd.DataFrame(peak_counts)

    # Calculate mean counts for BG and BM groups
    bg_samples = [col for col in counts_df.columns if col.startswith('BG')]
    bm_samples = [col for col in counts_df.columns if col.startswith('BM')]

    bg_mean = counts_df[bg_samples].mean(axis=1)
    bm_mean = counts_df[bm_samples].mean(axis=1)

    # Perform statistical test (Mann-Whitney U test)
    results = []
    for i in range(len(bg_mean)):
        bg_values = counts_df[bg_samples].iloc[i]
        bm_values = counts_df[bm_samples].iloc[i]
        try:
            statistic, pvalue = stats.mannwhitneyu(bg_values, bm_values, 
                                                 alternative='two-sided')
        except ValueError:  # Handle cases where test cannot be performed
            pvalue = 1.0
        
        # Avoid division by zero in log2 fold change
        if bg_mean[i] == 0:
            bg_mean[i] = 1e-10
        if bm_mean[i] == 0:
            bm_mean[i] = 1e-10
            
        fold_change = np.log2(bm_mean[i] / bg_mean[i])
        
        results.append({
            'chr': first_df['chr'].iloc[i],
            'start': first_df['start'].iloc[i],
            'end': first_df['end'].iloc[i],
            'bg_mean': bg_mean[i],
            'bm_mean': bm_mean[i],
            'log2_fold_change': fold_change,
            'pvalue': pvalue
        })

    # Create results DataFrame and add multiple testing correction
    results_df = pd.DataFrame(results)
    _, results_df['padj'] = multipletests(results_df['pvalue'], 
                                        method='fdr_bh')[:2]

    # Sort by adjusted p-value
    results_df = results_df.sort_values('padj')

    # Save results
    results_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Compare peak sizes between BG and BM samples')
    parser.add_argument('--peak-counts', nargs='+', required=True,
                        help='List of peak count files')
    parser.add_argument('--output', required=True,
                        help='Output comparison file')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    compare_peaks(args.peak_counts, args.output)

if __name__ == '__main__':
    main() 