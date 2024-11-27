import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os

def reanalyze_peaks(bg_counts, bm_counts, output_dir):
    """Perform a more thorough peak analysis."""
    
    # Read count files
    bg_df = pd.read_csv(bg_counts, sep='\t')
    bm_df = pd.read_csv(bm_counts, sep='\t')
    
    # Merge data
    merged_df = pd.merge(
        bg_df, bm_df,
        on=['chr', 'start', 'end'],
        suffixes=('_bg', '_bm')
    )
    
    # Calculate statistics
    merged_df['log2fc'] = np.log2((merged_df['count_bm'] + 1) / 
                                 (merged_df['count_bg'] + 1))
    
    # Perform statistical test
    pvalues = []
    for _, row in merged_df.iterrows():
        try:
            _, pval = stats.mannwhitneyu(
                [row['count_bg']], [row['count_bm']],
                alternative='two-sided'
            )
            pvalues.append(pval)
        except:
            pvalues.append(np.nan)
    
    merged_df['pvalue'] = pvalues
    merged_df['padj'] = multipletests(
        merged_df['pvalue'].fillna(1), 
        method='fdr_bh'
    )[1]
    
    # Create MA plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=merged_df,
        x=np.log2(merged_df['count_bg'] + merged_df['count_bm'])/2,
        y='log2fc',
        alpha=0.5
    )
    plt.title('MA Plot')
    plt.xlabel('Average log2 count')
    plt.ylabel('log2 fold change')
    plt.savefig(os.path.join(output_dir, 'ma_plot.png'))
    
    # Save results
    merged_df.to_csv(
        os.path.join(output_dir, 'reanalyzed_peaks.txt'),
        sep='\t',
        index=False
    )
    
    # Generate summary
    significant_peaks = merged_df[merged_df['padj'] < 0.05]
    summary = {
        'Total peaks': len(merged_df),
        'Significant peaks': len(significant_peaks),
        'Upregulated': len(significant_peaks[significant_peaks['log2fc'] > 0]),
        'Downregulated': len(significant_peaks[significant_peaks['log2fc'] < 0])
    }
    
    with open(os.path.join(output_dir, 'analysis_summary.txt'), 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    
    return merged_df

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--bg-counts', required=True, help='Background peak counts')
    parser.add_argument('--bm-counts', required=True, help='Treatment peak counts')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    reanalyze_peaks(args.bg_counts, args.bm_counts, args.output_dir)

if __name__ == '__main__':
    main() 