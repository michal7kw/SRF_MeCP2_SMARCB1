import argparse
import pandas as pd
import numpy as np
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def validate_input_data(count_files):
    """Validate input count files before processing."""
    for file in count_files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"Count file not found: {file}")
        
        df = pd.read_csv(file, sep='\t', names=['chr', 'start', 'end', 'count'])
        logger.info(f"File: {file}")
        logger.info(f"Shape: {df.shape}")
        logger.info(f"Missing values: {df.isnull().sum().sum()}")
        
        if df.empty or df['count'].isnull().all():
            raise ValueError(f"No valid count data in {file}")

def compare_peaks(peak_count_files, output_file, threads=1):
    """Compare normalized peak counts between conditions."""
    
    # Read and process peak counts
    peak_counts = {}
    first_df = None
    
    for count_file in peak_count_files:
        sample = os.path.basename(count_file).replace('_promoter_counts.txt', '')
        logger.info(f"Reading {count_file}...")
        
        df = pd.read_csv(count_file, sep='\t')
        if first_df is None:
            first_df = df[['chr', 'start', 'end']]
        peak_counts[sample] = df['count']  # Using normalized counts
    
    # Create counts DataFrame
    counts_df = pd.DataFrame(peak_counts)
    
    # Calculate means for BG and BM
    bg_mean = counts_df[[col for col in counts_df.columns if col.startswith('BG')]].mean(axis=1)
    bm_mean = counts_df[[col for col in counts_df.columns if col.startswith('BM')]].mean(axis=1)
    
    # Add pseudocount to avoid division by zero
    pseudocount = 1.0
    
    # Calculate log2 fold changes
    results_df = pd.DataFrame({
        'chr': first_df['chr'],
        'start': first_df['start'],
        'end': first_df['end'],
        'bg_mean': bg_mean,
        'bm_mean': bm_mean,
        'log2_fold_change': np.log2((bm_mean + pseudocount) / (bg_mean + pseudocount))
    })
    
    # Sort by absolute fold change
    results_df['abs_fc'] = abs(results_df['log2_fold_change'])
    results_df = results_df.sort_values('abs_fc', ascending=False)
    results_df = results_df.drop('abs_fc', axis=1)
    
    # Save results
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Log summary
    logger.info("Summary of comparisons:")
    logger.info(f"Mean BG signal: {bg_mean.mean():.2f}")
    logger.info(f"Mean BM signal: {bm_mean.mean():.2f}")
    logger.info(f"Peaks up in BM (log2FC > 1): {(results_df['log2_fold_change'] > 1).sum()}")
    logger.info(f"Peaks down in BM (log2FC < -1): {(results_df['log2_fold_change'] < -1).sum()}")

def main():
    parser = argparse.ArgumentParser(description='Compare peak sizes between BG and BM samples')
    parser.add_argument('--peak-counts', nargs='+', required=True,
                        help='List of peak count files')
    parser.add_argument('--output', required=True,
                        help='Output comparison file')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use')
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Run analysis
        compare_peaks(args.peak_counts, args.output, threads=args.threads)
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == '__main__':
    main() 