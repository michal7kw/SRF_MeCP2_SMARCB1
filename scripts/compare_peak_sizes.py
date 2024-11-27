import argparse
import pandas as pd
import numpy as np
import os
import logging
import matplotlib.pyplot as plt

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

def compare_peaks(peak_count_files, output_file, sample_name, threads=1):
    """Compare normalized peak counts between conditions."""
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read and process peak counts
    peak_counts = {}
    first_df = None
    
    for count_file in peak_count_files:
        sample = os.path.basename(count_file).replace('_promoter_counts.txt', '')
        logger.info(f"Reading {count_file}...")
        
        try:
            df = pd.read_csv(count_file, sep='\t')
            if first_df is None:
                first_df = df[['chr', 'start', 'end']]
            peak_counts[sample] = df['count']
        except Exception as e:
            logger.error(f"Error reading {count_file}: {str(e)}")
            raise
    
    # Create counts DataFrame
    counts_df = pd.DataFrame(peak_counts)
    
    # Calculate means for BG and BM
    bg_mean = counts_df[[col for col in counts_df.columns if col.startswith('BG')]].mean(axis=1)
    bm_mean = counts_df[[col for col in counts_df.columns if col.startswith('BM')]].mean(axis=1)
    
    # Calculate mean intensity for filtering
    mean_intensity = (bg_mean + bm_mean) / 2
    
    # Filter out low count regions
    min_count = 5  # Minimum count threshold
    mask = (bg_mean >= min_count) | (bm_mean >= min_count)
    
    # Calculate log2 fold changes with filtering and improved normalization
    results_df = pd.DataFrame({
        'chr': first_df['chr'],
        'start': first_df['start'],
        'end': first_df['end'],
        'bg_mean': bg_mean,
        'bm_mean': bm_mean,
        'mean_intensity': mean_intensity,
        'passed_filter': mask
    })
    
    # Calculate log2 fold changes only for peaks that pass filter
    pseudocount = 1.0
    results_df['log2_fold_change'] = np.where(
        results_df['passed_filter'],
        np.log2((results_df['bm_mean'] + pseudocount) / (results_df['bg_mean'] + pseudocount)),
        np.nan
    )
    
    # Add quality metrics
    results_df['coefficient_of_variation'] = np.where(
        results_df['mean_intensity'] > 0,
        np.sqrt(np.var([bg_mean, bm_mean], axis=0)) / results_df['mean_intensity'],
        np.nan
    )
    
    # Sort by absolute fold change for filtered peaks
    results_df['abs_fc'] = abs(results_df['log2_fold_change'])
    results_df = results_df.sort_values('abs_fc', ascending=False)
    results_df = results_df.drop(['abs_fc', 'passed_filter'], axis=1)
    
    # Save results
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Log summary statistics
    total_peaks = len(results_df)
    filtered_peaks = results_df['log2_fold_change'].notna().sum()
    
    logger.info("Summary of comparisons:")
    logger.info(f"Total peaks analyzed: {total_peaks}")
    logger.info(f"Peaks passing filters: {filtered_peaks} ({filtered_peaks/total_peaks*100:.1f}%)")
    logger.info(f"Mean BG signal: {bg_mean.mean():.2f}")
    logger.info(f"Mean BM signal: {bm_mean.mean():.2f}")
    logger.info(f"Peaks up in BM (log2FC > 1): {(results_df['log2_fold_change'] > 1).sum()}")
    logger.info(f"Peaks down in BM (log2FC < -1): {(results_df['log2_fold_change'] < -1).sum()}")
    
    # Generate QC plot
    plt.figure(figsize=(10, 5))
    plt.hist(results_df['coefficient_of_variation'].dropna(), bins=50)
    plt.xlabel('Coefficient of Variation')
    plt.ylabel('Count')
    plt.title('Distribution of Peak Variability')
    plt.savefig(f"{output_file}_qc.png")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare peak sizes between BG and BM samples')
    parser.add_argument('--peak-counts', nargs='+', required=True,
                        help='List of peak count files')
    parser.add_argument('--output', required=True,
                        help='Output comparison file')
    parser.add_argument('--sample-name', required=True,
                        help='Sample name')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use')
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Run analysis
        compare_peaks(args.peak_counts, args.output, args.sample_name, threads=args.threads)
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == '__main__':
    main() 