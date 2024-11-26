import pandas as pd
import argparse
import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def prepare_significant_peaks(comparison_file, output_file, fc_threshold=1):
    """
    Extract peaks with significant fold changes.
    fc_threshold: log2 fold change threshold (default=1, equivalent to FC=2)
    """
    # Read comparison file
    df = pd.read_csv(comparison_file, sep='\t')
    logger.info(f"Loaded {len(df)} total peaks")
    
    # Filter significant peaks
    significant = df[abs(df['log2_fold_change']) > fc_threshold].copy()
    
    # Sort by absolute fold change
    significant['abs_fc'] = abs(significant['log2_fold_change'])
    significant = significant.sort_values('abs_fc', ascending=False)
    significant = significant.drop('abs_fc', axis=1)
    
    # Save filtered peaks
    significant.to_csv(output_file, sep='\t', index=False)
    
    # Log summary
    up_regulated = (significant['log2_fold_change'] > fc_threshold).sum()
    down_regulated = (significant['log2_fold_change'] < -fc_threshold).sum()
    
    logger.info(f"Found {len(significant)} significant peaks:")
    logger.info(f"  Up-regulated in BM: {up_regulated}")
    logger.info(f"  Down-regulated in BM: {down_regulated}")
    
    return significant

def main():
    parser = argparse.ArgumentParser(description='Prepare significant peaks for annotation')
    parser.add_argument('--input', required=True,
                        help='Peak comparison file')
    parser.add_argument('--output', required=True,
                        help='Output significant peaks file')
    parser.add_argument('--fc-threshold', type=float, default=1.0,
                        help='Log2 fold change threshold (default=1, equivalent to FC=2)')
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        prepare_significant_peaks(args.input, args.output, args.fc_threshold)
    except Exception as e:
        logger.error(f"Failed to prepare significant peaks: {str(e)}")
        raise

if __name__ == '__main__':
    main() 