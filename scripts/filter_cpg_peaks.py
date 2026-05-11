"""
This script filters genomic peaks to identify those overlapping with CpG islands.

Key features:
- Takes BED format files for peaks and CpG islands as input
- Uses bedtools intersect to find overlapping regions
- Outputs filtered peaks that overlap CpG islands
- Provides statistics on filtering results
- Handles temporary files and cleanup
- Includes detailed logging

Input:
- Peak regions in BED format
- CpG islands in BED format 
- Output file path for filtered peaks

Output:
- BED file containing only peaks that overlap CpG islands
- Logging information with filtering statistics
"""

import argparse
import subprocess
import os
import logging
import pandas as pd

# Configure logging to show informational messages
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def filter_cpg_peaks(peaks_file, cpg_islands_file, output_file):
    """
    Filter peaks to only those overlapping with CpG islands.
    
    Args:
        peaks_file (str): Path to input peaks BED file
        cpg_islands_file (str): Path to CpG islands BED file
        output_file (str): Path to output filtered peaks file
        
    The function:
    1. Creates temporary working directory
    2. Uses bedtools to find peaks overlapping CpG islands
    3. Saves filtered peaks to output file
    4. Calculates and logs filtering statistics
    5. Cleans up temporary files
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create temporary directory for intermediate files
    temp_dir = os.path.join(os.path.dirname(output_file), 'tmp')
    os.makedirs(temp_dir, exist_ok=True)
    
    # Set path for temporary intersection results
    temp_intersect = os.path.join(temp_dir, "cpg_intersect.tmp")
    
    try:
        # Use bedtools intersect to find peaks overlapping CpG islands
        # -wa outputs the original peak entry
        # -u outputs each peak only once even if multiple overlaps
        logger.info("Finding peaks overlapping with CpG islands...")
        intersect_cmd = f"bedtools intersect -a {peaks_file} -b {cpg_islands_file} -wa -u > {temp_intersect}"
        subprocess.run(intersect_cmd, shell=True, check=True)
        
        # Read the filtered peaks into pandas DataFrame
        df = pd.read_csv(temp_intersect, sep='\t', header=None)
        
        # Write filtered peaks to output file
        df.to_csv(output_file, sep='\t', header=False, index=False)
        
        # Calculate and log filtering statistics
        total_peaks = sum(1 for _ in open(peaks_file))
        cpg_peaks = len(df)
        logger.info(f"Total peaks: {total_peaks}")
        logger.info(f"Peaks overlapping CpG islands: {cpg_peaks} ({cpg_peaks/total_peaks*100:.1f}%)")
        
    finally:
        # Clean up temporary files and directory
        if os.path.exists(temp_intersect):
            os.remove(temp_intersect)
        try:
            os.rmdir(temp_dir)
        except:
            pass

def main():
    """Parse command line arguments and run the filtering pipeline."""
    parser = argparse.ArgumentParser(description='Filter peaks for CpG island overlap')
    parser.add_argument('--peaks', required=True,
                        help='Input peaks file (BED format)')
    parser.add_argument('--cpg-islands', required=True,
                        help='CpG islands file (BED format)')
    parser.add_argument('--output', required=True,
                        help='Output filtered peaks file')
    
    args = parser.parse_args()
    
    try:
        filter_cpg_peaks(args.peaks, args.cpg_islands, args.output)
    except Exception as e:
        logger.error(f"Error filtering peaks: {str(e)}")
        raise

if __name__ == '__main__':
    main()

# Summary:
# This script filters genomic peak regions to identify those that overlap with CpG islands.
# It takes two BED format files as input - one containing peak regions and another with CpG islands.
# Using bedtools intersect, it identifies peaks that overlap with CpG islands and outputs them
# to a new BED file. The script provides detailed logging of the filtering process and statistics
# on how many peaks overlap with CpG islands. It handles temporary files and includes error handling
# for robustness. This filtering is useful for focusing analysis on peaks that occur in CpG-rich
# regions, which are often associated with gene regulatory elements.