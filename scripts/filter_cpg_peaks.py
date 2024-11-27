import argparse
import subprocess
import os
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def filter_cpg_peaks(peaks_file, cpg_islands_file, output_file):
    """Filter peaks to only those overlapping with CpG islands."""
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create temporary directory
    temp_dir = os.path.join(os.path.dirname(output_file), 'tmp')
    os.makedirs(temp_dir, exist_ok=True)
    
    # Temporary files
    temp_intersect = os.path.join(temp_dir, "cpg_intersect.tmp")
    
    try:
        # Find peaks that overlap with CpG islands
        logger.info("Finding peaks overlapping with CpG islands...")
        intersect_cmd = f"bedtools intersect -a {peaks_file} -b {cpg_islands_file} -wa -u > {temp_intersect}"
        subprocess.run(intersect_cmd, shell=True, check=True)
        
        # Read the filtered peaks
        df = pd.read_csv(temp_intersect, sep='\t', header=None)
        
        # Save filtered peaks
        df.to_csv(output_file, sep='\t', header=False, index=False)
        
        # Log statistics
        total_peaks = sum(1 for _ in open(peaks_file))
        cpg_peaks = len(df)
        logger.info(f"Total peaks: {total_peaks}")
        logger.info(f"Peaks overlapping CpG islands: {cpg_peaks} ({cpg_peaks/total_peaks*100:.1f}%)")
        
    finally:
        # Cleanup
        if os.path.exists(temp_intersect):
            os.remove(temp_intersect)
        try:
            os.rmdir(temp_dir)
        except:
            pass

def main():
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