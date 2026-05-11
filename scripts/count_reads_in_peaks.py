"""
This script counts and normalizes reads in genomic peak regions from BAM files.

Key features:
- Processes BAM files to count reads overlapping peak regions
- Normalizes counts to reads per million (RPM) 
- Performs quality control checks on read counts
- Handles temporary files and cleanup
- Provides detailed logging and error handling

Input:
- Peak regions in BED format
- Aligned reads in BAM format
- Output file path for normalized counts
- Sample name for identification
- Optional number of threads

Output:
- Tab-separated file with normalized read counts per peak
- Logging information and QC metrics
"""

import argparse
import subprocess
import os
import sys
import pandas as pd
import logging

# Set up logging configuration
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def count_reads(peaks_file, bam_file, output_file, sample_name, threads=1):
    """
    Count and normalize reads in peak regions from a BAM file.
    
    Args:
        peaks_file (str): Path to BED file containing peak regions
        bam_file (str): Path to BAM file containing aligned reads
        output_file (str): Path to output file for normalized counts
        sample_name (str): Name identifier for the sample
        threads (int): Number of threads to use (default: 1)
        
    The function:
    1. Creates temporary working directories
    2. Calculates total mapped reads for normalization
    3. Sorts peaks using genome coordinates
    4. Counts reads in peak regions
    5. Normalizes counts to reads per million
    6. Performs QC checks
    7. Cleans up temporary files
    """
    
    # Create output directory and temporary directory
    output_dir = os.path.dirname(output_file)
    temp_dir = os.path.join(output_dir, 'tmp')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Set up paths for temporary files
    temp_prefix = os.path.join(temp_dir, os.path.basename(output_file))
    genome_file = f"{temp_prefix}.genome"
    sorted_peaks = f"{temp_prefix}.sorted.tmp"
    counts_tmp = f"{temp_prefix}.counts.tmp"
    
    # Calculate total mapped reads for RPM normalization
    logger.info("Calculating total mapped reads...")
    total_reads_cmd = f"samtools view -c -F 4 {bam_file}"  # -F 4 excludes unmapped reads
    total_reads = int(subprocess.check_output(total_reads_cmd, shell=True))
    logger.info(f"Total mapped reads: {total_reads}")
    
    # Extract chromosome sizes from BAM header for proper sorting
    logger.info("Creating genome file from BAM header...")
    genome_cmd = f"samtools view -H {bam_file} | grep '^@SQ' | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > {genome_file}"
    subprocess.run(genome_cmd, shell=True, check=True)
    
    # Sort peak regions by genomic coordinates
    logger.info("Sorting BED file...")
    sort_cmd = f"bedtools sort -g {genome_file} -i {peaks_file} > {sorted_peaks}"
    subprocess.run(sort_cmd, shell=True, check=True)
    
    try:
        # Use bedtools coverage to count reads overlapping peaks
        logger.info("Counting reads in peaks...")
        cmd = (f"bedtools coverage -a {sorted_peaks} -b {bam_file} "
               f"-sorted -g {genome_file} "
               f"-counts > {counts_tmp}")
        
        subprocess.run(cmd, shell=True, check=True)
        
        # Load count data and add column names
        df = pd.read_csv(counts_tmp, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'gene', 'raw_count'])
        
        # Convert raw counts to reads per million
        df['count'] = df['raw_count'] * 1e6 / total_reads
        
        # Write normalized counts to output file
        df.to_csv(output_file, sep='\t', index=False)
        
        # Log summary statistics
        logger.info(f"Normalized counts saved to {output_file}")
        logger.info(f"Mean raw count: {df['raw_count'].mean():.2f}")
        logger.info(f"Mean normalized count: {df['count'].mean():.2f}")
        
        # Perform quality control checks
        zero_peaks = (df['raw_count'] == 0).sum()
        if zero_peaks > len(df) * 0.5:
            logger.warning(f"More than 50% of peaks have zero reads: {zero_peaks}/{len(df)}")
        
        if total_reads < 1000000:
            logger.warning(f"Low number of mapped reads: {total_reads}")
            
    finally:
        # Clean up all temporary files
        for tmp_file in [counts_tmp, genome_file, sorted_peaks]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        # Try to remove temp directory if empty
        try:
            os.rmdir(temp_dir)
        except:
            pass

def main():
    """Parse command line arguments and execute read counting."""
    parser = argparse.ArgumentParser(description='Count reads in peaks')
    parser.add_argument('--peaks', required=True,
                        help='Peaks bed file')
    parser.add_argument('--bam', required=True,
                        help='BAM file')
    parser.add_argument('--output', required=True,
                        help='Output counts file')
    parser.add_argument('--sample-name', required=True,
                        help='Sample name')
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads to use')
    
    args = parser.parse_args()
    
    try:
        count_reads(args.peaks, args.bam, args.output, args.sample_name, threads=args.threads)
    except Exception as e:
        logger.error(f"Error processing {args.bam}: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()

"""
Summary:
This script is designed for ChIP-seq data analysis, specifically counting reads in peak regions.
It takes aligned sequencing reads (BAM) and peak regions (BED) as input, then:

1. Processes the BAM file to count total mapped reads
2. Ensures proper sorting of peak regions using chromosome coordinates
3. Uses bedtools to count reads overlapping each peak region
4. Normalizes raw counts to reads per million (RPM) to account for sequencing depth
5. Performs quality control checks for:
   - Peaks with zero reads
   - Low total read counts
6. Outputs a tab-separated file containing:
   - Peak coordinates (chr, start, end)
   - Gene names (if provided)
   - Raw read counts
   - Normalized RPM values
7. Provides detailed logging throughout the process
8. Handles temporary files and cleanup automatically

The script is commonly used in ChIP-seq workflows to quantify protein binding
or histone modification levels at specific genomic regions.
"""