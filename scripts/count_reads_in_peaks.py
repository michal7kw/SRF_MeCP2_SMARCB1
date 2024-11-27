import argparse
import subprocess
import os
import sys
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def count_reads(peaks_file, bam_file, output_file, sample_name, threads=1):
    """Count reads in peaks with library size normalization"""
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # First get total mapped reads for normalization
    logger.info("Calculating total mapped reads...")
    total_reads_cmd = f"samtools view -c -F 4 {bam_file}"
    total_reads = int(subprocess.check_output(total_reads_cmd, shell=True))
    logger.info(f"Total mapped reads: {total_reads}")
    
    # Create genome file from BAM header for proper chromosome ordering
    logger.info("Creating genome file from BAM header...")
    genome_file = f"{peaks_file}.genome"
    genome_cmd = f"samtools view -H {bam_file} | grep '^@SQ' | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > {genome_file}"
    subprocess.run(genome_cmd, shell=True, check=True)
    
    # Create sorted BED file
    logger.info("Sorting BED file...")
    sorted_peaks = f"{peaks_file}.sorted.tmp"
    sort_cmd = f"bedtools sort -g {genome_file} -i {peaks_file} > {sorted_peaks}"
    subprocess.run(sort_cmd, shell=True, check=True)
    
    try:
        # Count reads in peaks
        logger.info("Counting reads in peaks...")
        cmd = (f"bedtools coverage -a {sorted_peaks} -b {bam_file} "
               f"-sorted -g {genome_file} "
               f"-counts > {output_file}.tmp")
        
        subprocess.run(cmd, shell=True, check=True)
        
        # Read counts and normalize
        df = pd.read_csv(f"{output_file}.tmp", sep='\t', header=None,
                         names=['chr', 'start', 'end', 'gene', 'raw_count'])
        
        # Normalize to reads per million (RPM)
        df['count'] = df['raw_count'] * 1e6 / total_reads
        
        # Save normalized counts
        df.to_csv(output_file, sep='\t', index=False)
        
        logger.info(f"Normalized counts saved to {output_file}")
        logger.info(f"Mean raw count: {df['raw_count'].mean():.2f}")
        logger.info(f"Mean normalized count: {df['count'].mean():.2f}")
        
        # Quality control checks
        zero_peaks = (df['raw_count'] == 0).sum()
        if zero_peaks > len(df) * 0.5:
            logger.warning(f"More than 50% of peaks have zero reads: {zero_peaks}/{len(df)}")
        
        if total_reads < 1000000:
            logger.warning(f"Low number of mapped reads: {total_reads}")
            
    finally:
        # Cleanup temporary files
        for tmp_file in [f"{output_file}.tmp", genome_file, sorted_peaks]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

def main():
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
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    try:
        count_reads(args.peaks, args.bam, args.output, args.sample_name, threads=args.threads)
    except Exception as e:
        logger.error(f"Error processing {args.bam}: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main() 