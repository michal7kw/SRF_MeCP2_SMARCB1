import argparse
import subprocess
import os
import sys
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor
import pysam

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def count_reads_chunk(chunk, bam_file):
    """Count reads in a chunk of peaks"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    
    for _, peak in chunk.iterrows():
        count = bam.count(peak['chr'], peak['start'], peak['end'])
        results.append([peak['chr'], peak['start'], peak['end'], peak['gene'], count])
    
    bam.close()
    return results

def count_reads(peaks_file, bam_file, output_file, sample_name, temp_dir=None, threads=1):
    """Count reads in peaks with parallel processing"""
    
    # Read peaks
    peaks_df = pd.read_csv(peaks_file, sep='\t', header=None,
                          names=['chr', 'start', 'end', 'gene'])
    
    # Split into chunks for parallel processing
    chunk_size = max(1, len(peaks_df) // (threads * 2))
    chunks = [peaks_df[i:i + chunk_size] for i in range(0, len(peaks_df), chunk_size)]
    
    # Process chunks in parallel
    all_results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(count_reads_chunk, chunk, bam_file)
            for chunk in chunks
        ]
        for future in futures:
            all_results.extend(future.result())
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results, 
                            columns=['chr', 'start', 'end', 'gene', 'raw_count'])
    
    # Save results
    if temp_dir:
        temp_file = os.path.join(temp_dir, f"temp_{sample_name}_counts.txt")
        results_df.to_csv(temp_file, sep='\t', index=False)
        os.system(f"sort -k1,1 -k2,2n {temp_file} > {output_file}")
        os.remove(temp_file)
    else:
        results_df.to_csv(output_file, sep='\t', index=False)

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
    
    try:
        count_reads(args.peaks, args.bam, args.output, args.sample_name, threads=args.threads)
    except Exception as e:
        logger.error(f"Error processing {args.bam}: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main() 