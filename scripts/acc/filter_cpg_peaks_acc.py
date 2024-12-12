import argparse
import subprocess
import os
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import pybedtools
import tempfile

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_chunk(chunk, cpg_islands):
    """Process a chunk of peaks in parallel"""
    return pybedtools.BedTool(chunk.values.tolist()).intersect(cpg_islands, u=True)

def filter_cpg_peaks(peaks_file, cpg_islands_file, output_file, temp_dir=None, threads=1):
    """Filter peaks to only those overlapping with CpG islands."""
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Set temporary directory
    if temp_dir:
        pybedtools.set_tempdir(temp_dir)
    
    # Read peaks and CpG islands
    peaks_df = pd.read_csv(peaks_file, sep='\t', header=None)
    cpg_islands = pybedtools.BedTool(cpg_islands_file)
    
    # Split data into chunks for parallel processing
    chunk_size = max(1, len(peaks_df) // (threads * 2))
    chunks = [peaks_df[i:i + chunk_size] for i in range(0, len(peaks_df), chunk_size)]
    
    # Process chunks in parallel
    filtered_peaks = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(process_chunk, chunk, cpg_islands)
            for chunk in chunks
        ]
        for future in futures:
            result = future.result()
            filtered_peaks.extend([line for line in result])
    
    # Write results
    with open(output_file, 'w') as f:
        for peak in filtered_peaks:
            f.write(str(peak))
    
    # Cleanup
    pybedtools.cleanup()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--peaks', required=True)
    parser.add_argument('--cpg-islands', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--temp-dir', help='Directory for temporary files')
    parser.add_argument('--threads', type=int, default=1)
    
    args = parser.parse_args()
    
    filter_cpg_peaks(
        args.peaks, 
        args.cpg_islands, 
        args.output,
        args.temp_dir,
        args.threads
    )

if __name__ == '__main__':
    main() 