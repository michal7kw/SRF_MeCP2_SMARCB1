import pandas as pd
import argparse
import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def download_gene_annotations():
    """Download gene annotations if not present."""
    if not os.path.exists("genes.bed"):
        logger.info("Downloading gene annotations...")
        # Download from UCSC
        cmd = """
        wget -O - "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz" | \
        gunzip -c | \
        awk 'BEGIN{OFS="\t"} {print $3,$5,$6,$13,".",".",$4}' > genes.bed
        """
        subprocess.run(cmd, shell=True, check=True)
        
        # Sort the file
        subprocess.run("sort -k1,1 -k2,2n genes.bed > genes.sorted.bed", 
                      shell=True, check=True)
    return "genes.sorted.bed"

def annotate_peaks(peaks_file, output_file, window_size=5000):
    """
    Annotate peaks with nearby genes.
    window_size: base pairs upstream/downstream to look for genes
    """
    # Ensure we have gene annotations
    genes_file = download_gene_annotations()
    
    # Create temporary sorted peaks file with proper BED format
    sorted_peaks = peaks_file + ".sorted"
    
    # Convert peaks file to proper BED format while preserving original data
    subprocess.run(f"""awk 'BEGIN{{OFS="\t"}} 
                       NR>1 {{print $1,$2,$3,$4";"$5";"$6,0,"."}}' {peaks_file} | \
                    sort -k1,1 -k2,2n > {sorted_peaks}""", 
                  shell=True, check=True)
    
    try:
        # Find nearest genes using bedtools
        cmd = f"""
        bedtools window -a {sorted_peaks} -b {genes_file} \
        -w {window_size} > {output_file}.tmp
        """
        subprocess.run(cmd, shell=True, check=True)
        
        # Read and process the results
        # Split the combined field back into original columns
        df = pd.read_csv(f"{output_file}.tmp", sep='\t', header=None,
                        names=['chr', 'start', 'end', 'combined_data', 'score', 'strand',
                              'gene_chr', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand', 'gene_id'])
        
        # Split the combined data back into original columns
        df[['bg_mean', 'bm_mean', 'log2_fold_change']] = df['combined_data'].str.split(';', expand=True)
        
        # Convert back to numeric
        df['bg_mean'] = pd.to_numeric(df['bg_mean'])
        df['bm_mean'] = pd.to_numeric(df['bm_mean'])
        df['log2_fold_change'] = pd.to_numeric(df['log2_fold_change'])
        
        # Group by peak and aggregate gene information
        peak_annotations = df.groupby(['chr', 'start', 'end', 'bg_mean', 'bm_mean', 'log2_fold_change']).agg({
            'gene_name': lambda x: ','.join(sorted(set(x))),
            'gene_id': lambda x: ','.join(sorted(set(x)))
        }).reset_index()
        
        # Add distance to TSS
        peak_annotations['peak_center'] = (peak_annotations['start'] + peak_annotations['end']) // 2
        
        # Save results
        peak_annotations.to_csv(output_file, sep='\t', index=False)
        
        # Generate summary
        logger.info(f"Total peaks annotated: {len(peak_annotations)}")
        logger.info(f"Peaks with gene associations: {(peak_annotations['gene_name'] != '').sum()}")
        
    finally:
        # Cleanup
        if os.path.exists(sorted_peaks):
            os.remove(sorted_peaks)
        if os.path.exists(f"{output_file}.tmp"):
            os.remove(f"{output_file}.tmp")

def main():
    parser = argparse.ArgumentParser(description='Annotate peaks with nearby genes')
    parser.add_argument('--input', required=True,
                        help='Significant peaks file (from visualize_peaks.py)')
    parser.add_argument('--output', required=True,
                        help='Output annotated peaks file')
    parser.add_argument('--window', type=int, default=5000,
                        help='Window size (bp) to search for nearby genes')
    
    args = parser.parse_args()
    
    try:
        annotate_peaks(args.input, args.output, args.window)
    except Exception as e:
        logger.error(f"Annotation failed: {str(e)}")
        raise

if __name__ == '__main__':
    main() 