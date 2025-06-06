"""
This script annotates promoter comparison results with gene information and calculates differential binding statistics.

The script performs the following steps:
1. Reads promoter comparison data containing read counts and fold changes
2. Merges this data with gene information from a provided gene list
3. Calculates absolute log2 fold changes for sorting
4. Saves the annotated results to a file
5. Identifies significantly changed promoters (|log2FC| > 1)
6. Generates summary statistics on up/down-regulated promoters
7. Saves significant changes to a separate file

Input files:
- Promoter comparison file with read counts and fold changes
- Gene list file with coordinates and gene names

Output files:
- Annotated promoter results with gene information
- Filtered file containing only significant changes
"""

import pandas as pd
import argparse
import logging
import os
import numpy as np

# Set up logging configuration
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def read_gene_coordinates(gene_list_file, gtf_file="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/data/gencode.vM10.annotation.gtf.gz"):
    """
    Read gene coordinates from GTF file, filtered by genes in gene_list_file.
    
    Args:
        gene_list_file (str): Path to file containing list of gene names to filter for
        gtf_file (str): Path to GTF annotation file (gzipped)
        
    Returns:
        pandas.DataFrame: DataFrame containing filtered gene coordinates with columns:
            - chr: Chromosome name
            - start: Gene start position 
            - end: Gene end position
            - gene_name: Gene name/ID
    """
    logger.info(f"Reading gene list from {gene_list_file}")
    # Read list of genes to filter for
    with open(gene_list_file) as f:
        wanted_genes = set(line.strip() for line in f if line.strip())
    
    logger.info(f"Reading gene coordinates from {gtf_file}")
    # Read only gene entries from GTF
    gene_coords = []
    import gzip
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            # Only process gene features
            if fields[2] != 'gene':
                continue
            
            # Parse gene name from GTF attributes field
            attrs = dict(x.strip().split(' ', 1) for x in fields[8].rstrip(';').split('; '))
            gene_name = attrs.get('gene_name', '').strip('"')
            
            # Only keep genes in our wanted list
            if gene_name in wanted_genes:
                gene_coords.append({
                    'chr': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'gene_name': gene_name
                })
    
    gene_coords = pd.DataFrame(gene_coords)
    
    # Ensure chromosome names are standardized with 'chr' prefix
    gene_coords['chr'] = gene_coords['chr'].apply(
        lambda x: f"chr{x}" if not str(x).startswith('chr') else x
    )
    
    # Log summary statistics
    logger.info(f"Found coordinates for {len(gene_coords)} genes")
    logger.info("\nFirst few rows of gene coordinates:")
    logger.info(gene_coords.head())
    logger.info("\nChromosome formats in gene list:")
    logger.info(gene_coords['chr'].value_counts())
    
    return gene_coords

def find_overlapping_gene(row, gene_coords):
    """
    Find gene that overlaps with the given genomic coordinates.
    
    Args:
        row (pandas.Series): Series containing genomic coordinates with 'chr', 'start', 'end'
        gene_coords (pandas.DataFrame): DataFrame of gene coordinates to search
        
    Returns:
        str: Name of overlapping gene, or None if no overlap found
    """
    # Ensure chromosome format matches
    chr_name = row['chr'] if str(row['chr']).startswith('chr') else f"chr{row['chr']}"
    
    # Add debug logging for this region
    logger.debug(f"\nLooking for overlap with region: {chr_name}:{row['start']}-{row['end']}")
    
    # Check if chromosome exists in gene_coords
    if chr_name not in gene_coords['chr'].values:
        logger.debug(f"Chromosome {chr_name} not found in gene coordinates")
        return None
    
    # Get genes on this chromosome
    chr_genes = gene_coords[gene_coords['chr'] == chr_name]
    logger.debug(f"Found {len(chr_genes)} genes on chromosome {chr_name}")
    
    # Find genes that overlap with our region
    overlapping = chr_genes[
        (chr_genes['start'] <= row['end']) &
        (chr_genes['end'] >= row['start'])
    ]
    
    if len(overlapping) == 0:
        logger.debug(f"No overlapping genes found")
        # Log nearby genes for debugging
        nearby = chr_genes[
            (abs(chr_genes['start'] - row['start']) < 10000) |
            (abs(chr_genes['end'] - row['end']) < 10000)
        ]
        if not nearby.empty:
            logger.debug("Nearby genes (within 10kb):")
            logger.debug(nearby[['chr', 'start', 'end', 'gene_name']])
        return None
    elif len(overlapping) == 1:
        gene_name = overlapping.iloc[0]['gene_name']
        logger.debug(f"Found single overlapping gene: {gene_name}")
        return gene_name
    else:
        # If multiple genes overlap, take the one with the most overlap
        overlaps = []
        for _, gene in overlapping.iterrows():
            overlap_start = max(row['start'], gene['start'])
            overlap_end = min(row['end'], gene['end'])
            overlap_length = overlap_end - overlap_start
            overlaps.append((overlap_length, gene['gene_name']))
        gene_name = max(overlaps)[1]
        logger.debug(f"Found multiple overlapping genes ({len(overlapping)}), selected: {gene_name}")
        return gene_name

def annotate_promoters(comparison_file, gene_list_file, output_file, sample_name):
    """
    Annotate the promoter comparison results with gene information.
    
    Args:
        comparison_file (str): Path to input file with promoter comparison results
        gene_list_file (str): Path to file with gene list
        output_file (str): Path to save annotated results
        sample_name (str): Name of the sample for output files
        
    This function:
    1. Reads promoter comparison results
    2. Calculates missing log2 fold changes
    3. Matches promoter regions to genes
    4. Adds gene annotations and additional statistics
    5. Saves full results and significant changes
    """
    # Set logging to DEBUG for more detailed output
    logger.setLevel(logging.DEBUG)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read comparison results
    logger.info(f"Reading comparison results from {comparison_file}")
    comp_df = pd.read_csv(comparison_file, sep='\t')
    
    # Handle missing or invalid log2_fold_change values
    comp_df['log2_fold_change'] = pd.to_numeric(comp_df['log2_fold_change'], errors='coerce')
    
    # Calculate log2_fold_change for rows where it's missing
    # Add small pseudocount (0.01) to avoid division by zero
    mask = comp_df['log2_fold_change'].isna()
    comp_df.loc[mask, 'log2_fold_change'] = np.log2(
        (comp_df.loc[mask, 'bm_mean'] + 0.01) / (comp_df.loc[mask, 'bg_mean'] + 0.01)
    )
    
    logger.info(f"Processed {len(comp_df)} rows, including {mask.sum()} rows with recalculated fold changes")
    
    # Read gene coordinates
    gene_coords = read_gene_coordinates(gene_list_file)
    logger.info(f"Read {len(gene_coords)} gene coordinates")
    
    # Log example data for debugging
    logger.info("Example promoter region:")
    logger.info(comp_df.iloc[0].to_dict())
    logger.info("Example gene coordinates:")
    logger.info(gene_coords.iloc[0].to_dict())
    
    # Annotate results with gene names
    logger.info("Matching genomic coordinates to genes")
    annotated_df = comp_df.copy()
    
    # Find overlapping genes for each region
    annotated_df['gene'] = annotated_df.apply(
        lambda row: find_overlapping_gene(row, gene_coords), axis=1
    )
    
    # Handle unmatched regions
    unmatched = annotated_df['gene'].isna().sum()
    if unmatched > 0:
        logger.warning(f"{unmatched} regions could not be matched to genes")
        # Mark unmatched regions as "Unknown" instead of dropping them
        annotated_df['gene'] = annotated_df['gene'].fillna('Unknown')
    
    # Calculate additional statistics
    annotated_df['fold_change'] = 2 ** annotated_df['log2_fold_change']
    
    # Sort by absolute fold change
    annotated_df['abs_log2fc'] = abs(annotated_df['log2_fold_change'])
    annotated_df = annotated_df.sort_values('abs_log2fc', ascending=False)
    annotated_df = annotated_df.drop('abs_log2fc', axis=1)
    
    # Save results
    logger.info(f"Saving annotated results to {output_file}")
    annotated_df.to_csv(output_file, sep='\t', index=False)
    
    # Calculate and log summary statistics
    up_regulated = (annotated_df['log2_fold_change'] > 1).sum()
    down_regulated = (annotated_df['log2_fold_change'] < -1).sum()
    
    logger.info(f"Summary:")
    logger.info(f"  Total promoters analyzed: {len(annotated_df)}")
    logger.info(f"  Up-regulated in BM (log2FC > 1): {up_regulated}")
    logger.info(f"  Down-regulated in BM (log2FC < -1): {down_regulated}")
    
    # Save significant changes to separate file
    output_dir = os.path.dirname(output_file)
    top_changes = annotated_df[abs(annotated_df['log2_fold_change']) > 1].copy()
    top_changes.to_csv(os.path.join(output_dir, f'significant_changes_{sample_name}.tsv'), 
                      sep='\t', index=False)
    
    logger.info(f"Saved {len(top_changes)} significant changes to significant_changes_{sample_name}.tsv")

def main():
    """
    Main function to parse command line arguments and run promoter annotation.
    """
    parser = argparse.ArgumentParser(description='Annotate promoter comparison results')
    parser.add_argument('--input', required=True,
                        help='Promoter comparison file')
    parser.add_argument('--gene-list', required=True,
                        help='Gene coordinates file (chr, start, end, gene_name)')
    parser.add_argument('--output', required=True,
                        help='Output annotated file')
    parser.add_argument('--sample-name', required=True,
                        help='Sample name') 
    
    args = parser.parse_args()
    
    try:
        annotate_promoters(args.input, args.gene_list, args.output, args.sample_name)
    except Exception as e:
        logger.error(f"Annotation failed: {str(e)}")
        raise

if __name__ == '__main__':
    main() 


