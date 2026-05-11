"""
This script creates promoter regions in BED format for a list of genes.

Key features:
- Reads gene list from input file
- Extracts gene coordinates from GTF annotation file
- Creates promoter regions with configurable upstream/downstream distances
- Outputs promoter regions in BED format
- Handles both gzipped and uncompressed GTF files
- Provides detailed logging

Input:
- Gene list file with one gene name per line
- GTF annotation file (gzipped or uncompressed)
- Upstream and downstream distances for promoter regions
- Output file path for BED format promoters

Output:
- BED file containing promoter regions with columns:
  chr, start, end, gene_name
- Logging information about processing
"""

import pandas as pd
import argparse
import logging
import gzip
import os

# Set up logging configuration
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_gtf_line(line):
    """
    Parse a single line from a GTF file and extract relevant fields.
    
    Args:
        line (str): A line from GTF file
        
    Returns:
        dict: Dictionary containing parsed fields (chrom, start, end, strand, gene_name, feature)
             Returns None if line cannot be parsed
    """
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return None
    
    chrom, source, feature, start, end, score, strand, frame, attributes = fields
    
    # Parse GTF attributes field into dictionary
    attr_dict = {}
    for attr in attributes.split(';'):
        if attr.strip():
            key, value = attr.strip().split(' ', 1)
            attr_dict[key] = value.strip('"')
    
    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end),
        'strand': strand,
        'gene_name': attr_dict.get('gene_name', ''),
        'feature': feature
    }

def get_gene_coordinates(gene, gtf_file):
    """
    Extract coordinates for a specific gene from GTF file.
    
    Args:
        gene (str): Gene name to search for
        gtf_file (str): Path to GTF annotation file
        
    Returns:
        dict: Gene coordinates and metadata if found, None otherwise
    """
    # Handle gzipped or regular files
    if gtf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            
            parsed = parse_gtf_line(line)
            if parsed and parsed['feature'] == 'gene' and parsed['gene_name'] == gene:
                return parsed
    
    return None

def create_promoter_regions(gene_list_file, output_file, upstream=2000, downstream=500, 
                          gtf_file="/beegfs/datasets/genomes/mm10/annotation/gencode.vM25.annotation.gtf.gz"):
    """
    Create BED file containing promoter regions for specified genes.
    
    Args:
        gene_list_file (str): Path to file containing gene names
        output_file (str): Path to output BED file
        upstream (int): Base pairs upstream of TSS to include
        downstream (int): Base pairs downstream of TSS to include
        gtf_file (str): Path to GTF annotation file
        
    Raises:
        FileNotFoundError: If GTF file not found
        Exception: If no promoter regions could be created
    """
    
    # Validate GTF file exists
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Read gene list from file
    with open(gene_list_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Processing {len(genes)} genes using GTF file: {gtf_file}")
    
    # Process each gene to get promoter coordinates
    promoters = []
    for gene in genes:
        logger.info(f"Getting coordinates for {gene}")
        coords = get_gene_coordinates(gene, gtf_file)
        if coords:
            # Calculate promoter region based on strand
            if coords['strand'] == '+':
                prom_start = max(0, coords['start'] - upstream)
                prom_end = coords['start'] + downstream
            else:
                prom_start = max(0, coords['end'] - downstream)
                prom_end = coords['end'] + upstream
            
            promoters.append([
                coords['chrom'],
                prom_start,
                prom_end,
                gene
            ])
        else:
            logger.warning(f"Could not find coordinates for {gene}")
    
    # Write promoter regions to BED file
    if promoters:
        df = pd.DataFrame(promoters, columns=['chr', 'start', 'end', 'gene'])
        df.to_csv(output_file, sep='\t', index=False, header=False)
        logger.info(f"Created promoter regions for {len(promoters)} genes")
    else:
        raise Exception("No promoter regions could be created!")

def main():
    """Parse command line arguments and run promoter region creation."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-list', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--upstream', type=int, default=2000)
    parser.add_argument('--downstream', type=int, default=500)
    parser.add_argument('--gtf', default="/beegfs/datasets/genomes/mm10/annotation/gencode.vM25.annotation.gtf.gz",
                      help="Path to GTF file")
    args = parser.parse_args()
    
    create_promoter_regions(args.gene_list, args.output, args.upstream, args.downstream, args.gtf)

if __name__ == '__main__':
    main()

"""
Summary:
This script takes a list of genes and creates a BED file containing their promoter regions.
It uses a GTF annotation file to look up the genomic coordinates and strand information
for each gene. For genes on the + strand, the promoter is defined as a region upstream
and downstream of the transcription start site (TSS). For genes on the - strand, the
promoter is defined relative to the transcription end site (TES). The script handles
both gzipped and uncompressed GTF files, provides detailed logging, and includes error
handling for missing files or failed promoter creation.

The output BED file can be used for downstream analyses like:
- Analyzing ChIP-seq peaks in promoter regions
- Calculating promoter methylation levels
- Identifying transcription factor binding sites
- Studying regulatory elements near genes of interest
"""