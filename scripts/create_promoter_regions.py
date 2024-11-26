import pandas as pd
import argparse
import logging
import gzip
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_gtf_line(line):
    """Parse a GTF line and return relevant information."""
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return None
    
    chrom, source, feature, start, end, score, strand, frame, attributes = fields
    
    # Parse attributes
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
    """Get gene coordinates from GTF file."""
    if gtf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parsed = parse_gtf_line(line)
            if parsed and parsed['feature'] == 'gene' and parsed['gene_name'] == gene:
                return parsed
    
    return None

def create_promoter_regions(gene_list_file, output_file, upstream=2000, downstream=500, 
                          gtf_file="/beegfs/datasets/genomes/mm10/annotation/gencode.vM25.annotation.gtf.gz"):
    """Create BED file of promoter regions for the specified genes."""
    
    # Check if GTF file exists
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Read gene list
    with open(gene_list_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Processing {len(genes)} genes using GTF file: {gtf_file}")
    
    # Get gene coordinates
    promoters = []
    for gene in genes:
        logger.info(f"Getting coordinates for {gene}")
        coords = get_gene_coordinates(gene, gtf_file)
        if coords:
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
    
    # Save to BED file
    if promoters:
        df = pd.DataFrame(promoters, columns=['chr', 'start', 'end', 'gene'])
        df.to_csv(output_file, sep='\t', index=False, header=False)
        logger.info(f"Created promoter regions for {len(promoters)} genes")
    else:
        raise Exception("No promoter regions could be created!")

def main():
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