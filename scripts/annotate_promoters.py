import pandas as pd
import argparse
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def annotate_promoters(comparison_file, gene_list_file, output_file, sample_name):
    """
    Annotate the promoter comparison results with gene information.
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read comparison results
    logger.info(f"Reading comparison results from {comparison_file}")
    comp_df = pd.read_csv(comparison_file, sep='\t')
    
    # Read gene list
    logger.info(f"Reading gene list from {gene_list_file}")
    with open(gene_list_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    # Create gene information dictionary
    gene_info = pd.DataFrame({
        'gene': genes,
        'index': range(len(genes))
    })
    
    # Merge comparison results with gene information
    logger.info("Annotating results")
    annotated_df = comp_df.copy()
    
    # Add gene names based on the order in the gene list
    # (assuming the order matches between the BED file and gene list)
    if len(annotated_df) != len(genes):
        logger.warning(f"Number of peaks ({len(annotated_df)}) doesn't match number of genes ({len(genes)})")
    
    annotated_df['gene'] = genes[:len(annotated_df)]
    
    # Calculate additional statistics
    annotated_df['fold_change'] = 2 ** annotated_df['log2_fold_change']
    
    # Sort by absolute fold change
    annotated_df['abs_log2fc'] = abs(annotated_df['log2_fold_change'])
    annotated_df = annotated_df.sort_values('abs_log2fc', ascending=False)
    annotated_df = annotated_df.drop('abs_log2fc', axis=1)
    
    # Save results
    logger.info(f"Saving annotated results to {output_file}")
    annotated_df.to_csv(output_file, sep='\t', index=False)
    
    # Log summary statistics
    up_regulated = (annotated_df['log2_fold_change'] > 1).sum()
    down_regulated = (annotated_df['log2_fold_change'] < -1).sum()
    
    logger.info(f"Summary:")
    logger.info(f"  Total promoters analyzed: {len(annotated_df)}")
    logger.info(f"  Up-regulated in BM (log2FC > 1): {up_regulated}")
    logger.info(f"  Down-regulated in BM (log2FC < -1): {down_regulated}")
    
    # Save top changes to a separate file
    output_dir = os.path.dirname(output_file)
    top_changes = annotated_df[abs(annotated_df['log2_fold_change']) > 1].copy()
    top_changes.to_csv(os.path.join(output_dir, f'significant_changes_{sample_name}.tsv'), 
                      sep='\t', index=False)
    
    logger.info(f"Saved {len(top_changes)} significant changes to significant_changes_{sample_name}.tsv")

def main():
    parser = argparse.ArgumentParser(description='Annotate promoter comparison results')
    parser.add_argument('--input', required=True,
                        help='Promoter comparison file')
    parser.add_argument('--gene-list', required=True,
                        help='Original gene list file')
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