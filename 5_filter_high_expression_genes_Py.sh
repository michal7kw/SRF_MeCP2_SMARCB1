#!/bin/bash
#SBATCH --job-name=5_filter_high_expression_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="logs/5_filter_high_expression_genes.err"
#SBATCH --output="logs/5_filter_high_expression_genes.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Set the baseMean threshold
THRESHOLD1=100.0
THRESHOLD2=1000.0

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python 5_filter_high_expression_genes.py --threshold ${THRESHOLD1}
python 5_filter_high_expression_genes.py --threshold ${THRESHOLD2}

echo "5_filter_high_expression_genes completed!" 