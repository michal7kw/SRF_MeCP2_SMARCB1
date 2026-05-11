#!/bin/bash
#SBATCH --job-name=4_extract_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="logs/4_extract_genes.err"
#SBATCH --output="logs/4_extract_genes.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python 4_extract_genes.py

echo "4_extract_genes completed!" 