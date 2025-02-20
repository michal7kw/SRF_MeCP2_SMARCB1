#!/bin/bash
#SBATCH --job-name=process_gene_lists_Py
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/process_gene_lists_Py.err"
#SBATCH --output="logs/process_gene_lists_Py.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python process_gene_lists.py

echo "process_gene_lists.py completed!" 