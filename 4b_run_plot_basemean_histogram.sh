#!/bin/bash
#SBATCH --job-name=SMARCB1_make_expression_histogram
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="logs/make_expression_histogram.err"
#SBATCH --output="logs/make_expression_histogram.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python plot_basemean_histogram.py

echo "Expression histogram completed!" 