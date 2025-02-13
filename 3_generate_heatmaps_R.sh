#!/bin/bash
#SBATCH --job-name=3_generate_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/3_generate_heatmaps.err"
#SBATCH --output="logs/3_generate_heatmaps.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/heatmaps

Rscript 3_run_generate_heatmaps.R

echo "3_generate_heatmaps completed!" 