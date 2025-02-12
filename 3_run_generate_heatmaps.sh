#!/bin/bash
#SBATCH --job-name=SMARCB1_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/heatmaps.err"
#SBATCH --output="logs/heatmaps.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/heatmaps

Rscript plot_binding_heatmaps.R

echo "Heatmap generation completed!" 