#!/bin/bash
#SBATCH --job-name=8c_compare_bivalent_nonbivalent_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/8c_compare_bivalent_nonbivalent_heatmaps.err"
#SBATCH --output="logs/8c_compare_bivalent_nonbivalent_heatmaps.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/heatmaps_comparison_R

python 8c_compare_bivalent_nonbivalent_heatmaps.py

echo "8c_compare_bivalent_nonbivalent_heatmaps completed!" 