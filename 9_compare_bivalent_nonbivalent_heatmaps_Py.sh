#!/bin/bash
#SBATCH --job-name=9_compare_bivalent_nonbivalent_heatmaps_Py
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/9_compare_bivalent_nonbivalent_heatmaps_Py.err"
#SBATCH --output="logs/9_compare_bivalent_nonbivalent_heatmaps_Py.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/metaprofiles_comparison_R

python 9_compare_bivalent_nonbivalent_heatmaps.py

echo "9_compare_bivalent_nonbivalent_heatmaps.py completed!" 