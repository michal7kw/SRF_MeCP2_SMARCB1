#!/bin/bash
#SBATCH --job-name=8b_compare_bivalent_nonbivalent_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/8b_compare_bivalent_nonbivalent_metaprofiles.err"
#SBATCH --output="logs/8b_compare_bivalent_nonbivalent_metaprofiles.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/metaprofiles_comparison_R

Rscript 8b_compare_bivalent_nonbivalent_metaprofiles.R

echo "8b_compare_bivalent_nonbivalent_metaprofiles completed!" 