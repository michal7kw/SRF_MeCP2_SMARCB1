#!/bin/bash
#SBATCH --job-name=8e_combined_metaprofile_comparisons_R
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/8e_combined_metaprofile_comparisons_R.err"
#SBATCH --output="logs/8e_combined_metaprofile_comparisons_R.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for plots
mkdir -p results/metaprofiles_comparison_R

Rscript 8e_combined_metaprofile_comparisons.R

echo "8e_combined_metaprofile_comparisons_R completed!" 