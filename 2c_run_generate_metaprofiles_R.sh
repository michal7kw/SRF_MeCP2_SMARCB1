#!/bin/bash
#SBATCH --job-name=generate_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/generate_metaprofiles.err"
#SBATCH --output="logs/generate_metaprofiles.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/metaprofiles_R

Rscript generate_metaprofiles.R

echo "generate_metaprofiles completed!" 