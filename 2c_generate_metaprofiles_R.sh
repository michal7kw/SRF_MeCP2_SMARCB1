#!/bin/bash
#SBATCH --job-name=2c_generate_metaprofiles_R
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/2c_generate_metaprofiles_R.err"
#SBATCH --output="logs/2c_generate_metaprofiles_R.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory for heatmaps
mkdir -p results/metaprofiles_R

Rscript 2c_run_generate_metaprofiles.R

echo "2c_generate_metaprofiles_R completed!" 