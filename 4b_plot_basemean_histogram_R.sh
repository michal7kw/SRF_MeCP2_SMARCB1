#!/bin/bash
#SBATCH --job-name=4b_plot_basemean_histogram
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="logs/4b_plot_basemean_histogram.err"
#SBATCH --output="logs/4b_plot_basemean_histogram.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python 4b_plot_basemean_histogram.py

echo "4b_plot_basemean_histogram completed!" 