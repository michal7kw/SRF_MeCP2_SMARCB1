#!/bin/bash
#SBATCH --job-name=6_find_expressed_bivalent_targets
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="logs/6_find_expressed_bivalent_targets.err"
#SBATCH --output="logs/6_find_expressed_bivalent_targets.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Set the target file
TARGET_FILE="Gene_lists/targets/high_expression_targets2_1000.0.csv"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python 6_find_expressed_bivalent_targets.py --target_file ${TARGET_FILE}

echo "6_find_expressed_bivalent_targets completed!" 