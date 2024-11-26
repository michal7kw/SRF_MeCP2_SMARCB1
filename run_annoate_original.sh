#!/bin/bash
#SBATCH --job-name=Annotate
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/annotate.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/annotate.out"

set -e  # Exit on error
set -x  # Print commands being executed

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set number of cores for parallel processing
export NCORES=32

# Enable parallel processing for samtools
export SAMTOOLS_THREADS=$((NCORES/2))

BG_SAMPLE="BG3"
BM_SAMPLE="BM3"

# Prepare significant peaks
# python scripts/prepare_significant_peaks.py \
#     --input results/peak_analysis_alt/peak_size_comparison.txt \
#     --output results/peak_analysis_alt/significant_peaks.txt \
#     --fc-threshold 1.0

# Annotate significant peaks with nearby genes
python scripts/annotate_peaks.py \
    --input results/peak_analysis_alt/significant_peaks.txt \
    --output results/peak_analysis_alt/annotated_peaks.txt \
    --window 5000