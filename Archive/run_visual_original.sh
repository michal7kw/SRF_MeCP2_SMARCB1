#!/bin/bash
#SBATCH --job-name=visual
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/visual.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/visual.out"

set -e  # Exit on error
set -x  # Print commands being executed

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

python scripts/visualize_peaks.py \
    --input results/peak_analysis_alt/peak_size_comparison.txt \
    --output-dir results/peak_analysis_alt/plots

    