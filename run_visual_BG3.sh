#!/bin/bash
#SBATCH --job-name=SMARCB1_Visual
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/visual.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/visual.out"

set -e
set -x

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

python scripts/visualize_promoters.py \
    --input results/peak_analysis_alt/annotated_promoters.txt \
    --output-dir results/peak_analysis_alt/plots \
    --sample-name BG3
    