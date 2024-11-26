#!/bin/bash
#SBATCH --job-name=SMARCB1_Annotate
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/annotate.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/annotate.out"

set -e
set -x

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Annotate promoter regions with SMARCB1 binding changes
python scripts/annotate_promoters.py \
    --input results/peak_analysis_alt/promoter_comparison.txt \
    --gene-list scripts/gene_list.txt \
    --output results/peak_analysis_alt/annotated_promoters.txt