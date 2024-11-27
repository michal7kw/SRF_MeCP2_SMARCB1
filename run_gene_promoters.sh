#!/bin/bash
#SBATCH --job-name=PromoterRegions
#SBATCH --account=kubacki.michal
#SBATCH --mem=256GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/gene_promoters.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/gene_promoters.out"

set -e
set -x

# Load required environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set parameters
export NCORES=32
export SAMTOOLS_THREADS=$((NCORES/2))

BG_SAMPLE="BG1"
BM_SAMPLE="BM1"

# Step 2: Create promoter regions
echo "Creating promoter regions..."
python scripts/create_promoter_regions.py \
    --gene-list scripts/gene_list.txt \
    --output results/gene_promoters.bed \
    --upstream 2000 \
    --downstream 500 \
    --gtf /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/gencode.vM10.annotation.gtf.gz

