#!/bin/bash
#SBATCH --job-name=SMARCB1_Analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=256GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.out"

set -e
set -x

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

export NCORES=32
export SAMTOOLS_THREADS=$((NCORES/2))

BG_SAMPLE="BG3"
BM_SAMPLE="BM3"

# Create gene promoter regions file (2kb upstream of TSS)
python scripts/create_promoter_regions.py \
    --gene-list scripts/gene_list.txt \
    --output results/gene_promoters.bed \
    --upstream 2000 \
    --downstream 500

# Count reads in promoter regions
for sample in ${BG_SAMPLE} ${BM_SAMPLE}; do
    python scripts/count_reads_in_peaks.py \
        --peaks results/gene_promoters.bed \
        --bam results/bowtie2_alt/${sample}.sorted.bam \
        --output results/peak_analysis_alt/${sample}_promoter_counts.txt \
        --threads $NCORES
done

# Compare SMARCB1 binding in promoters
python scripts/compare_peak_sizes.py \
    --peak-counts results/peak_analysis_alt/${BG_SAMPLE}_promoter_counts.txt \
                  results/peak_analysis_alt/${BM_SAMPLE}_promoter_counts.txt \
    --output results/peak_analysis_alt/promoter_comparison.txt \
    --threads $NCORES