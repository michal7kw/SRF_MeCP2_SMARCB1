#!/bin/bash
#SBATCH --job-name=SMARCB1_Analysis_BG1_BM2
#SBATCH --account=kubacki.michal
#SBATCH --mem=256GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/smarcb1_BG1_BM2.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/smarcb1_BG1_BM2.out"

set -e
set -x

# Load required environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set parameters
export NCORES=32
export SAMTOOLS_THREADS=$((NCORES/2))

BG_SAMPLE="BG1"
BM_SAMPLE="BM2"

# # Step 2: Create promoter regions
# echo "Creating promoter regions..."
# python scripts/create_promoter_regions.py \
#     --gene-list scripts/gene_list.txt \
#     --output results/gene_promoters.bed \
#     --upstream 2000 \
#     --downstream 500 \
#     --gtf /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/gencode.vM25.annotation.gtf.gz

# Step 3: Count reads in promoter regions
echo "Counting reads in promoter regions..."
echo "Processing ${BG_SAMPLE}..."
python scripts/count_reads_in_peaks.py \
    --peaks results/gene_promoters.bed \
    --bam results/bowtie2_alt/${BG_SAMPLE}.sorted.bam \
    --output results/peak_analysis_alt/${BG_SAMPLE}_promoter_counts.txt \
    --sample-name ${BG_SAMPLE} \
    --threads $NCORES

echo "Counting reads in promoter regions..."
echo "Processing ${BM_SAMPLE}..."
python scripts/count_reads_in_peaks.py \
    --peaks results/gene_promoters.bed \
    --bam results/bowtie2_alt/${BM_SAMPLE}.sorted.bam \
    --output results/peak_analysis_alt/${BM_SAMPLE}_promoter_counts.txt \
    --sample-name ${BM_SAMPLE} \
    --threads $NCORES

# Step 4: Compare SMARCB1 binding
echo "Comparing SMARCB1 binding..."
python scripts/compare_peak_sizes.py \
    --peak-counts results/peak_analysis_alt/${BG_SAMPLE}_promoter_counts.txt \
                  results/peak_analysis_alt/${BM_SAMPLE}_promoter_counts.txt \
    --output results/peak_analysis_alt/promoter_comparison_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE} \
    --threads $NCORES

# Step 5: Annotate results
echo "Annotating results..."
python scripts/annotate_promoters.py \
    --input results/peak_analysis_alt/promoter_comparison_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --gene-list scripts/gene_list.txt \
    --output results/peak_analysis_alt/annotated_promoters_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE}

# Step 6: Create visualizations
echo "Creating visualizations..."
python scripts/visualize_promoters.py \
    --input results/peak_analysis_alt/annotated_promoters_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --output-dir results/peak_analysis_alt/plots \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE}

echo "Analysis complete. Check results in results/peak_analysis_alt/" 