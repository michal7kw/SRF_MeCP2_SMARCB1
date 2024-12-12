#!/bin/bash
#SBATCH --job-name=SMARCB1_Analysis_BG1_BM1
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="logs/smarcb1_BG1_BM1.err"
#SBATCH --output="logs/smarcb1_BG1_BM1.out"

set -e
set -x

# Load required environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

# Set parameters with optimized thread allocation
export NCORES=32
export SAMTOOLS_THREADS=$((NCORES/2))
export PYTHONUNBUFFERED=1
export TMPDIR="/tmp"

BG_SAMPLE="BG1"
BM_SAMPLE="BM1"

# Add variables for file paths
GENE_LIST_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/scripts/genes_lists/iterative_alternative_results_alternative"
GENE_LIST="${GENE_LIST_DIR}/genes_width_weighted_NSC.csv"

RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/analysis/results_iara"
PROMOTERS_BED="${RESULTS_DIR}/gene_promoters.bed"
PROMOTERS_CPG_BED="${RESULTS_DIR}/gene_promoters_cpg.bed"

DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/data"
ANALYSIS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/results"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta

mkdir -p $RESULTS_DIR
mkdir -p $TMPDIR/smarcb1_analysis_${BG_SAMPLE}_${BM_SAMPLE}
TEMP_DIR="$TMPDIR/smarcb1_analysis_${BG_SAMPLE}_${BM_SAMPLE}"

# Run promoter region creation and CpG filtering in parallel
echo "Creating promoter regions and filtering CpG islands in parallel..."
(
    python scripts/create_promoter_regions.py \
        --gene-list $GENE_LIST \
        --output $PROMOTERS_BED \
        --upstream 2000 \
        --downstream 500 \
        --gtf $DATA_DIR/gencode.vM10.annotation.gtf.gz \
        --temp-dir $TEMP_DIR \
        --threads $((NCORES/2)) &

    python scripts/filter_cpg_peaks.py \
        --peaks $PROMOTERS_BED \
        --cpg-islands $DATA_DIR/cpg_islands.bed \
        --output $PROMOTERS_CPG_BED \
        --temp-dir $TEMP_DIR \
        --threads $((NCORES/2)) &

    wait
)

# Process both samples in parallel
echo "Processing samples in parallel..."
(
    echo "Processing ${BG_SAMPLE}..."
    python scripts/count_reads_in_peaks.py \
        --peaks $PROMOTERS_CPG_BED \
        --bam ${ANALYSIS_DIR}/bowtie2_alt/${BG_SAMPLE}.sorted.bam \
        --output ${RESULTS_DIR}/peak_analysis_alt/${BG_SAMPLE}_promoter_counts.txt \
        --sample-name ${BG_SAMPLE} \
        --temp-dir $TEMP_DIR \
        --threads $((NCORES/2)) &

    echo "Processing ${BM_SAMPLE}..."
    python scripts/count_reads_in_peaks.py \
        --peaks $PROMOTERS_CPG_BED \
        --bam ${ANALYSIS_DIR}/bowtie2_alt/${BM_SAMPLE}.sorted.bam \
        --output ${RESULTS_DIR}/peak_analysis_alt/${BM_SAMPLE}_promoter_counts.txt \
        --sample-name ${BM_SAMPLE} \
        --temp-dir $TEMP_DIR \
        --threads $((NCORES/2)) &

    wait
)

# Compare and annotate results
echo "Comparing and annotating results..."
python scripts/compare_peak_sizes.py \
    --peak-counts ${RESULTS_DIR}/peak_analysis_alt/${BG_SAMPLE}_promoter_counts.txt \
                  ${RESULTS_DIR}/peak_analysis_alt/${BM_SAMPLE}_promoter_counts.txt \
    --output ${RESULTS_DIR}/peak_analysis_alt/promoter_comparison_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE} \
    --temp-dir $TEMP_DIR \
    --threads $NCORES

python scripts/annotate_promoters.py \
    --input ${RESULTS_DIR}/peak_analysis_alt/promoter_comparison_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --gene-list $GENE_LIST \
    --output ${RESULTS_DIR}/peak_analysis_alt/annotated_promoters_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE} \
    --temp-dir $TEMP_DIR \
    --threads $NCORES

# Create visualizations (can be parallelized if multiple plots)
echo "Creating visualizations..."
python scripts/visualize_promoters.py \
    --input ${RESULTS_DIR}/peak_analysis_alt/annotated_promoters_${BG_SAMPLE}_${BM_SAMPLE}.txt \
    --output-dir ${RESULTS_DIR}/peak_analysis_alt/plots \
    --sample-name ${BG_SAMPLE}_${BM_SAMPLE} \
    --temp-dir $TEMP_DIR \
    --threads $NCORES

# Cleanup temporary directory
rm -rf $TEMP_DIR

echo "Analysis complete. Check results in results/peak_analysis_alt/" 