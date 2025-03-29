#!/bin/bash
#SBATCH --job-name=1b_resume_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=1-4
#SBATCH --error="logs/1b_resume_pipeline.err"
#SBATCH --output="logs/1b_resume_pipeline.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define sample names array
declare -a SAMPLES=(
    "BG1"
    "BG2"
    "BG3"
    "BM3"
)

# Get current sample name from array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing sample: $SAMPLE"

# Define directories
RESULTS_DIR="results"

echo "=== Resuming pipeline for ${SAMPLE} ==="

# 5. Add read groups to header
echo "Adding read groups to header..."
# Create temporary header file
samtools view -H ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam > temp_header.sam
# Add read group line if not present
if ! grep -q "^@RG" temp_header.sam; then
    echo -e "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA" >> temp_header.sam
fi
# Apply new header and add RG tags to reads
samtools reheader temp_header.sam ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam | \
    samtools view -h | \
    awk -v sample=$SAMPLE 'BEGIN {OFS="\t"} /^@/ {print} !/^@/ {print $0 "\tRG:Z:"sample}' | \
    samtools view -b > ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
rm temp_header.sam

# Index the BAM with read groups
samtools index ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam

# Remove original sorted BAM to save space
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam

# 6. Remove duplicates
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam \
    OUTPUT=${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/${SAMPLE}.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

# Remove intermediate BAM file with read groups
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam.bai

# 7. Generate bigWig files
echo "Generating bigWig files..."
# RPKM normalized
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${RESULTS_DIR}/bigwig/${SAMPLE}_RPKM.bw \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 32 \
    --extendReads \
    --ignoreDuplicates

# CPM normalized
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${RESULTS_DIR}/bigwig/${SAMPLE}_CPM.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 32 \
    --extendReads \
    --ignoreDuplicates

# 8. Call peaks
echo "Calling peaks..."
macs2 callpeak \
    -t ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    -f BAMPE \
    -g mm \
    -n ${SAMPLE} \
    --outdir ${RESULTS_DIR}/peaks \
    --nomodel \
    --extsize 200 \
    --keep-dup all \
    --qvalue 0.05 \
    --call-summits

# 9. Generate QC metrics
echo "Generating QC metrics..."
{
    echo "=== QC Metrics for ${SAMPLE} ==="
    echo -e "\nPost-deduplication stats:"
    samtools flagstat ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam
    echo -e "\nNumber of peaks:"
    wc -l ${RESULTS_DIR}/peaks/${SAMPLE}_peaks.narrowPeak
} > ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

echo "=== Pipeline resumed and completed successfully for ${SAMPLE}! ===" 