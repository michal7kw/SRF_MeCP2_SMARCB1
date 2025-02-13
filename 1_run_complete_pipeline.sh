#!/bin/bash

# This script runs a complete ChIP-seq analysis pipeline for SMARCB1 samples
# It performs the following steps:
# 1. Initial FastQC quality control
# 2. Adapter and quality trimming with Trimmomatic
# 3. Alignment with Bowtie2
# 4. Read group addition
# 5. PCR duplicate removal
# 6. BigWig file generation
# 7. Peak calling with MACS2

#SBATCH --job-name=1_run_complete_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=1-4                          # Process 4 samples in parallel
#SBATCH --error="logs/1_run_complete_pipelin.err"
#SBATCH --output="logs/1_run_complete_pipeline.out"

# Set working directory and paths
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
FASTQ_DIR="90-1102945428/00_fastq"
RESULTS_DIR="results"
cd ${WORKING_DIR}

# Load required software environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define sample names to process
declare -a SAMPLES=(
    "BG1"
    "BG2"
    "BG3"
    "BM3"
)

# Get current sample based on array task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing sample: $SAMPLE"

# Create output directory structure
mkdir -p ${RESULTS_DIR}/{trimmed,filtered,bowtie2,peaks,bigwig,fastqc}

# Generate chromosome sizes file for mm10 genome if not present
CHROM_SIZES="${WORKING_DIR}/mm10/mm10.chrom.sizes"
if [ ! -f "$CHROM_SIZES" ]; then
    echo "Generating chromosome sizes file..."
    fetchChromSizes mm10 > ${CHROM_SIZES}
fi

echo "=== Starting pipeline for ${SAMPLE} ==="

# Step 1: Initial quality control with FastQC
echo "Running initial FastQC..."
fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

# Step 2: Trim adapters and low quality bases with Trimmomatic
echo "Starting adapter trimming..."
# Download Illumina adapter sequences if not present
if [ ! -f "TruSeq3-PE.fa" ]; then
    wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi

# Run Trimmomatic in paired-end mode
trimmomatic PE -threads 32 \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

# Step 3: Align reads to mm10 genome with Bowtie2
echo "Starting alignment..."
# Align with sensitive parameters and filter for high quality paired reads
bowtie2 -p 32 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --maxins 1000 \
    -x ${WORKING_DIR}/mm10/mm10 \
    -1 ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
    -2 ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz 2> ${RESULTS_DIR}/bowtie2_alt/${SAMPLE}_align.log | \
    samtools view -q 30 -F 1804 -f 2 -b | \
    samtools sort -@ 32 -o ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam -

# Step 4: Add read groups to BAM header
echo "Adding read groups to header..."
# Create temporary header with read group info
samtools view -H ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam > temp_header.sam
if ! grep -q "^@RG" temp_header.sam; then
    echo -e "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA" >> temp_header.sam
fi
# Apply new header and add RG tags to reads
samtools reheader temp_header.sam ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam | \
    samtools view -h | \
    awk -v sample=$SAMPLE 'BEGIN {OFS="\t"} /^@/ {print} !/^@/ {print $0 "\tRG:Z:"sample}' | \
    samtools view -b > ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
rm temp_header.sam

# Index BAM file
samtools index ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam

# Remove intermediate BAM to save space
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.sorted.bam

# Step 5: Remove PCR duplicates with Picard
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam \
    OUTPUT=${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/${SAMPLE}.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

# Generate comprehensive QC metrics report
echo "Generating QC metrics..."
{
    echo "=== QC Metrics for ${SAMPLE} ==="
    echo "Initial read counts:"
    echo "Raw reads: $(zcat ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz | wc -l | awk '{print $1/4}')"
    echo "After trimming: $(zcat ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz | wc -l | awk '{print $1/4}')"
    
    echo -e "\nAlignment stats (pre-deduplication):"
    samtools flagstat ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
    
    echo -e "\nPost-deduplication stats:"
    samtools flagstat ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam
    
    echo -e "\nDuplication metrics:"
    cat ${RESULTS_DIR}/filtered/${SAMPLE}.metrics.txt
    
    echo -e "\nNumber of peaks:"
    wc -l ${RESULTS_DIR}/peaks/${SAMPLE}_peaks.narrowPeak
} > ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

# Remove intermediate BAM files
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
rm ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam.bai

# Step 6: Generate coverage tracks as bigWig files
echo "Generating bigWig files..."
# Generate RPKM normalized coverage
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${RESULTS_DIR}/bigwig/${SAMPLE}_RPKM.bw \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 32 \
    --extendReads \
    --ignoreDuplicates

# Generate CPM normalized coverage
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${RESULTS_DIR}/bigwig/${SAMPLE}_CPM.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 32 \
    --extendReads \
    --ignoreDuplicates

# Step 7: Call peaks with MACS2
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

echo "=== Pipeline completed successfully for ${SAMPLE}! ==="