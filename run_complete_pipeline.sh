#!/bin/bash
#SBATCH --job-name=SMARCB1_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=1-4
#SBATCH --error="logs/%A_%a.err"
#SBATCH --output="logs/%A_%a.out"

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
FASTQ_DIR="90-1102945428/00_fastq"
RESULTS_DIR="results"
mkdir -p ${RESULTS_DIR}/{trimmed,filtered,bowtie2,peaks,bigwig,fastqc}

# Generate chromosome sizes file if not already present
CHROM_SIZES="${WORKING_DIR}/mm10/mm10.chrom.sizes"
if [ ! -f "$CHROM_SIZES" ]; then
    echo "Generating chromosome sizes file..."
    fetchChromSizes mm10 > ${CHROM_SIZES}
fi

echo "=== Starting pipeline for ${SAMPLE} ==="

# # 1. Initial FastQC
# echo "Running initial FastQC..."
# fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
#     ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
#     ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

# 2. Trim adapters and low quality bases
echo "Starting adapter trimming..."
# Download adapter file if it doesn't exist
if [ ! -f "TruSeq3-PE.fa" ]; then
    wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi

trimmomatic PE -threads 32 \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

# 3. Align with bowtie2
echo "Starting alignment..."
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

# 4. Add read groups to header
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

# 5. Remove duplicates
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

# 6. Index deduplicated BAM (not needed since CREATE_INDEX=true above)
# samtools index ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam

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
    echo "Initial read counts:"
    echo "Raw reads: $(zcat ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz | wc -l | awk '{print $1/4}')"
    echo "After trimming: $(zcat ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz | wc -l | awk '{print $1/4}')"
    echo -e "\nAlignment stats:"
    samtools flagstat ${RESULTS_DIR}/bowtie2/${SAMPLE}.with_rg.bam
    echo -e "\nPost-deduplication stats:"
    samtools flagstat ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam
    echo -e "\nNumber of peaks:"
    wc -l ${RESULTS_DIR}/peaks/${SAMPLE}_peaks.narrowPeak
} > ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

echo "=== Pipeline completed successfully for ${SAMPLE}! ==="