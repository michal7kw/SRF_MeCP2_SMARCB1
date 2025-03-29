#!/bin/bash
#SBATCH --job-name=CutTag_array
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=1-4
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/%A_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/%A_%a.out"

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

# Define input/output directories
FASTQ_DIR="90-1102945428/00_fastq"
RESULTS_DIR="results"
mkdir -p ${RESULTS_DIR}/{trimmed,filtered,bowtie2_alt,peaks_alt}

# 1. Initial FastQC
echo "Running initial FastQC..."
mkdir -p ${RESULTS_DIR}/fastqc
fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

# 2. Trim adapters and low quality bases with Cut&Tag optimized parameters
echo "Starting adapter trimming..."
trimmomatic PE -threads 32 \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

# 3. Align with bowtie2 using Cut&Tag specific parameters
echo "Starting alignment..."
bowtie2 -p 32 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --maxins 1000 \
    -x /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/mm10/mm10 \
    -1 ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
    -2 ${RESULTS_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz 2> ${RESULTS_DIR}/bowtie2_alt/${SAMPLE}_align.log | \
    samtools view -h | \
    awk -v sample=$SAMPLE 'BEGIN {OFS="\t"} /^@/ {print} !/^@/ {print $0 "\tRG:Z:"sample}' | \
    samtools view -q 30 -F 1804 -f 2 -b | \
    samtools sort -@ 32 -o ${RESULTS_DIR}/bowtie2_alt/${SAMPLE}.sorted.bam -

# 4. Index the BAM file
echo "Indexing BAM file..."
samtools index ${RESULTS_DIR}/bowtie2_alt/${SAMPLE}.sorted.bam

# 5. Remove PCR duplicates
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2_alt/${SAMPLE}.sorted.bam \
    OUTPUT=${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/${SAMPLE}.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# 6. Index the deduplicated BAM
samtools index ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam

# 7. Call peaks using MACS2 with Cut&Tag specific parameters
echo "Calling peaks..."
macs2 callpeak \
    -t ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    -f BAMPE \
    -g mm \
    -n ${SAMPLE} \
    --outdir ${RESULTS_DIR}/peaks_alt \
    --nomodel \
    --extsize 200 \
    --keep-dup all \
    --qvalue 0.05 \
    --call-summits

# 8. Generate basic QC metrics
echo "Generating QC metrics..."
# Count initial reads
echo "Initial read counts:" > ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt
zcat ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz | echo "Raw reads: $((`wc -l`/4))" >> ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

# Count reads after trimming
zcat ${RESULTS_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz | echo "After trimming: $((`wc -l`/4))" >> ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

# Get alignment stats
samtools flagstat ${RESULTS_DIR}/bowtie2_alt/${SAMPLE}.sorted.bam >> ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

# Get final dedup stats
samtools flagstat ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam >> ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

# Count final peaks
wc -l ${RESULTS_DIR}/peaks_alt/${SAMPLE}_peaks.narrowPeak | cut -d' ' -f1 >> ${RESULTS_DIR}/${SAMPLE}_qc_metrics.txt

echo "Pipeline completed successfully for ${SAMPLE}!"