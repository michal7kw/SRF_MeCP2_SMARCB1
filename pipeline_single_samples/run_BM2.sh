#!/bin/bash
#SBATCH --job-name=BM2
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM2.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM2.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define input/output directories
FASTQ_DIR="failed_data"
RESULTS_DIR="results"
mkdir -p ${RESULTS_DIR}/{trimmed,filtered,bowtie2_alt,peaks_alt}

# # 1. Initial FastQC (commented out as requested)
# echo "Running initial FastQC..."
# mkdir -p ${RESULTS_DIR}/fastqc
# fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
#     ${FASTQ_DIR}/BM2_R1_001.fastq.gz \
#     ${FASTQ_DIR}/BM2_R2_001.fastq.gz

# 1. Trim adapters and low quality bases with Cut&Tag optimized parameters
echo "Starting adapter trimming..."
trimmomatic PE -threads 32 \
    ${FASTQ_DIR}/BM2_R1_001.fastq.gz ${FASTQ_DIR}/BM2_R2_001.fastq.gz \
    ${RESULTS_DIR}/trimmed/BM2_R1_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/BM2_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/BM2_R2_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/BM2_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

# 2. Align with bowtie2 using Cut&Tag specific parameters
echo "Starting alignment..."
bowtie2 -p 32 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --maxins 1000 \
    -x /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/mm10/mm10 \
    -1 ${RESULTS_DIR}/trimmed/BM2_R1_trimmed.fastq.gz \
    -2 ${RESULTS_DIR}/trimmed/BM2_R2_trimmed.fastq.gz 2> ${RESULTS_DIR}/bowtie2_alt/BM2_align.log | \
    samtools view -h | \
    awk 'BEGIN {OFS="\t"} /^@/ {print} !/^@/ {print $0 "\tRG:Z:BM2"}' | \
    samtools view -q 30 -F 1804 -f 2 -b | \
    samtools sort -@ 32 -o ${RESULTS_DIR}/bowtie2_alt/BM2.sorted.bam -

# 3. Index the BAM file
echo "Indexing BAM file..."
samtools index ${RESULTS_DIR}/bowtie2_alt/BM2.sorted.bam

# 4. Remove PCR duplicates
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2_alt/BM2.sorted.bam \
    OUTPUT=${RESULTS_DIR}/filtered/BM2.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/BM2.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# 5. Index the deduplicated BAM
samtools index ${RESULTS_DIR}/filtered/BM2.dedup.bam

# 6. Call peaks using MACS2 with Cut&Tag specific parameters
echo "Calling peaks..."
macs2 callpeak \
    -t ${RESULTS_DIR}/filtered/BM2.dedup.bam \
    -f BAMPE \
    -g mm \
    -n BM2 \
    --outdir ${RESULTS_DIR}/peaks_alt \
    --nomodel \
    --extsize 200 \
    --keep-dup all \
    --qvalue 0.05 \
    --call-summits

# 7. Generate basic QC metrics
echo "Generating QC metrics..."
# Count initial reads
echo "Initial read counts:" > ${RESULTS_DIR}/BM2_qc_metrics.txt
zcat ${FASTQ_DIR}/BM2_R1_001.fastq.gz | echo "Raw reads: $((`wc -l`/4))" >> ${RESULTS_DIR}/BM2_qc_metrics.txt

# Count reads after trimming
zcat ${RESULTS_DIR}/trimmed/BM2_R1_trimmed.fastq.gz | echo "After trimming: $((`wc -l`/4))" >> ${RESULTS_DIR}/BM2_qc_metrics.txt

# Get alignment stats
samtools flagstat ${RESULTS_DIR}/bowtie2_alt/BM2.sorted.bam >> ${RESULTS_DIR}/BM2_qc_metrics.txt

# Get final dedup stats
samtools flagstat ${RESULTS_DIR}/filtered/BM2.dedup.bam >> ${RESULTS_DIR}/BM2_qc_metrics.txt

# Count final peaks
wc -l ${RESULTS_DIR}/peaks_alt/BM2_peaks.narrowPeak | cut -d' ' -f1 >> ${RESULTS_DIR}/BM2_qc_metrics.txt

echo "Pipeline completed successfully!"