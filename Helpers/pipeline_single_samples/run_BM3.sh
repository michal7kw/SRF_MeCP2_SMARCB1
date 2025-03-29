#!/bin/bash
#SBATCH --job-name=BM3
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM3.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM3.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define input/output directories
FASTQ_DIR="90-1102945428/00_fastq"
RESULTS_DIR="results"

# Create directory structure
# We create separate directories for each step of the pipeline:
# - fastqc: Quality control reports
# - trimmed: Adapter-trimmed reads
# - filtered: PCR duplicate removed files
# - bowtie2_alt: Alignment files
# - peaks_alt: MACS2 peak calls
# - seacr_peaks: SEACR peak calls (alternative peak caller for CUT&Tag)
mkdir -p ${RESULTS_DIR}/{fastqc,trimmed,filtered,bowtie2_alt,peaks_alt,seacr_peaks}

# Directory check and workspace information
if [ ! -d "${RESULTS_DIR}/filtered" ]; then
    echo "Error: Failed to create filtered directory"
    exit 1
fi

echo "Working directory: $(pwd)"
echo "Contents of results directory:"
ls -l ${RESULTS_DIR}

# 1. Initial FastQC
# Important for CUT&Tag QC as it helps identify:
# - Quality of base calls
# - Presence of adapters
# - Sequence complexity
# - GC content distribution (important for CUT&Tag)
echo "Running initial FastQC..."
fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
    ${FASTQ_DIR}/BM3_R1_001.fastq.gz \
    ${FASTQ_DIR}/BM3_R2_001.fastq.gz

# 2. Trim adapters
# CUT&Tag typically produces cleaner data than ChIP-seq, so we use less stringent parameters:
# - LEADING/TRAILING:3 - Remove bases with quality below 3
# - SLIDINGWINDOW:4:15 - Moderate quality filtering
# - MINLEN:20 - Keep reads ≥20bp (shorter than ChIP-seq as CUT&Tag produces cleaner data)
# - ILLUMINACLIP - Remove Illumina adapters while keeping both reads when possible
echo "Starting adapter trimming..."
trimmomatic PE -threads 32 \
    ${FASTQ_DIR}/BM3_R1_001.fastq.gz ${FASTQ_DIR}/BM3_R2_001.fastq.gz \
    ${RESULTS_DIR}/trimmed/BM3_R1_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/BM3_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/BM3_R2_trimmed.fastq.gz ${RESULTS_DIR}/trimmed/BM3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

# 3. Post-trim FastQC
# Verify the effectiveness of trimming and ensure read quality
fastqc -t 32 -o ${RESULTS_DIR}/fastqc \
    ${RESULTS_DIR}/trimmed/BM3_R1_trimmed.fastq.gz \
    ${RESULTS_DIR}/trimmed/BM3_R2_trimmed.fastq.gz

# 4. Align with bowtie2
# CUT&Tag specific parameters:
# --end-to-end: Use full read length (recommended for CUT&Tag)
# --very-sensitive: Thorough alignment for high accuracy
# --no-mixed/--no-discordant: Only keep properly paired reads
# -I 10 -X 700: Accept insert sizes between 10-700bp (typical for CUT&Tag)
# -q 30: High-quality alignments only
# -F 1804: Remove unmapped, mate unmapped, secondary alignments, and duplicates
# -f 2: Keep properly paired reads
echo "Starting alignment..."
bowtie2 -p 32 \
    --end-to-end \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 -X 700 \
    -x /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/mm10/mm10 \
    -1 ${RESULTS_DIR}/trimmed/BM3_R1_trimmed.fastq.gz \
    -2 ${RESULTS_DIR}/trimmed/BM3_R2_trimmed.fastq.gz 2> ${RESULTS_DIR}/bowtie2_alt/BM3_align.log | \
    samtools view -h | \
    awk 'BEGIN {OFS="\t"} /^@/ {print} !/^@/ {print $0 "\tRG:Z:BM3"}' | \
    samtools view -q 30 -F 1804 -f 2 -b | \
    samtools sort -@ 32 -o ${RESULTS_DIR}/bowtie2_alt/BM3.sorted.bam -

# 5-8. Fragment analysis
# Critical for CUT&Tag as fragment size distribution helps verify successful tagmentation
# Convert BAM to BEDPE for easier fragment analysis
echo "Performing fragment analysis..."
samtools index ${RESULTS_DIR}/bowtie2_alt/BM3.sorted.bam

# Convert to BEDPE format
bedtools bamtobed -i ${RESULTS_DIR}/bowtie2_alt/BM3.sorted.bam -bedpe > ${RESULTS_DIR}/bowtie2_alt/BM3.bedpe

# Filter fragments: Keep only proper pairs (<1000bp)
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${RESULTS_DIR}/bowtie2_alt/BM3.bedpe > ${RESULTS_DIR}/bowtie2_alt/BM3.clean.bedpe

# Calculate fragment lengths for QC
awk '{print $6-$2}' ${RESULTS_DIR}/bowtie2_alt/BM3.clean.bedpe > ${RESULTS_DIR}/bowtie2_alt/BM3.fragmentLen.txt

# 9. PCR Duplicate Removal
# Optional for CUT&Tag as duplicates might be genuine signal
# However, removing them can increase confidence in peak calls
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2_alt/BM3.sorted.bam \
    OUTPUT=${RESULTS_DIR}/filtered/BM3.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/BM3.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# Error checking and indexing
if [ ! -f "${RESULTS_DIR}/filtered/BM3.dedup.bam" ]; then
    echo "Error: Deduplication failed"
    exit 1
fi

samtools index ${RESULTS_DIR}/filtered/BM3.dedup.bam

# 11. Generate normalized coverage
# Create bedgraph for visualization and SEACR input
# -bg: bedgraph format
# -pc: calculate coverage of paired-end fragments
bedtools genomecov -bg -pc -ibam ${RESULTS_DIR}/filtered/BM3.dedup.bam > ${RESULTS_DIR}/filtered/BM3.bedgraph

# 12. Peak Calling with MACS2
# CUT&Tag specific parameters:
# --nomodel: Don't build shifting model (different from ChIP-seq)
# --shift -100: Adjust for Tn5 binding characteristics
# --extsize 200: Typical fragment size for CUT&Tag
# --keep-dup all: Consider all reads (important for CUT&Tag)
# --qvalue 0.01: Stringent peak calling
echo "Calling peaks with MACS2..."
macs2 callpeak \
    -t ${RESULTS_DIR}/filtered/BM3.dedup.bam \
    -f BAMPE \
    -g mm \
    -n BM3 \
    --outdir ${RESULTS_DIR}/peaks_alt \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    --keep-dup all \
    --qvalue 0.01 \
    --call-summits

# 13. Quality Control Metrics
# Comprehensive QC including:
# - Read counts at each step
# - Alignment statistics
# - Fragment length distribution
# - Peak statistics
echo "Generating QC metrics..."
echo "Initial read counts:" > ${RESULTS_DIR}/BM3_qc_metrics.txt
zcat ${FASTQ_DIR}/BM3_R1_001.fastq.gz | echo "Raw reads: $((`wc -l`/4))" >> ${RESULTS_DIR}/BM3_qc_metrics.txt
zcat ${RESULTS_DIR}/trimmed/BM3_R1_trimmed.fastq.gz | echo "After trimming: $((`wc -l`/4))" >> ${RESULTS_DIR}/BM3_qc_metrics.txt

# Detailed alignment statistics
samtools flagstat ${RESULTS_DIR}/bowtie2_alt/BM3.sorted.bam >> ${RESULTS_DIR}/BM3_qc_metrics.txt
samtools flagstat ${RESULTS_DIR}/filtered/BM3.dedup.bam >> ${RESULTS_DIR}/BM3_qc_metrics.txt

# Fragment length distribution statistics
# Important for verifying successful CUT&Tag
awk '{sum+=$1; sumsq+=$1*$1} END {print "Mean fragment length: "sum/NR; print "StdDev fragment length: "sqrt(sumsq/NR - (sum/NR)**2)}' \
    ${RESULTS_DIR}/bowtie2_alt/BM3.fragmentLen.txt >> ${RESULTS_DIR}/BM3_qc_metrics.txt

# Final peak count
wc -l ${RESULTS_DIR}/peaks_alt/BM3_peaks.narrowPeak | cut -d' ' -f1 >> ${RESULTS_DIR}/BM3_qc_metrics.txt

echo "Pipeline completed successfully!"