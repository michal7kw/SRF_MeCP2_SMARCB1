#!/bin/bash
#SBATCH --job-name=BM1_continue
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM1_continue.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM1_continue.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

RESULTS_DIR="results"

# Add read groups to existing BAM file
echo "Adding read groups..."
samtools addreplacerg \
    -r 'ID:BM1' -r 'SM:BM1' -r 'LB:lib1' -r 'PL:ILLUMINA' -r 'PU:unit1' \
    -o ${RESULTS_DIR}/bowtie2_alt/BM1.sorted.withRG.bam \
    ${RESULTS_DIR}/bowtie2_alt/BM1.sorted.bam

# 4. Remove PCR duplicates
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${RESULTS_DIR}/bowtie2_alt/BM1.sorted.withRG.bam \
    OUTPUT=${RESULTS_DIR}/filtered/BM1.dedup.bam \
    METRICS_FILE=${RESULTS_DIR}/filtered/BM1.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

# 6. Call peaks using MACS2 with Cut&Tag specific parameters
echo "Calling peaks..."
macs2 callpeak \
    -t ${RESULTS_DIR}/filtered/BM1.dedup.bam \
    -f BAMPE \
    -g mm \
    -n BM1 \
    --outdir ${RESULTS_DIR}/peaks_alt \
    --nomodel \
    --extsize 200 \
    --keep-dup all \
    --qvalue 0.05 \
    --call-summits

# 7. Generate basic QC metrics
echo "Generating QC metrics..."
echo "Initial read counts:" > ${RESULTS_DIR}/BM1_qc_metrics.txt

# Skip the initial fastq counts since we don't need them
echo "=== Alignment Stats ===" >> ${RESULTS_DIR}/BM1_qc_metrics.txt
samtools flagstat ${RESULTS_DIR}/bowtie2_alt/BM1.sorted.withRG.bam >> ${RESULTS_DIR}/BM1_qc_metrics.txt

echo "=== After Deduplication Stats ===" >> ${RESULTS_DIR}/BM1_qc_metrics.txt
samtools flagstat ${RESULTS_DIR}/filtered/BM1.dedup.bam >> ${RESULTS_DIR}/BM1_qc_metrics.txt

echo "=== Peak Counts ===" >> ${RESULTS_DIR}/BM1_qc_metrics.txt
wc -l ${RESULTS_DIR}/peaks_alt/BM1_peaks.narrowPeak | cut -d' ' -f1 >> ${RESULTS_DIR}/BM1_qc_metrics.txt

echo "Pipeline completed successfully!" 