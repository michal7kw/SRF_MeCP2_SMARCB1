#!/bin/bash
#SBATCH --job-name=Python_Serial
#SBATCH --account=kubacki.michal
#SBATCH --mem=256GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python_serial.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python_serial.out"

set -e  # Exit on error
set -x  # Print commands being executed

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

BG_SAMPLE="BG3"
BM_SAMPLE="BM3"

cp results/peaks/${BG_SAMPLE}_peaks.narrowPeak results/peak_analysis/${BG_SAMPLE}_merged.bed
cp results/peaks/${BM_SAMPLE}_peaks.narrowPeak results/peak_analysis/${BM_SAMPLE}_merged.bed

# Sort input files before processing
for sample in ${BG_SAMPLE} ${BM_SAMPLE}; do
    # Sort BED files
    sort -k1,1 -k2,2n results/peak_analysis/${sample}_merged.bed > results/peak_analysis/${sample}_merged.sorted.bed
    mv results/peak_analysis/${sample}_merged.sorted.bed results/peak_analysis/${sample}_merged.bed

    # Ensure BAM files are sorted
    if ! samtools view -H results/bowtie2/${sample}.sorted.bam | grep -q '@HD.*SO:coordinate'; then
        samtools sort results/bowtie2/${sample}.sorted.bam -o results/bowtie2/${sample}.sorted.tmp.bam
        mv results/bowtie2/${sample}.sorted.tmp.bam results/bowtie2/${sample}.sorted.bam
        samtools index results/bowtie2/${sample}.sorted.bam
    fi
done

# Find common peaks in CpG islands
python scripts/find_common_cpg_peaks_serial.py \
    --bg-peaks results/peak_analysis/${BG_SAMPLE}_merged.bed \
    --bm-peaks results/peak_analysis/${BM_SAMPLE}_merged.bed \
    --cpg-islands cpg_islands.bed \
    --output results/peak_analysis/common_cpg_peaks.bed

# Count reads for each sample
for bam in results/bowtie2/${BG_SAMPLE}.sorted.bam results/bowtie2/${BM_SAMPLE}.sorted.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename "$bam" .sorted.bam)
        echo "Processing $sample..."
        python scripts/count_reads_in_peaks_serial.py \
            --peaks results/peak_analysis/common_cpg_peaks.bed \
            --bam "$bam" \
            --output "results/peak_analysis/${sample}_peak_counts.txt"
    fi
done

# Wait for all count jobs to finish
wait

# Check if any count files were created
count_files=(results/peak_analysis/*_peak_counts.txt)
if [ ${#count_files[@]} -eq 0 ]; then
    echo "No count files were created. Exiting."
    exit 1
fi

# Compare peak sizes
python scripts/compare_peak_sizes_serial.py \
    --peak-counts results/peak_analysis/${BG_SAMPLE}_peak_counts.txt results/peak_analysis/${BM_SAMPLE}_peak_counts.txt \
    --output results/peak_analysis/peak_size_comparison.txt
