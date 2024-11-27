#!/bin/bash
#SBATCH --job-name=Python
#SBATCH --account=kubacki.michal
#SBATCH --mem=256GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.out"

set -e  # Exit on error
set -x  # Print commands being executed

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set number of cores for parallel processing
export NCORES=32

# Enable parallel processing for samtools
export SAMTOOLS_THREADS=$((NCORES/2))

BG_SAMPLE="BG3"
BM_SAMPLE="BM3"

# Validate input files and create directories
for file in "results/peaks_alt/${BG_SAMPLE}_peaks.narrowPeak" "results/peaks_alt/${BM_SAMPLE}_peaks.narrowPeak" "cpg_islands.bed"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required input file $file not found"
        exit 1
    fi
done

# Create directories with proper error checking
for dir in "results/peak_analysis_alt" "results/peaks_alt" "results/bowtie2_alt"; do
    mkdir -p "$dir" || { echo "Error creating directory: $dir"; exit 1; }
done

# Process each sample's peaks
for sample in ${BG_SAMPLE} ${BM_SAMPLE}; do
    echo "Processing ${sample}..."
    
    # Count reads in peaks
    if [ -f "results/bowtie2_alt/${sample}.sorted.bam" ]; then
        python scripts/count_reads_in_peaks.py \
            --peaks results/peak_analysis_alt/common_cpg_peaks.bed \
            --bam results/bowtie2_alt/${sample}.sorted.bam \
            --output results/peak_analysis_alt/${sample}_peak_counts.txt \
            --threads $NCORES
            
        # Validate output
        if [ ! -s "results/peak_analysis_alt/${sample}_peak_counts.txt" ]; then
            echo "Error: Empty or missing count file for ${sample}"
            exit 1
        fi
    else
        echo "Error: BAM file not found for ${sample}"
        exit 1
    fi
done

# Compare peak sizes with proper error handling
if [ -f "results/peak_analysis_alt/${BG_SAMPLE}_peak_counts.txt" ] && \
   [ -f "results/peak_analysis_alt/${BM_SAMPLE}_peak_counts.txt" ]; then
    python scripts/compare_peak_sizes.py \
        --peak-counts results/peak_analysis_alt/${BG_SAMPLE}_peak_counts.txt \
                      results/peak_analysis_alt/${BM_SAMPLE}_peak_counts.txt \
        --output results/peak_analysis_alt/peak_size_comparison.txt \
        --threads $NCORES
else
    echo "Error: Missing count files for comparison"
    exit 1
fi