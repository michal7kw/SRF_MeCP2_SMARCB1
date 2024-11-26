#!/bin/bash
#SBATCH --job-name=Python
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/python.out"

set -e  # Exit on error
set -x  # Print commands being executed

# Create necessary directories
mkdir -p results/peak_analysis_alt

cd /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta

# # First merge peaks for each group
# python scripts/merge_peak_groups.py \
#     --input-peaks results/peaks/BG*_peaks.narrowPeak \
#     --output results/peak_analysis_alt/BG_merged.bed

# python scripts/merge_peak_groups.py \
#     --input-peaks results/peaks_alt/BM*_peaks.narrowPeak \
#     --output results/peak_analysis_alt/BM_merged.bed

# # Find common peaks in CpG islands
# python scripts/find_common_cpg_peaks.py \
#     --bg-peaks results/peak_analysis_alt/BG_merged.bed \
#     --bm-peaks results/peak_analysis_alt/BM_merged.bed \
#     --cpg-islands cpg_islands.bed \
#     --output results/peak_analysis_alt/common_cpg_peaks.bed

# Count reads for each sample - process one file at a time
for bam in results/bowtie2_alt/{BG,BM}*.sorted.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename "$bam" .sorted.bam)
        echo "Processing $sample..."
        python scripts/count_reads_in_peaks.py \
            --peaks results/peak_analysis_alt/common_cpg_peaks.bed \
            --bam "$bam" \
            --output "results/peak_analysis_alt/${sample}_peak_counts.txt"
    fi
done

# # Wait for all count jobs to finish
# wait

# # Check if any count files were created
# count_files=(results/peak_analysis_alt/*_peak_counts.txt)
# if [ ${#count_files[@]} -eq 0 ]; then
#     echo "No count files were created. Exiting."
#     exit 1
# fi

# # Compare peak sizes (explicitly list the files to avoid shell expansion issues)
# python scripts/compare_peak_sizes.py \
#     --peak-counts results/peak_analysis_alt/BG*_peak_counts.txt results/peak_analysis_alt/BM*_peak_counts.txt \
#     --output results/peak_analysis_alt/peak_size_comparison.txt
