#!/bin/bash

#SBATCH --job-name=1c_fix_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/1c_fix_bigwig.err"
#SBATCH --output="logs/1c_fix_bigwig.out"

# Description:
# This script processes BAM files to generate normalized bigWig coverage files.
# It checks BAM file validity, creates indices, and generates RPKM-normalized 
# bigWig files for visualization.
#
# Inputs:
# - BAM files: results/filtered/<sample>.dedup.bam
#   Deduplicated BAM files containing aligned reads
#
# Outputs:
# - BAM indices: results/filtered/<sample>.dedup.bam.bai
#   Index files for the input BAM files
# - BigWig files: results/bigwig/<sample>_RPKM.bw
#   RPKM-normalized coverage files in bigWig format
#
# Dependencies:
# - samtools: For BAM file operations
# - deepTools: For bamCoverage
# - UCSC tools: For bigWigInfo

# Exit on error
set -e

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create output directory for bigwig files if it doesn't exist
mkdir -p results/bigwig

# Define samples
SAMPLES="BG1 BG2 BG3 BM3"

# Function to check BAM file validity
check_bam() {
    local bam_file=$1
    echo "Checking BAM file: $bam_file"
    
    # Check if file exists and is not empty
    if [ ! -s "$bam_file" ]; then
        echo "Error: BAM file $bam_file does not exist or is empty"
        return 1
    fi
    
    # Check if file is a valid BAM
    if ! samtools quickcheck "$bam_file"; then
        echo "Error: $bam_file is not a valid BAM file"
        return 1
    fi
    
    # Get basic stats
    echo "BAM file statistics:"
    samtools flagstat "$bam_file"
    
    return 0
}

# First, check and index all BAM files
echo "Checking and indexing BAM files..."
for sample in $SAMPLES; do
    bam_file="results/filtered/${sample}.dedup.bam"
    if [ -f "$bam_file" ]; then
        echo "Processing ${sample}.dedup.bam..."
        
        # Check BAM validity
        if ! check_bam "$bam_file"; then
            echo "Error: Invalid BAM file for ${sample}"
            exit 1
        fi
        
        # Index BAM file
        echo "Indexing ${sample}.dedup.bam..."
        samtools index "$bam_file"
        
        # Verify index was created
        if [ ! -f "${bam_file}.bai" ]; then
            echo "Error: Failed to create index for ${sample}"
            exit 1
        fi
    else
        echo "Error: BAM file for ${sample} not found!"
        exit 1
    fi
done

# Then regenerate bigWig files
echo "Generating bigWig files..."
for sample in $SAMPLES; do
    bam_file="results/filtered/${sample}.dedup.bam"
    if [ -f "$bam_file" ] && [ -f "${bam_file}.bai" ]; then
        echo "Processing ${sample}..."
        
        # Remove old bigWig file if it exists
        rm -f "results/bigwig/${sample}_RPKM.bw"
        
        # Generate new RPKM normalized coverage with verbose output
        echo "Running bamCoverage for ${sample}..."
        bamCoverage \
            --bam "$bam_file" \
            --outFileName "results/bigwig/${sample}_RPKM.bw" \
            --outFileFormat bigwig \
            --normalizeUsing RPKM \
            --binSize 10 \
            --numberOfProcessors 16 \
            --extendReads \
            --ignoreDuplicates \
            --verbose
        
        # Check if bamCoverage succeeded
        if [ $? -ne 0 ]; then
            echo "Error: bamCoverage failed for ${sample}"
            exit 1
        fi
        
        # Verify the generated bigWig file
        echo "Verifying bigWig file for ${sample}..."
        if ! bigWigInfo "results/bigwig/${sample}_RPKM.bw" > /dev/null 2>&1; then
            echo "Error: Failed to create valid bigWig file for ${sample}"
            exit 1
        fi
        
        echo "Successfully created bigWig file for ${sample}"
    else
        echo "Error: Required files for ${sample} not found!"
        exit 1
    fi
done

echo "BigWig file generation completed successfully!"

# Print final verification
echo -e "\nFinal verification of generated files:"
for sample in $SAMPLES; do
    echo "Checking ${sample}_RPKM.bw..."
    bigWigInfo "results/bigwig/${sample}_RPKM.bw" || echo "Warning: Issue with ${sample}_RPKM.bw"
done 