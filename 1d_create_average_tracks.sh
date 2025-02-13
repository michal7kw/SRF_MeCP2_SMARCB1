#!/bin/bash

# This script creates average signal tracks from bigWig files for BG and BM sample groups

#SBATCH --job-name=1d_create_average_tracks
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/1d_create_average_tracks.err"
#SBATCH --output="logs/1d_create_average_tracks.out"

# Exit on error
set -e

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create output directory for metaprofiles
mkdir -p results/metaprofiles

# Define sample groups for analysis
BG_SAMPLES="BG1 BG2 BG3"  # Control samples
BM_SAMPLES="BM3"          # Treatment sample

BE_DIR="bigwig"

# Create lists of bigwig files for each sample group
BG_BIGWIGS=""
for sample in $BG_SAMPLES; do
    BG_BIGWIGS="$BG_BIGWIGS results/${BE_DIR}/${sample}_CPM.bw"
done

BM_BIGWIGS=""
for sample in $BM_SAMPLES; do
    BM_BIGWIGS="$BM_BIGWIGS results/${BE_DIR}/${sample}_CPM.bw"
done

# Function to check if bigwig files exist and are valid
check_bigwig() {
    local bw_file=$1
    if [ ! -f "$bw_file" ]; then
        echo "Error: BigWig file $bw_file not found"
        return 1
    fi
    
    # Try to read the bigwig file header using bigWigInfo
    if ! bigWigInfo $bw_file > /dev/null 2>&1; then
        echo "Error: BigWig file $bw_file is not valid"
        return 1
    fi
    return 0
}

echo "Creating average bigWig files..."

# Check input bigwig files
for bw in $BG_BIGWIGS $BM_BIGWIGS; do
    if ! check_bigwig "$bw"; then
        echo "Error: Invalid input bigwig file: $bw"
        exit 1
    fi
done

# Create average bigwig files with error checking
echo "Averaging BG samples..."
if ! bigwigAverage -b $BG_BIGWIGS -o results/metaprofiles/BG_average.bw; then
    echo "Error: Failed to create BG average bigwig"
    exit 1
fi

echo "Averaging BM samples..."
if ! bigwigAverage -b $BM_BIGWIGS -o results/metaprofiles/BM_average.bw; then
    echo "Error: Failed to create BM average bigwig"
    exit 1
fi

# Verify the created files
if ! check_bigwig "results/metaprofiles/BG_average.bw" || ! check_bigwig "results/metaprofiles/BM_average.bw"; then
    echo "Error: Failed to create valid average bigwig files"
    exit 1
fi

# Create flag file to indicate successful completion
echo "BigWig averaging completed successfully"
touch results/metaprofiles/.bigwig_averaging_complete

echo "Average signal tracks have been created successfully!" 