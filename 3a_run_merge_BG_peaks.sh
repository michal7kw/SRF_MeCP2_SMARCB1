#!/bin/bash

# SLURM submission script for merging BG peaks
# This script runs the R script for identifying consistent peaks across BG replicates

#SBATCH --job-name=3a_run_merge_BG_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error="logs/3a_run_merge_BG_peaks.err"
#SBATCH --output="logs/3a_run_merge_BG_peaks.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required software environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Check if required R packages are installed
echo "Checking required R packages..."
Rscript -e '
required_packages <- c("GenomicRanges", "rtracklayer", "dplyr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse=", "), "\n")
  cat("Installing missing packages...\n")
  if(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")
  BiocManager::install(missing_packages)
} else {
  cat("All required packages are installed.\n")
}'

# Check if input peak files exist
echo "Checking if peak files exist..."
for file in results/peaks/BG{1,2,3}_peaks.narrowPeak; do
    if [ ! -f "$file" ]; then
        echo "Error: $file not found!"
        exit 1
    fi
done

# Run the peak merging script
echo "Starting peak merging..."
Rscript 3a_merge_BG_peaks.R

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Peak merging completed successfully!"
    
    if [ -f "results/peaks/BG_consistent_peaks.narrowPeak" ]; then
        num_peaks=$(wc -l < results/peaks/BG_consistent_peaks.narrowPeak)
        echo "Generated BG_consistent_peaks.narrowPeak with $num_peaks peaks"
    else
        echo "Error: Output file was not created"
        exit 1
    fi
else
    echo "Error: Peak merging failed. Check the error log for details."
    exit 1
fi 