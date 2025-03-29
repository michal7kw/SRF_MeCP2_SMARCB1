#!/bin/bash

# SLURM submission script for advanced peak visualization
# This script runs the R script for creating customized pie charts and additional visualizations

#SBATCH --job-name=3c_run_advanced_visualization
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error="logs/3c_run_advanced_visualization.err"
#SBATCH --output="logs/3c_run_advanced_visualization.out"

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
required_packages <- c("ggplot2", "dplyr", "gridExtra", "reshape2", 
                      "RColorBrewer", "ChIPseeker", "GenomicFeatures")
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

# Check if annotation files exist
echo "Checking if annotation files exist..."
if [ ! -f "results/annotations/BG1_peak_annotation.csv" ]; then
    echo "Annotation files not found. Please run 2_peak_annotation.R first."
    echo "Running 2_peak_annotation.R now..."
    Rscript 2_peak_annotation.R
    
    if [ $? -ne 0 ]; then
        echo "Error: Peak annotation failed. Cannot proceed with advanced visualization."
        exit 1
    fi
fi

# Run the advanced visualization script
echo "Starting advanced visualization..."
Rscript 3c_advanced_peak_visualization.R

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Advanced visualization completed successfully!"
    
    # Count the number of files generated
    num_pdf_files=$(find results/annotations/advanced_plots -name "*.pdf" | wc -l)
    num_png_files=$(find results/annotations/advanced_plots -name "*.png" | wc -l)
    
    echo "Generated $num_pdf_files PDF files and $num_png_files PNG files."
    echo "Results are available in the 'results/annotations/advanced_plots' directory."
else
    echo "Error: Advanced visualization failed. Check the error log for details."
    exit 1
fi 