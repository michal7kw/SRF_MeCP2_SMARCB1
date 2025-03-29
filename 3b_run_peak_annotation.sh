#!/bin/bash

# SLURM submission script for peak annotation
# This script runs the R script for annotating ChIP-seq peaks and generating pie charts

#SBATCH --job-name=3b_run_peak_annotation
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/3b_run_peak_annotation.err"
#SBATCH --output="logs/3b_run_peak_annotation.out"

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
required_packages <- c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                      "org.Mm.eg.db", "GenomicFeatures", "rtracklayer", 
                      "ggplot2", "dplyr", "gridExtra")
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

# Run the peak annotation script
echo "Starting peak annotation..."
Rscript 3b_peak_annotation.R

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Peak annotation completed successfully!"
    
    # Count the number of files generated
    num_csv_files=$(find results/annotations -name "*.csv" | wc -l)
    num_pdf_files=$(find results/annotations/plots -name "*.pdf" | wc -l)
    
    echo "Generated $num_csv_files CSV files and $num_pdf_files PDF files."
    echo "Results are available in the 'results/annotations' directory."
else
    echo "Error: Peak annotation failed. Check the error log for details."
    exit 1
fi 