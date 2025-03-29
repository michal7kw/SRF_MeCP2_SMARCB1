#!/bin/bash

# Master script to run the entire peak annotation pipeline
# This script submits both the annotation and visualization jobs

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Create logs directory if it doesn't exist
mkdir -p logs

# Make all scripts executable
chmod +x 2_run_peak_annotation.sh
chmod +x 3_run_advanced_visualization.sh

echo "=== Starting Peak Annotation Pipeline ==="
echo "This pipeline will:"
echo "1. Annotate ChIP-seq peaks using the Gencode annotation"
echo "2. Generate pie charts and other visualizations"
echo "3. Create advanced customized visualizations"

# Submit the peak annotation job
echo "Submitting peak annotation job..."
annotation_job=$(sbatch 2_run_peak_annotation.sh | awk '{print $4}')
echo "Peak annotation job submitted with ID: $annotation_job"

# Submit the advanced visualization job with dependency on the annotation job
echo "Submitting advanced visualization job (will start after annotation completes)..."
viz_job=$(sbatch --dependency=afterok:$annotation_job 3_run_advanced_visualization.sh | awk '{print $4}')
echo "Advanced visualization job submitted with ID: $viz_job"

echo "=== Pipeline Submission Complete ==="
echo "To check job status, use: squeue -u $(whoami)"
echo "Results will be available in:"
echo "- results/annotations/ (annotation files)"
echo "- results/annotations/plots/ (standard plots)"
echo "- results/annotations/advanced_plots/ (customized visualizations)"

# Create a README file with information about the results
cat > results/annotations/README.txt << EOF
Peak Annotation Results
======================

This directory contains the results of the ChIP-seq peak annotation pipeline.

Directory Structure:
-------------------
- annotations/
  - *.csv: Peak annotation files with genomic feature information
  - plots/
    - *_genomic_feature_pie.pdf: Pie charts showing peak distribution
    - *_genomic_feature_bar.pdf: Bar plots showing peak distribution
    - *_TSS_distribution.pdf: Distribution of peaks relative to TSS
    - combined_*.pdf: Combined plots for all samples
  - advanced_plots/
    - *_custom_pie.pdf/png: Customized pie charts with percentages
    - all_samples_pie.pdf/png: Faceted pie charts for all samples
    - all_samples_stacked_bar.pdf/png: Stacked bar chart for comparison
    - peak_distribution_heatmap_data.csv: Data for creating heatmaps

Files:
------
- peak_annotation_summary.csv: Summary table with peak counts per feature for all samples

Analysis was performed using ChIPseeker with the Gencode vM10 annotation.
EOF

echo "A README file has been created in results/annotations/ with information about the results." 