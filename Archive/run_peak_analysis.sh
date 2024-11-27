#!/bin/bash

# Set paths
PEAK_COMPARISON="results/peak_analysis_alt/peak_size_comparison.txt"
BG_COUNTS="results/peak_analysis_alt/${BG_SAMPLE}_peak_counts.txt"
BM_COUNTS="results/peak_analysis_alt/${BM_SAMPLE}_peak_counts.txt"
OUTPUT_DIR="results/peak_analysis_alt/reanalysis"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run diagnostic analysis
python scripts/diagnose_peaks.py \
    --input "$PEAK_COMPARISON" \
    --output-dir "$OUTPUT_DIR/diagnostics"

# Run reanalysis
python scripts/reanalyze_peaks.py \
    --bg-counts "$BG_COUNTS" \
    --bm-counts "$BM_COUNTS" \
    --output-dir "$OUTPUT_DIR/results" 