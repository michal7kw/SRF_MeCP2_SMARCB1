#!/bin/bash

# This script generates metaprofiles for ChIP-seq data around transcription start sites (TSS)
# It performs the following steps:
# 1. Averages bigWig files for sample groups (BG and BM)
# 2. Processes gene lists to create TSS bed files
# 3. Generates matrices and plots showing ChIP-seq signal around TSS

#SBATCH --job-name=2_run_generate_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/2_run_generate_metaprofiles.err"
#SBATCH --output="logs/2_run_generate_metaprofiles.out"

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

# Create lists of bigwig files for each sample group
# These will be used to generate average signal tracks
BG_BIGWIGS=""
for sample in $BG_SAMPLES; do
    BG_BIGWIGS="$BG_BIGWIGS results/bigwig/${sample}_RPKM.bw"
done

BM_BIGWIGS=""
for sample in $BM_SAMPLES; do
    BM_BIGWIGS="$BM_BIGWIGS results/bigwig/${sample}_RPKM.bw"
done

# Average bigwig files using deepTools to create composite signal tracks
bigwigAverage -b $BG_BIGWIGS -o results/metaprofiles/BG_average.bw
bigwigAverage -b $BM_BIGWIGS -o results/metaprofiles/BM_average.bw

# Function to process gene lists and create bed files with TSS coordinates
# Parameters:
#   $1: input file containing gene names
#   $2: output bed file path
process_gene_list() {
    local input_file=$1
    local output_bed=$2
    local temp_bed="${output_bed}.temp"
    
    # Process each gene and write to temporary file
    while read gene; do
        # Extract gene coordinates from GTF file and create BED format
        zgrep "gene_name \"$gene\"" data/gencode.vM10.annotation.gtf.gz | \
        awk '$3=="gene"' | \
        awk -v OFS="\t" '{
            if($7=="+") {
                # For positive strand: TSS is start position
                tss = $4
            } else {
                # For negative strand: TSS is end position
                tss = $5
            }
            # Clean gene name field and ensure proper BED format
            gsub(/[";]/, "", $10)  # Remove quotes and semicolons from gene name
            # Create window around TSS regardless of strand
            printf "%s\t%d\t%d\t%s\t.\t%s\n", $1, tss-2500, tss+2500, $10, $7
        }' > "$temp_bed"
    done < "$input_file"
    
    # Sort the bed file by chromosome and position, remove duplicates
    sort -k1,1 -k2,2n "$temp_bed" | uniq > "$output_bed"
    rm "$temp_bed"
}

# Process each gene list and generate corresponding metaprofiles
# Lists contain genes with different expression patterns
for list in Gene_lists/by_regulation_type/enriched_down_regulated.csv \
           Gene_lists/by_regulation_type/enriched_not_disregulated.csv \
           Gene_lists/by_regulation_type/enriched_up_regulated.csv; do
    base_name=$(basename $list .csv)
    echo "Processing $base_name..."
    
    # Create TSS bed file for the gene list
    process_gene_list "$list" "results/metaprofiles/${base_name}_TSS.bed"
    
    # Generate matrix of ChIP-seq signal around TSS regions
    computeMatrix reference-point \
        --referencePoint TSS \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --scoreFileName results/metaprofiles/BG_average.bw results/metaprofiles/BM_average.bw \
        --regionsFileName results/metaprofiles/${base_name}_TSS.bed \
        --numberOfProcessors 16 \
        --skipZeros \
        -o results/metaprofiles/${base_name}_matrix.gz
    
    # Create profile plot showing average signal distribution
    plotProfile \
        --matrixFile results/metaprofiles/${base_name}_matrix.gz \
        --outFileName results/metaprofiles/${base_name}_profile.pdf \
        --plotTitle "SMARCB1 binding around TSS - ${base_name}" \
        --regionsLabel "${base_name}" \
        --samplesLabel "BG" "BM" \
        --colors "#1f77b4" "#d62728" \
        --plotHeight 6 \
        --plotWidth 8 \
        --yAxisLabel "Average RPKM"
done

echo "Metaprofile generation completed!"