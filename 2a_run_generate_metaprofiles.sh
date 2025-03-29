#!/bin/bash

# This script generates metaprofiles for ChIP-seq data around transcription start sites (TSS)
# It processes gene lists to create TSS bed files and generates matrices and plots showing ChIP-seq signal around TSS
# Note: This script requires average signal tracks to be generated first using 1d_create_average_tracks.sh

#SBATCH --job-name=2a_run_generate_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-2
#SBATCH --error="logs/2a_run_generate_metaprofiles%a.err"
#SBATCH --output="logs/2a_run_generate_metaprofiles%a.out"

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

# Check if average bigWig files exist and are valid
if ! check_bigwig "results/metaprofiles/BG_average.bw" || ! check_bigwig "results/metaprofiles/BM_average.bw"; then
    echo "Error: Average bigWig files not found or invalid. Please run 1d_create_average_tracks.sh first."
    exit 1
fi

# Function to process gene lists and create bed files with TSS coordinates
process_gene_list() {
    local input_file=$1
    local output_bed=$2
    local temp_bed="${output_bed}.temp"
    
    # Check if input file exists
    if [ ! -f "$input_file" ]; then
        echo "Error: Input file $input_file not found"
        exit 1
    }
    
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
    
    # Check if output bed file was created successfully
    if [ ! -s "$output_bed" ]; then
        echo "Error: Failed to create BED file $output_bed"
        exit 1
    fi
}

# Define array of gene lists
GENE_LISTS=(
    "Gene_lists/by_regulation_type/enriched_down_regulated.csv"
    "Gene_lists/by_regulation_type/enriched_not_disregulated.csv"
    "Gene_lists/by_regulation_type/enriched_up_regulated.csv"
)

# Process the gene list corresponding to the array task ID
list="${GENE_LISTS[$SLURM_ARRAY_TASK_ID]}"
base_name=$(basename $list .csv)
echo "Processing $base_name..."

# Create TSS bed file for the gene list
process_gene_list "$list" "results/metaprofiles/${base_name}_TSS.bed"

# Generate matrix of ChIP-seq signal around TSS regions
echo "Generating matrix for $base_name..."
computeMatrix reference-point \
    --referencePoint TSS \
    --beforeRegionStartLength 2500 \
    --afterRegionStartLength 2500 \
    --scoreFileName results/metaprofiles/BG_average.bw results/metaprofiles/BM_average.bw \
    --regionsFileName results/metaprofiles/${base_name}_TSS.bed \
    --numberOfProcessors 16 \
    --skipZeros \
    -o results/metaprofiles/${base_name}_matrix.gz

# Check if matrix was created successfully
if [ ! -s "results/metaprofiles/${base_name}_matrix.gz" ]; then
    echo "Error: Failed to create matrix file for $base_name"
    exit 1
fi

# Create profile plot showing average signal distribution
echo "Creating profile plot for $base_name..."
plotProfile \
    --matrixFile results/metaprofiles/${base_name}_matrix.gz \
    --outFileName results/metaprofiles/${base_name}_profile.pdf \
    --plotTitle "SMARCB1 binding around TSS - ${base_name}" \
    --regionsLabel "${base_name}" \
    --samplesLabel "BG" "BM" \
    --colors "#1f77b4" "#d62728" \
    --plotHeight 6 \
    --plotWidth 8 \
    --yAxisLabel "Average CPM"

echo "Metaprofile generation completed for $base_name!"