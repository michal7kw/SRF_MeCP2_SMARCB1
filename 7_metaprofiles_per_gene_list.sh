#!/bin/bash
#SBATCH --job-name=7_metaprofiles_per_gene_list
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/7_metaprofiles_per_gene_list.err"
#SBATCH --output="logs/7_metaprofiles_per_gene_list.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create output directory for metaprofiles if not exists
mkdir -p results/metaprofiles

# Function to process gene lists and create BED files with TSS regions
process_gene_list() {
    local input_file=$1
    local output_bed=$2
    local temp_bed="${output_bed}.temp"
    
    # Process each gene and extract TSS regions from GTF
    while read -r gene; do
        # Extract gene entries from GTF and process TSS coordinates
        zgrep "gene_name \"${gene}\"" data/gencode.vM10.annotation.gtf.gz | \
        awk '$3 == "gene"' | \
        awk -v OFS="\t" '{
            if ($7 == "+") {
                # For positive strand: TSS is start position
                print $1, $4-2500, $4+2500, $10, ".", $7
            } else {
                # For negative strand: TSS is end position
                print $1, $5-2500, $5+2500, $10, ".", $7
            }
        }' | tr -d '";' >> "${temp_bed}"
    done < "${input_file}"
    
    # Sort and deduplicate the BED file
    sort -k1,1 -k2,2n "${temp_bed}" | uniq > "${output_bed}"
    rm "${temp_bed}"
}

# Process each gene list and generate TSS BED files
gene_lists=(
    "Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv"
    "Gene_lists/targets/high_expression_targets1_1000.0.csv"
    "Gene_lists/targets/high_expression_targets2_1000.0.csv"
)

for list in "${gene_lists[@]}"; do
    base_name=$(basename "${list}" .csv)
    echo "Processing ${base_name}..."
    
    # Create TSS BED file if it doesn't exist
    if [ ! -f "results/metaprofiles/${base_name}_TSS.bed" ]; then
        echo "Generating ${base_name}_TSS.bed..."
        process_gene_list "${list}" "results/metaprofiles/${base_name}_TSS.bed"
    else
        echo "TSS BED file for ${base_name} exists, skipping."
    fi
    
    # Run computeMatrix only if matrix file doesn't exist
    if [ ! -f results/metaprofiles/${base_name}_matrix.gz ]; then
    echo "Computing matrix for ${base_name}..."
    computeMatrix reference-point \
        --referencePoint TSS \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --scoreFileName results/metaprofiles/BG_average.bw results/metaprofiles/BM_average.bw \
        --regionsFileName results/metaprofiles/${base_name}_TSS.bed \
        --numberOfProcessors 16 \
        --skipZeros \
        -o results/metaprofiles/${base_name}_matrix.gz
    else
        echo "Matrix file for ${base_name} exists, skipping computeMatrix."
    fi
    
    # Run plotProfile only if PDF doesn't exist
    if [ ! -f results/metaprofiles/${base_name}_profile.pdf ]; then
    echo "Generating plot for ${base_name}..."
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
    else
        echo "Profile plot for ${base_name} exists, skipping plotProfile."
    fi
done

echo "Rerun of failed or missing tasks completed!" 