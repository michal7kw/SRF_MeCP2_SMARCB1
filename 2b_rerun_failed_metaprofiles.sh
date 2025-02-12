#!/bin/bash
#SBATCH --job-name=SMARCB1_rerun_metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="logs/rerun_metaprofiles.err"
#SBATCH --output="logs/rerun_metaprofiles.out"

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create output directory for metaprofiles if not exists
mkdir -p results/metaprofiles

# Define sample groups
BG_SAMPLES="BG1 BG2 BG3"
BM_SAMPLES="BM3"

# Create bigwig file lists for averaging
BG_BIGWIGS=""
for sample in $BG_SAMPLES; do
    BG_BIGWIGS="$BG_BIGWIGS results/bigwig/${sample}_RPKM.bw"
done

BM_BIGWIGS=""
for sample in $BM_SAMPLES; do
    BM_BIGWIGS="$BM_BIGWIGS results/bigwig/${sample}_RPKM.bw"
done

# Compute averaged bigwigs only if not already computed
if [ ! -f results/metaprofiles/BG_average.bw ]; then
    echo "Computing BG_average.bw..."
    bigwigAverage -b $BG_BIGWIGS -o results/metaprofiles/BG_average.bw
else
    echo "BG_average.bw exists, skipping."
fi

if [ ! -f results/metaprofiles/BM_average.bw ]; then
    echo "Computing BM_average.bw..."
    bigwigAverage -b $BM_BIGWIGS -o results/metaprofiles/BM_average.bw
else
    echo "BM_average.bw exists, skipping."
fi

# Function to process gene lists and create bed files with TSS
process_gene_list() {
    local input_file=$1
    local output_bed=$2
    local temp_bed="${output_bed}.temp"
    
    # Process each gene and write to temporary file
    while read gene; do
        zgrep "gene_name \"$gene\"" data/gencode.vM10.annotation.gtf.gz | \
        awk '$3=="gene"' | \
        awk -v OFS="\t" '{
            if($7=="+") {
                print $1, $4-2500, $4+2500, $10, ".", $7
            } else {
                print $1, $5-2500, $5+2500, $10, ".", $7
            }
        }' | tr -d '";' >> "$temp_bed"
    done < "$input_file"
    
    # Sort the bed file and remove potential duplicates
    sort -k1,1 -k2,2n "$temp_bed" | uniq > "$output_bed"
    rm "$temp_bed"
}

# Process each gene list and run subsequent steps conditionally
for list in enriched_down_regulated.csv enriched_not_disregulated.csv enriched_up_regulated.csv; do
    base_name=$(basename $list .csv)
    echo "Processing $base_name..."
    
    # Create TSS bed file if not exists
    if [ ! -f results/metaprofiles/${base_name}_TSS.bed ]; then
        echo "Generating ${base_name}_TSS.bed..."
        process_gene_list $list "results/metaprofiles/${base_name}_TSS.bed"
    else
        echo "TSS bed file for ${base_name} exists, skipping."
    fi
    
    # Run computeMatrix only if matrix file doesn't exist
    # if [ ! -f results/metaprofiles/${base_name}_matrix.gz ]; then
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
    # else
    #     echo "Matrix file for ${base_name} exists, skipping computeMatrix."
    # fi
    
    # Run plotProfile only if PDF doesn't exist
    # if [ ! -f results/metaprofiles/${base_name}_profile.pdf ]; then
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
        --yAxisLabel "Average RPKM"
    # else
    #     echo "Profile plot for ${base_name} exists, skipping plotProfile."
    # fi
done

echo "Rerun of failed or missing tasks completed!" 