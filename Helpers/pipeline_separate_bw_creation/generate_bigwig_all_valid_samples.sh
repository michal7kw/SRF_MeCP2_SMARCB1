#!/bin/bash
#SBATCH --job-name=bigwig_array
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=1-4
#SBATCH --error="logs/bigwig_%A_%a.err"
#SBATCH --output="logs/bigwig_%A_%a.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1"
cd ${WORKING_DIR}

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define sample names array
declare -a SAMPLES=(
    "BG1"
    "BG2"
    "BG3"
    "BM3"
)

# Get current sample name from array index
# SLURM_ARRAY_TASK_ID starts from 1, so we subtract 1 for zero-based array indexing
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing sample: $SAMPLE"

# Define directories
RESULTS_DIR="results"
BIGWIG_DIR="${RESULTS_DIR}/bigwig"
mkdir -p ${BIGWIG_DIR}

# Generate chromosome sizes file if not already present
CHROM_SIZES="${WORKING_DIR}/mm10/mm10.chrom.sizes"
if [ ! -f "$CHROM_SIZES" ]; then
    echo "Generating chromosome sizes file..."
    fetchChromSizes mm10 > ${CHROM_SIZES}
fi

echo "Converting BAM to bigWig for ${SAMPLE}..."

# Generate RPKM-normalized bigWig file
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${BIGWIG_DIR}/${SAMPLE}_RPKM.bw \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates

# Generate CPM-normalized bigWig file
bamCoverage \
    --bam ${RESULTS_DIR}/filtered/${SAMPLE}.dedup.bam \
    --outFileName ${BIGWIG_DIR}/${SAMPLE}_CPM.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates

echo "BigWig generation completed successfully for ${SAMPLE}!"