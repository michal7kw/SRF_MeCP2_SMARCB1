#!/bin/bash
#SBATCH --job-name=BM2
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM2.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM2.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# 1. First, run bowtie2 alignment for BM3
bowtie2 -p 32 -x /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/mm10/mm10 \
    -1 failed_data/BM2_R1_001.fastq.gz \
    -2 failed_data/BM2_R2_001.fastq.gz | \
samtools sort -@ 32 -o results/bowtie2_alt/BM2.sorted.bam - 

# 2. Index the BAM file
samtools index results/bowtie2_alt/BM2.sorted.bam

# 3. Call peaks using MACS2
macs2 callpeak \
    -t results/bowtie2_alt/BM2.sorted.bam \
    -f BAMPE \
    -g mm \
    -n BM2 \
    --outdir results/peaks_alt \
    --nomodel \
    --extsize 200