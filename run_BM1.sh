#!/bin/bash
#SBATCH --job-name=BM1
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM1.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/BM1.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# 1. First, run bowtie2 alignment for BM3
bowtie2 -p 32 -x /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/mm10/mm10 \
    -1 failed_data/BM1_R1_001.fastq.gz \
    -2 failed_data/BM1_R2_001.fastq.gz | \
samtools sort -@ 32 -o results/bowtie2_alt/BM1.sorted.bam - 

# 2. Index the BAM file
samtools index results/bowtie2_alt/BM1.sorted.bam

# 3. Call peaks using MACS2
macs2 callpeak \
    -t results/bowtie2_alt/BM1.sorted.bam \
    -f BAMPE \
    -g mm \
    -n BM1 \
    --outdir results/peaks_alt \
    --nomodel \
    --extsize 200