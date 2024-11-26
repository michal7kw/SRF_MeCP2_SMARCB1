#!/bin/bash
#SBATCH --job-name=Azenta
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/snakemake.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta/logs/snakemake.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

cd /beegfs/scratch/ric.broccoli/kubacki.michal/Azenta

# Create logs directory
mkdir -p logs/slurm

# Clear snakemake cache
rm -rf .snakemake

# Unlock the working directory in case it's locked
# snakemake --unlock

# First do a dry run to check what will be executed
# snakemake -n -p

# Force rerun all rules by removing all output files and rerunning
snakemake --cores 32 --forceall

# snakemake --rerun-incomplete --cores 32