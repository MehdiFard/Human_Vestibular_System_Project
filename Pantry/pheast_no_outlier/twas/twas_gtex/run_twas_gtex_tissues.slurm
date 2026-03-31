#!/bin/bash

# TWAS - GTEx Tissues

#SBATCH --job-name=twas_gtex
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00  
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org


# Activating env
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast


# WD
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast_no_outlier/twas/twas_gtex


# Run Snakemake
snakemake -s run_twas_gtex_tissue.smk \
                --profile slurm \
                --jobs 250 \
                --keep-going \
                --rerun-incomplete

### END