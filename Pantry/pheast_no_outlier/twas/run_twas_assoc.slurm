#!/bin/bash
### Running TWAS association test

#SBATCH --job-name=twas_assoc
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org


## Go to WD
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast_no_outlier/twas/


## Activate environment
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3//etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast


# Run Snakemake
snakemake -s run_twas_with_tissue.smk \
            --default-resources mem_mb=32000 cpus=8 \
            --jobs 10 \
            --printshellcmds \
            --rerun-incomplete \
            --keep-going \
            --latency-wait 60

### END