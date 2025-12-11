#!/bin/bash
### Pantry - Pheast

#SBATCH --job-name=pheast_vest
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00  
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org

# Activating env to use apptainer
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast

# cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast_no_outlier

# Run snakemake
snakemake --profile slurm \
                --config dataset=vestibular \
                --configfile config.yml \
                --printshellcmds \
                --keep-going \
                --rerun-incomplete \
                -j 20

### END