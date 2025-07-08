#!/bin/bash
### Mapping by STAR

#SBATCH --job-name=mapping
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org


## Go to WD
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry

## Activate environment
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate Env1


# Run Snakemake
snakemake   --cores 34 \
            --default-resources mem_mb=128000 cpus=34 \
            --printshellcmds \
            --latency-wait 60 \
            --scheduler greedy


# --use-conda \
# --rerun-incomplete \
# --keep-going \

### END
