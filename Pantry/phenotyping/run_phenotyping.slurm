#!/bin/bash
### Running phenotyping by pantry

#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-gdml-investor
#SBATCH --account=cpu-gdml-investor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org

# cpu-core-sponsored
## Go to WD
# cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry


## Activate environment
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry

# ========== Check DATASET ==========
: "${DATASET:?Environment variable DATASET not set. Please submit with --export=ALL,DATASET=your_dataset_name}"

snakemake   --configfile config.yml \
            --config dataset=$DATASET \
            --jobs 34 \
            --default-resources mem_mb=256000 cpus=34 \
            --printshellcmds \
            --latency-wait 60 \
            --rerun-incomplete \
            --keep-going


### END