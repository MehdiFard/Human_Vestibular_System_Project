#!/bin/bash
### Pantry - Pheast - wrapper

#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=gpu-core-sponsored
#SBATCH --account=gpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --gres=gpu:1
#SBATCH --mem=256G
#SBATCH --time=48:00:00  
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org

# Activating env to use apptainer
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate apptainer_env

# Activating docker container and running the scripts
CONTAINER="/data/hps/assoc/private/posit_workbench/container/posit-base-20250410-ppm-testing.sif"
apptainer exec --nv --bind /data/hps/assoc/ $CONTAINER ./run_pheast_gpu.sh $DATASET


### END