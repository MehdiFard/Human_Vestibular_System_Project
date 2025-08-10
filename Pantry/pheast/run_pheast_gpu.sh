#!/bin/bash
# Running Pantry - Pheast on GPU
DATASET=$1

# Activating conda env in the docker container
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3//etc/profile.d/mamba.sh"
source $MAMBA
mamba activate pantry_pheast

# Start time
echo "Starting workflow at: $(date)"

# Check
echo "==== GPU Info ===="
nvidia-smi
echo ""

echo "==== CUDA Compiler Version ===="
nvcc --version
echo ""

echo "==== Running Snakemake workflow ===="
echo ""

# Running snakemake workflow
snakemake   --configfile config.yml \
            --config dataset=$DATASET \
            --jobs 34 \
            --default-resources mem_mb=256000 cpus=34 \
            --resources gpu=1 \
            --printshellcmds \
            --latency-wait 60 \
            --keep-going


# --keep-going
# -R tensorqtl_trans
# --rerun-incomplete


# Finish time
echo "Finished workflow at: $(date)"


### END