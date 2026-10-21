#!/bin/bash
# Run pheast (QTL mapping, heritability, TWAS weights) on one GPU node with the pantry_pheast environment
# (tensorQTL, torch, R + qvalue). Called by ../Snakefile (rule pheast), which sets the variables below from the config.
# Optional: PHEAST_CONTAINER = an Apptainer image to run inside (site-specific, e.g. when the host system libraries do
# not work with the environment); empty = run directly in the environment.
# By hand:
#   MAMBA_SH=<miniforge>/etc/profile.d/mamba.sh sbatch --gres=gpu:1 -c 34 --mem=256G -t 48:00:00 \
#     --wrap "bash analyses/05_pantry/pheast/run_pheast_gpu.sh"
#   bash run_pheast_gpu.sh [configfile=config.yml] [dataset=vestibular] [extra snakemake args / targets]
set -eo pipefail
cd "$(dirname "$0")"
CFG=${1:-config.yml}; DS=${2:-vestibular}; shift 2 2>/dev/null || shift $#
MAMBA=${MAMBA_SH:?set MAMBA_SH (miniforge mamba.sh; config mamba_sh)}
ENV=${PHEAST_ENV:-pantry_pheast}
RUN="snakemake --configfile $CFG --config dataset=$DS --jobs 34 --default-resources mem_mb=256000 cpus=34 \
     --resources gpu=1 --printshellcmds --latency-wait 60 --keep-going $*"
echo "pheast ($CFG, $DS) start $(date)"
source "$MAMBA"; mamba activate "$ENV"
if [ -n "${PHEAST_CONTAINER:-}" ]; then
  apptainer exec --nv ${APPTAINER_BIND_DIR:+--bind "$APPTAINER_BIND_DIR"} "$PHEAST_CONTAINER" bash -c "
    source $MAMBA && mamba activate $ENV && cd $PWD && (nvidia-smi -L || true) && $RUN"
else
  (nvidia-smi -L || true); eval "$RUN"
fi
echo "pheast done $(date)"
