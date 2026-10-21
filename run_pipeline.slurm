#!/bin/bash
# Run the pipeline as a Slurm controller job; the profile submits every rule as its own job.
# Before use: set the #SBATCH account/partition below, the <...> values in profiles/slurm/config.yaml, and the paths
# in config/config.yaml (data, references, conda_prefix). Arguments are passed to snakemake:
#   sbatch run_pipeline.slurm                                  # everything (main Snakefile)
#   sbatch run_pipeline.slurm 08_expression_gtex               # one analysis
#   sbatch run_pipeline.slurm -s manuscript/figures/Snakefile fig3     # one figure
#SBATCH --job-name=vestibular_pipeline
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x-%j.log
set -eo pipefail
# environment with Snakemake >= 9 and snakemake-executor-plugin-slurm (envs/pantry_pheast.yml)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV:-pantry_pheast}"
snakemake --profile "$(pwd)/profiles/slurm" "$@"
