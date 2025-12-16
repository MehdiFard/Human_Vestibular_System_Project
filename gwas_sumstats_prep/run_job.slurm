#!/bin/bash
### Harmonizing GWAS sumstats and making TWAS input

#SBATCH --job-name=harmonizing_sumstats
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-core-sponsored
#SBATCH --nodes=3
#SBATCH --ntasks=22
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org


## Go to WD
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/gwas_sumstats_prep


## Activate environment
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3//etc/profile.d/mamba.sh"
source $MAMBA
mamba activate Env1


# Run Snakemake
snakemake -s Snakefile \
            --default-resources mem_mb=64000 cpus=8 \
            --jobs 22 \
            --printshellcmds \
            --scheduler greedy \
            --latency-wait 60 \
            --rerun-incomplete

### END