#!/bin/bash
### DEG

#SBATCH --job-name=DEA_sub
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=cpu-gdml-investor
#SBATCH --account=cpu-gdml-investor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mehdi.esmaeilifard@seattlechildrens.org

# cpu-core-sponsored
## Go to WD
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry


## Activate environment
MAMBA="/data/hps/assoc/private/gdml/user/sesmae/miniforge3/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate Env1

# Run DGE
Rscript Pantry/phenotyping/scripts/differential_gene_expression.R