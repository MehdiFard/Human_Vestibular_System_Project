# 06_twas/gtex_weights: TWAS with GTEx tissue weights (side analysis)

`run_twas_gtex_tissue.smk` runs FUSION TWAS of every GWAS trait with the Pantry weights of the GTEx tissues
(run from this folder: `snakemake -s run_twas_gtex_tissue.smk --profile ../../../profiles/slurm`). Not in git:
- `gwas_sumstats/`: `<trait>_hg38_twas_input.chr_b38.sumstats.gz` (04 TWAS input with chr-prefixed IDs; `scripts/prep_inputs.R`)
  and `gwas_metadata.tsv`
- `weights/`: Pantry TWAS weights of the GTEx tissues (https://github.com/PejLab/Pantry, GTEx v8 pheast output)
- `ld_ref/`: 1000 Genomes LD reference (as resources/ldref_1kg_hg38, chr-prefixed)
- `results/`, `logs/`: FUSION output per tissue, modality and trait
