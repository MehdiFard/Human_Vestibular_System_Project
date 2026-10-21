# 10 sfFDR: surrogate-functional FDR for the dizziness GWAS

Author: Thomas Willis. Merged with its full history from
[PejLab/sffdr_for_dizziness](https://github.com/PejLab/sffdr_for_dizziness) (October 2026). This folder keeps that
repository's own layout, config, Slurm profile and container.

Steps (`workflow/rules/`):
- `gwas.smk`: lead SNPs of the dizziness/balance meta-analysis (distance clumping, nearest gene, missense and QTL annotation)
- `pred_ld.smk`: PRED-LD imputation of the informative GWAS (tinnitus, hearing, anxiety, depression, frailty) to the
  1000 Genomes EUR panel
- `ldsc.smk`: LDSC genetic correlations with the candidate informative traits, expected chi-squared
- `sffdr.smk`: sfFDR input, sfFDR run (`results/sffdr/results.tsv.gz`), lead SNPs, LocusZoom
- `pub.smk`: publication plots (sfFDR Manhattan, chr18 LocusZoom) and the lead-SNP table (`results/pub/`)

The 1000 Genomes reference comes from the `twillis209/1kGP_pipeline` Snakemake module; tools run in the container
`docker://twillis209/sffdr_for_dizziness` (built from `Dockerfile`).

## Inputs and results

- GWAS summary statistics: `analyses/04_gwas_sumstats/results/formatted_sumstats` (`harmonised_gwas_root` in
  `config/config.yaml`); informative traits and sample sizes are also in `config/config.yaml`.
- `results/` (not in git): the workflow's output (~650 GB for the full run). Fig 2 reads `results/sffdr/results.tsv.gz`;
  it can be a link to an existing run (with `sffdr: rerun: false` in the project config, only its presence is checked).

## Running (compute node or Slurm)

Needs Snakemake 9.25 (`pathvars`) with polars and pandas: conda env `sffdr_smk` (`envs/sffdr_smk.yml`). Set the Slurm
account/partition in `profile/config.yaml`; site-specific values (paths outside the project) go in `config/config_lab.yaml` (not in git), which
`workflow/Snakefile` loads after `config/config.yaml`.

LDlink token: the chr18 LocusZoom plot (`draw_locuszoomr_plot_for_chr18_associations`) colours variants by LD from the
NIH LDlink API, which needs a personal token. Request your own at https://ldlink.nih.gov/?tab=apiaccess and set it as
`apis: ldlink: token:` in `config/config_lab.yaml` (never in `config/config.yaml`, which is in git).

```bash
cd analyses/10_sffdr
snakemake -s workflow/Snakefile --profile profile -n compile_rg draw_locuszoomr_plot_for_chr18_associations \
    draw_sffdr_manhattan_plot_for_publication write_tables_to_xlsx
```

With `results/` linked to the existing run nothing can be written there. To re-run here: `rm results` (the link
only), `mkdir results`, then run the command above without `-n`, or set `sffdr: rerun: true` in the project
`config/config.yaml` and run the main Snakefile.

