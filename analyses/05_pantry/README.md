# 05 Pantry

Multimodal RNA phenotyping and xQTL / heritability / TWAS mapping of human otolith organs with Pantry
(Munro et al. 2024). Pantry's own pipelines run inside `phenotyping/` and `pheast/`; `Snakefile` here runs them in order. TWAS association with the GWAS traits is a separate analysis (`../06_twas`).

```
analyses/05_pantry/
  Snakefile      runs phenotyping -> pheast inputs -> pheast -> xQTL summary tables
  src/           Pantry-specific inputs (editing sites, tx2gene, 107 analysed individuals, GTEx v8 WGS VEP annotation);
                 genome + GTF in ../../resources
  phenotyping/   Pantry phenotyping: 7 modalities (expression, isoforms, splicing, alt TSS, alt polyA, stability,
                 RNA editing) for utricle, saccule and vestibular (= both organs). Snakefile, steps/, scripts/, config.yml.
                 input/ref/ links to src/ and ../../resources (Ensembl 110 genome + GTF).
  pheast/        Pantry pheast: cis-QTL (tensorQTL; independent + nominal), heritability (GCTA), TWAS weights (FUSION).
                 prepare_inputs.sh builds input/ from ../phenotyping/output (see below); input/LDREF -> ../../resources/ldref_1kg_hg38.
  scripts/       xQTL summary tables for Fig 4 / Supp Fig 3 (from the Pantry paper's figure scripts):
                 xqtl_tables.R (phenotypes per gene, top associations, cross-modality xQTLs), xqtl_positions.R
                 (position along the gene), xqtl_func_enrichment.R (functional annotation enrichment)
  results/xqtl_summary/   their output (protein-coding and lncRNA genes)
```

## RNA editing in this pipeline

- **0/0 fix** (`phenotyping/scripts/assemble_bed.py`): samples below the coverage cutoff (`0/0` in the edit matrix)
  are set to missing instead of 100% edited, so the >40%-missing filter and mean imputation apply.
- **Host-gene expression regressed out** (`pheast/scripts/regress_rna_editing_host_expression.py`, Li et al. 2022):
  each site's normalised editing level is regressed on log2(TPM+1) of its host gene(s); residuals are rank-INT
  normalised. Cross-modality is rebuilt with the regressed editing (`pheast/scripts/build_cross_modality.py`).
  The unregressed inputs stay in `pheast/input/vestibular/phenotypes_unregressed/`.
- TWAS weights of RNA editing were computed after the FUSION sign fix (`scripts/fusion_twas/FUSION.compute_weights.R`),
  those of the other modalities before it; `../06_twas/scripts/annotate_twas_hits.R` corrects the sign of the latter.

## Run (compute node, never the login node)

```bash
# from the project root, all stages (phenotyping -> pheast inputs -> pheast on a GPU node -> xQTL summary tables):
snakemake -s analyses/05_pantry/Snakefile --profile profiles/slurm
# or by hand, inside analyses/05_pantry:
cd phenotyping && snakemake --config dataset=vestibular --profile ../../../profiles/slurm && cd ..
cd pheast && bash prepare_inputs.sh vestibular     # cross-modality, 107 individuals, editing regression, genotypes, covariates
MAMBA_SH=<mamba.sh> sbatch -p <gpu_partition> -A <account> --gres=gpu:1 -c 34 --mem=256G -t 48:00:00 \
  --wrap "bash run_pheast_gpu.sh config.yml vestibular"   # pantry_pheast env (tensorQTL, torch, R + qvalue)
```
TWAS association, coloc and hits: `analyses/06_twas`.

Final results: `pheast/output/vestibular/{qtl,heritability,twas}` (TWAS weights; association in `../06_twas`);
summary tables for the figures in `results/xqtl_summary/`.
Dataset `vestibular` is the one reported (n = 107 donors); utricle and saccule were run as well.

FUSION (`pheast/scripts/fusion_twas/`, https://github.com/gusevlab/fusion_twas): its two binaries `gemma` and
`gcta_nr_robust` are not in git; copy them from the FUSION repository into that folder.
