# 02_imputation/results (not in git)

- `phased_vcf/`: genotypes phased with Eagle 2.4.1, per chromosome
- `imputed_vcf/`: imputed with Minimac4 4.1.6 against 1000 Genomes phase 3, per chromosome (filtered on R2 >= config
  `imputation.min_r2`) and merged (`geno_lifted_phased_imputed_hg38_allchrs.vcf.gz`); genotype input of 03 and 05
