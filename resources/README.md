# resources/ (not in git, except gtex_colors.json)

Reference files shared by several analyses and figures (analysis-specific references are in `analyses/*/src/`).

```
resources/
  genome/ensembl110/       Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa (+ .fai, .dict)   Ensembl release 110
                           https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/
  annotation/ensembl110/   Homo_sapiens.GRCh38.110.chr.gtf                                   Ensembl release 110
                           https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/
  ldref_1kg_hg38/          1kg_ldref_hg38_chr{1..22}.{bed,bim,fam}: 1000 Genomes phase 3 (all populations), GRCh38,
                           PLINK; LD reference for FUSION TWAS (06), fine-mapping (07) and the LD panels of Fig 5 / Supp Fig 4
  genome/ensembl110/index_star/  STAR 2.7.11b index of that genome (refs.star_index): RNA-seq mapping for sample matching (03)
  1kgp_phase3_hg38_phased/ ALL.chr{1..22,X,Y}_GRCh38.genotypes.20170504.vcf.gz (+ .tbi): 1000 Genomes phase 3 phased,
                           GRCh38 (refs.kgp_phase3_url); phasing/imputation reference (02)
  1kg_eur_pruned_1000_1_0_2/  merged.prune.{in,out}: LD-pruned 1000 Genomes EUR SNPs (1000 kb window, r2 0.2) from the
                           twillis209/1kGP_pipeline module; sfFDR input (10)
  igsr_samples.tsv         1000 Genomes sample -> population (IGSR data portal); LD by ancestry, Fig 5 / Supp Fig 4
  ucsc_recomb/hg38/        recomb1000GAvg.bw: UCSC hg38 recombination rate (1000 Genomes average), Fig 5b
  gtex_pantry/             GTEx v8 tissues run through Pantry (Munro et al. 2024); only the files read here:
                             phenos/<tissue>/output/unnorm/expression.bed                    (08, Fig 3, Supp Fig 5)
                             phenos/<tissue>/intermediate/expression/<sample>/abundance.tsv  (08: GTEx counts for the PCA)
                             pheast/<tissue>/output/qtl/expression.cis_independent_qtl.txt.gz  (08)
                             pheast/<tissue>/intermediate/twas/                                (06 TWAS with GTEx weights)
                             GTEx/pheast/<tissue>/input/phenotypes/expression.bed.gz           (08)
  tools/                   STAR 2.7.11b, Eagle 2.4.1, Minimac4 4.1.6, LDSC (config tools:)
  gtex_colors.json         GTEx tissue colours (in git)
```
