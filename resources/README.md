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
  ucsc_recomb/hg38/        recomb1000GAvg.bw: UCSC hg38 recombination rate (1000 Genomes average), Fig 5b
  gtex_colors.json         GTEx tissue colours (in git)
```
