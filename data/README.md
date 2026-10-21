# data/ (not in git)

Raw inputs. Individual-level data (genotypes, RNA-seq reads, sample sheets) stay on the analysis server. GWAS summary statistics are public or obtained from the collaborators.

```
data/
  genotype/            genotyping array data, hg19, PLINK text format, one folder per batch (links to the raw data,
                       config.yaml data.genotype_raw_dir):
                         batch12/20240129_hg19.{ped,map}   batches 1-2
                         batch3/20250515_hg19.{ped,map}    batch 3
                         genotype_metadata_GDA_human_samples_b1b2b3.xlsx   sample sheet (sheet 1: batches 1-2, sheet 2: batch 3)
  rna_fq/              paired-end RNA-seq FASTQ of utricle and saccule, one folder per sequencing batch (batch1-3,
                       gold_standard), named <library>_S<n>_L<lane>_R{1,2}_001.fastq.gz (links to the raw data,
                       config.yaml data.fastq_raw_dir); rnaseq_metadata_human_b1b2b3_gs.xlsx = library sheet
  gwas_raw_sumstats/   GWAS summary statistics as received (file names, traits, builds and N in
                       analyses/04_gwas_sumstats/src/gwas_map.tsv)
```

## gwas_raw_sumstats/

| Trait (gwas_map.tsv) | File | Build | Source |
|---|---|---|---|
| balance_meta1 (dizziness, cross-ancestry meta-analysis UKB + MVP) | balance_meta1.tbl.sorted.filtered.gz | GRCh37 | this study (collaborators) |
| balance_meta_eur (dizziness, European ancestry) | balance_meta_EUR.gz | GRCh37 | this study (collaborators) |
| dizziness, UK Biobank (Fig 1) | UKB_EUR_balance.txt.gz | GRCh37 | this study (collaborators) |
| vertigo | vertigo_GRCH38_sorted_cln.tsv.gz | GRCh38 | Skuladottir et al. 2023, Neurobiol Aging (DOI as given in Supp Table S6.1: 10.1016/j.neurobiolaging.2023.03.0…, incomplete) |
| tinnitus_eur | eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma_ours.gz | GRCh37 | MVP + UKB meta-analysis (ndkp.hugeamp.org downloads) |
| hearing_hleur | hearing_HLeur_WellsFlip_STDERR_allchr.hearbroad.maf01.fuma_1.tbl.filtered.gz | GRCh37 | doi:10.1038/s41588-021-00921-z |
| MVP_hearing_loss | All_MVP_Trpchevska_De-Angelis_BBJ.N.txt.gz | GRCh37 | MVP + Trpchevska + De Angelis + BBJ meta-analysis |
| ukbb_falls | UKBB_falls_gwas_MAF001.txt.gz | GRCh37 | UK Biobank falls GWAS (Neurobiol Aging 2023, S0197458023000507) |
| frailty_gwas | Frailty_GWAS_maf0.01_FUMA.tsv | GRCh37 | doi:10.1007/s11357-023-00771-z |
| fluid_intelligence_eur | EUR_fluid_intelligence.sumstats.txt.gz | GRCh37 | GWAS Atlas trait 3413 (PMID 31427789) |
| parkinson_gp2_eur | GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz | GRCh38 | GP2 European meta-analysis 2024 |
| alzheimer | GCST90027158_buildGRCh38.gz.tsv | GRCh38 | GWAS Catalog GCST90027158 (PMID 35379992) |
| depression_mdd | pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv | GRCh37 | PGC MDD 2025 (figshare 27061255) |
| anx_eur | ANX_EUR.txt.gz (+ README_anxiety_disorders.txt) | GRCh37 | anxiety genomic-SEM factor (PMID 39294497) |
| cerebellar / cerebral / subcortical volume | Tissinketal_CommBio2022_*volume.txt.gz | GRCh37 | Tissink et al. 2022, doi:10.1038/s42003-022-03672-7 |
| cognitive_skills_ext, non_cognitive_skills_ext | *cognitive_skills_ext.gz | GRCh37 | cognitive / non-cognitive skills GWAS |
| icp/mcp/scp fa/md | {ICP,MCP,SCP}_{FA,MD}.txt.gz | GRCh37 | cerebellar peduncle diffusion MRI GWAS |

Also here: `magma_all_sets.csv.gz` (MAGMA gene-set results, Fig 2a) and the original files of a few traits
(`balance_meta*.tbl.sorted*`, `vertigo_*`, `GP2_*`) kept as received. Harmonised hg38 versions are made by
`analyses/04_gwas_sumstats`.
