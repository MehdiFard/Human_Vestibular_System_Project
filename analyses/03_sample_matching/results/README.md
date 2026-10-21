# 03_sample_matching/results (not in git)

- `bam/`: STAR 2.7.11b alignments of every RNA-seq library (`<library>.Aligned.sortedByCoord.out.bam`); `bam_list.txt`
- `genotype/`: imputed genotypes filtered for the matching test (common, exonic SNVs)
- `ase_count/`: GATK ASEReadCounter allele counts per library
- `rna_to_genotype_similarity.tsv`, `rna_to_genotype_summary.tsv`, `mixup_table_forcheck.tsv`: RNA vs genotype
  concordance; each library matched to its donor
- `sample_sheets/` (made by `scripts/make_metadata.R`): `metadata_geno.tsv`, `metadata_rna.tsv`, and per tissue set
  (`vestibular`, `utricle`, `saccule`) `*_covars.tsv` (age, sex) and `*_sample_ids.tsv`. Participant-level.
