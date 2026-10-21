# 01_liftover/results (not in git)

- `genotype/`: genotype batches converted from PLINK and merged (hg19)
- `geno_lifted_hg38.vcf.gz` (+ `.tbi`, `.idx`): merged genotypes lifted to GRCh38 (GATK LiftoverVcf)
- `rejected_variants_to_lift.vcf`: variants that could not be lifted
- `geno_lifted_hg38_deduped_renamed.bcf` (+ `.csi`): duplicates removed, IDs set to chr_pos_ref_alt; input of 02
- `geno_metadata.tsv`: genotyped samples
