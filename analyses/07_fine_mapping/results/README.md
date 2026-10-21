# 07_fine_mapping/results (not in git)

SuSiE credible sets and PIPs per locus (`susie_results_<GWAS>_<chr>_<region>.tsv`) for the dizziness GWAS loci, and for
the ZNF91 isoform / cross-modality QTLs.

`*_old.tsv`: results before the October 2026 update (made from the pheast run before the RNA editing 0/0 fix and host-gene
expression regression; the cross-modality QTL fine-mapping changed: lead rs426564 -> rs607168). The GWAS-locus results
reproduce within floating-point precision. UKB chr13 has no 95% credible set (empty table).
The nominal isoform-QTL fine-mapping (sus_isoQTL_chr19.R) finds no credible set; it writes all-SNP PIPs (*_all_snps.tsv).
