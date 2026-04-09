library(dplyr)
library(readr)

save.image('convert.RData')

# Read input
gwas <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# Filter to specific chromosome
chr <- as.integer(gsub("chr", "", snakemake@wildcards$chr))
gwas_chr <- gwas %>%
  filter(!!sym(snakemake@params$chr_col) == chr)

# Convert to PRED-LD format
pred_ld_input <- gwas_chr %>%
  select(
    snp = !!sym(snakemake@params$rsid_col),
    chr = !!sym(snakemake@params$chr_col),
    pos = !!sym(snakemake@params$bp_col),
    A1 = !!sym(snakemake@params$alt_col),  # effect allele
    A2 = !!sym(snakemake@params$ref_col),  # other allele
    beta = !!sym(snakemake@params$beta_col),
    SE = !!sym(snakemake@params$se_col)
  ) %>%
  filter(!is.na(snp), !is.na(beta), !is.na(SE))

# Write output in tab-separated format
write_tsv(pred_ld_input, snakemake@output[[1]])

cat("Converted", nrow(pred_ld_input), "variants for chromosome", chr, "\n")
