library(dplyr)
library(readr)

# Read PRED-LD imputed data
imputed <- read_tsv(snakemake@input$imputed, show_col_types = FALSE)

# Read original GWAS to get sample size and other metadata
original <- read_tsv(snakemake@input$original, show_col_types = FALSE)

# Get sample size from original (assuming it's consistent)
sample_size <- original[[grep("N", names(original), ignore.case = TRUE, value = TRUE)[1]]][1]
if (is.na(sample_size) && "sample_size" %in% names(original)) {
  sample_size <- original$sample_size[1]
}

# PRED-LD output typically has: snp, chr, pos, A1, A2, beta, SE, and possibly p-value
# Convert back to harmonised format
harmonised <- imputed %>%
  mutate(
    # Calculate p-value from beta and SE if not present
    z_score = beta / SE,
    p_value = if ("p_value" %in% names(imputed)) p_value else 2 * pnorm(-abs(z_score))
  ) %>%
  select(
    !!snakemake@params$chr_col := chr,
    !!snakemake@params$bp_col := pos,
    !!snakemake@params$ref_col := A2,
    !!snakemake@params$alt_col := A1,
    !!snakemake@params$rsid_col := snp,
    !!snakemake@params$beta_col := beta,
    !!snakemake@params$se_col := SE,
    !!snakemake@params$p_col := p_value
  )

# Add sample size if available
if (!is.na(sample_size)) {
  harmonised$sample_size <- sample_size
}

# Sort by chromosome and position
harmonised <- harmonised %>%
  arrange(!!sym(snakemake@params$chr_col), !!sym(snakemake@params$bp_col))

# Write output
write_tsv(harmonised, snakemake@output[[1]])

cat("Converted", nrow(harmonised), "variants to harmonised format\n")
cat("Output written to", snakemake@output[[1]], "\n")
