### Updating GWAS sumstats SNP positions to hg38
setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast/lift_gwas")
rm(list = ls())

library(data.table)
library(dplyr)
library(readr)
library(stringr)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript update_gwas_pos.R <bed_hg38> <gwas_sumstat> <output_file>")
}

bed_file <- args[1]
gwas_file <- args[2]
output_file <- args[3]

# Read inputs
bed_hg38 <- fread(bed_file)
names(bed_hg38) <- c("chr", "rsID", "pos_hg38", "ref", "alt")

gwas_sumstat <- fread(gwas_file)

# Join and update
gwas_updated_pos <- gwas_sumstat %>%
    inner_join(bed_hg38, by = c("MarkerName" = "rsID")) %>%
    mutate(chr_pos_ref_alt = str_c(Chromosome, pos_hg38, ref, alt, sep = "_")) %>%
    rename(rsID = MarkerName, MarkerName = chr_pos_ref_alt) %>%
    select(Chromosome, Position, MarkerName, everything(), -chr, -rsID)

# Write output
write_tsv(gwas_updated_pos, output_file, col_names = TRUE)


### END
