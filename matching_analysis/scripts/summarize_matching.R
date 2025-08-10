# Making a summary table to look at matched and unmatched strain IDs

# setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/matching_analysis")
rm(list = ls())

library(dplyr)
library(readr)


# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript update_gwas_pos.R <bed_hg38> <gwas_sumstat> <output_file>")
}

mixup_result <- args[1]
metadata <- args[2]
output_file <- args[3]

# Read mix-up data (ASE output)
mix_up <- read.table(mixup_result, header = T, sep = "\t")

# Reading meta data
rna_metadata <- read.table(metadata, header = T, sep = "\t")

# Libraries FTL_500 and FTL_501 are from 1 sample and 1 tissue, so I concatenated them when was mapped them, so merge them now.
indx1 <- which(rna_metadata$library_id == "FTL_500")
indx2 <- which(rna_metadata$library_id == "FTL_501")
rna_metadata[indx1, "library_id"] <- "FTL_500_501"
rna_metadata <- rna_metadata[-indx2, ]

# Merging metadata and mix-up table
mix_meta <- mix_up %>%
    # select(-top_variant_count) %>%
    inner_join(rna_metadata[, c("library_id", "sample_id_geno", "tissue")],
        by = c("RNA_ID" = "library_id"))

# Check the number of matched and unmatched IDs
table(mix_meta$top_simil_ID == mix_meta$sample_id_geno)
not_match <- mix_meta[!mix_meta$top_simil_ID == mix_meta$sample_id_geno, ]

# # Number of not matched ids - unique
# not_match %>% filter(sample_id_geno != "no_genotype") %>%
#     distinct(top_simil_ID) %>%
#     summarise(n())


# Save the results
write_tsv(mix_meta[, c(1,2,8,3,4,6,7,9)], file = output_file)


### END