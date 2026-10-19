library(dplyr)
library(readr)

# Read all chromosome files
chr_files <- snakemake@input
cat("Merging", length(chr_files), "chromosome files\n")

# Read and combine all files
all_data <- lapply(chr_files, function(f) {
  cat("Reading", f, "\n")
  read_tsv(f, show_col_types = FALSE)
}) %>%
  bind_rows()

cat("Total variants after merging:", nrow(all_data), "\n")

# Write compressed output
write_tsv(all_data, snakemake@output[[1]])

cat("Merged output written to", snakemake@output[[1]], "\n")
