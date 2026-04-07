library(data.table)

lead_snps <- fread(snakemake@input$lead_snps, header = T, sep = "\t")
sffdr_results <- fread(snakemake@input$sffdr, header = T, sep = "\t")

merged <- merge(lead_snps, sffdr_results[, -c("rsid", "p_value", "prune_in", "fp", "beta", "standard_error")], by = c("chromosome", "base_pair_location", "other_allele", "effect_allele"), all.x = TRUE)

fwrite(merged, file = snakemake@output[[1]], sep = "\t")
