library(data.table)

lead_snps <- fread(snakemake@input$lead_snps, header = T, sep = "\t")
sffdr_results <- fread(snakemake@input$sffdr, header = T, sep = "\t")

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]
ref_col <- snakemake@config[["col_names"]][["other_allele"]]
alt_col <- snakemake@config[["col_names"]][["effect_allele"]]
p_col <- snakemake@config[["col_names"]][["p_value"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]

sffdr_results[, c(rsid_col, p_col, "prune_in", "fp") := NULL]
sffdr_results[, (chr_col) := as.character(get(chr_col))]
sffdr_results[, (ref_col) := as.character(get(ref_col))]
sffdr_results[, (alt_col) := as.character(get(alt_col))]

lead_snps[, (chr_col) := as.character(get(chr_col))]
lead_snps[, (ref_col) := as.character(get(ref_col))]
lead_snps[, (alt_col) := as.character(get(alt_col))]

lead_snps[, genes := ""]

lead_snps[missense_gene != "" & qtl_genes == "", genes := missense_gene]
lead_snps[missense_gene == "" & qtl_genes != "", genes := qtl_genes]
lead_snps[missense_gene != "" & qtl_genes != "", genes := mapply(function(missense_gene, qtl_genes) {
  paste(unique(sort(c(missense_gene, unlist(strsplit(qtl_genes, split = ","))))), collapse = ",")
}, missense_gene, qtl_genes, SIMPLIFY = TRUE)]

lead_snps[genes == "" & nearest_gene_name != "", genes := nearest_gene_name]
lead_snps[genes == "" & nearest_gene_name == "", genes := nearest_gene_id]

merged <- merge(lead_snps, sffdr_results, by = c(chr_col, bp_col, ref_col, alt_col), all.x = TRUE)

fwrite(merged, file = snakemake@output[[1]], sep = "\t")
