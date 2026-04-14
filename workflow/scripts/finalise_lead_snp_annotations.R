library(data.table)
library(stringr)

dat <- fread(snakemake@input[[1]], sep = "\t", header = T)
dat[, genes := ""]

dat[missense_gene != "" & qtl_genes == "", genes := missense_gene]
dat[missense_gene == "" & qtl_genes != "", genes := qtl_genes]
dat[missense_gene != "" & qtl_genes != "", genes := mapply(function(missense_gene, qtl_genes) {
  paste(unique(sort(c(missense_gene, unlist(strsplit(qtl_genes, split = ","))))), collapse = ",")
}, missense_gene, qtl_genes, SIMPLIFY = TRUE)]

if (!is.null(snakemake@params$missing_missense_variants)) {
  dat[rsid %in% names(snakemake@params$missing_missense_variants), genes := snakemake@params$missing_missense_variants[rsid]]
}

dat[genes == "" & nearest_gene_name != "", genes := nearest_gene_name]
dat[genes == "" & nearest_gene_name == "", genes := nearest_gene_id]

fwrite(dat, snakemake@output[[1]], sep = "\t")
