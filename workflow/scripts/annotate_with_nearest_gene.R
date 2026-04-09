library(GenomicRanges)
library(ensembldb)

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]

edb <- EnsDb(snakemake@input$edb)

daf <- read.table(snakemake@input$lead, header = TRUE, sep = '\t')

daf[[chr_col]] <- as.character(daf[[chr_col]])

daf[daf[[chr_col]] == '23', chr_col] <- 'X'

snps_gr <- GRanges(
  seqnames = daf[[chr_col]],
  ranges = IRanges(start = daf[[bp_col]], width = 1),
  rsid = daf[[rsid_col]]
)

genes_gr <- genes(edb, filter = GeneBiotypeFilter("protein_coding"))

suppressWarnings(nearest_hits <- distanceToNearest(snps_gr, genes_gr))

nearest_genes <- data.frame(
  rsid = mcols(snps_gr[queryHits(nearest_hits)])$rsid,
  nearest_gene_id        = names(genes_gr[subjectHits(nearest_hits)]),
  nearest_gene_name      = genes_gr[subjectHits(nearest_hits)]$gene_name,
  distance_bp    = mcols(nearest_hits)$distance
)

names(nearest_genes)[names(nearest_genes) == 'rsid'] <- rsid_col

merged <- merge(daf, nearest_genes, by = rsid_col, all.x = TRUE)

merged[merged[[chr_col]] == 'X', chr_col] <- '23'

write.table(merged, file = snakemake@output[[1]], sep = '\t', row.names = FALSE, quote = FALSE)
