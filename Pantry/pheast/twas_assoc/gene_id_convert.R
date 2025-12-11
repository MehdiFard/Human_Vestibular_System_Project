### Adding gene symbols to the TWAS output
rm(list = ls())

library(data.table)
library(rtracklayer)
library(dplyr)

# Reading TWAS hits
twas_hits <- fread("Pantry/pheast/twas_assoc/vestibular/output/expression/twas_hits.tsv")
twas_hits$ENSG <- sub("_.*$", "", twas_hits$ID)

# Import GTF as a GRanges object
gtf <- import("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# Keep only gene entries
genes_gtf <- gtf[gtf$type == "gene"]

# Extract metadata columns into a data.frame
anno_gtf <- mcols(genes_gtf) %>%
    as.data.frame() %>%
    transmute(
        gene_id    = sub("\\..*$", "", gene_id),  # drop version .1, .2, ...
        gene_name  = gene_name,
        gene_type  = if ("gene_biotype" %in% names(.)) gene_biotype else NA_character_
    ) %>%
    distinct(gene_id, .keep_all = TRUE)

# Merging
twas_annot <- twas_hits %>%
    left_join(anno_gtf, by = c("ENSG" = "gene_id"), relationship = "many-to-many") %>%
    relocate(ENSG, gene_name, gene_type, .after = ID)

# Inspect
twas_annot %>%
    select(ID, ENSG, gene_name, gene_type, TWAS.Z, TWAS.P) %>%
    head()

fwrite(twas_annot, "Pantry/pheast/twas_assoc/vestibular/twas_hits_annot.tsv", sep = "\t", row.names = F, quote = F)


# ###-------------------------- Using org.Hs.eg.db
# library(org.Hs.eg.db)
# library(AnnotationDbi)

# # Extracting IDs
# gene_map <- AnnotationDbi::select(
#     org.Hs.eg.db,
#     keys       = twas_hits$ENSG,
#     columns    = c("SYMBOL", "GENENAME", "GENETYPE", "ENSEMBL"),
#     keytype    = "ENSEMBL"
# )

# # Merging
# twas_annot <- twas_hits %>%
#     left_join(gene_map, by = c("ENSG" = "ENSEMBL"), relationship = "many-to-many")


### END