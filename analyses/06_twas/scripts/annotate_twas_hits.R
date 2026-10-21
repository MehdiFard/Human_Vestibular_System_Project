### Annotate the combined TWAS hits (<dataset>/output/twas_hits.tsv): TWAS.Z sign correction + gene symbols.
### FUSION.compute_weights.R read genotype dosages for the wrong allele (A2) for the enet and top1 models; the script
### was fixed (scripts/fusion_twas/FUSION.compute_weights.R: bed <- 2 - bed). Weights computed before the fix give the
### opposite TWAS.Z sign for enet/top1, so those rows are flipped for the modalities listed as "old-weight" (in this
### run: all modalities except RNA_editing, whose weights were recomputed after the fix). For a run where all
### weights were made with the fixed script, pass "none".
### Usage (from analyses/06_twas):
###   Rscript scripts/annotate_twas_hits.R [dataset=vestibular] [old_weight_modalities=alt_polyA,alt_TSS,expression,isoforms,splicing,stability] [out]
suppressPackageStartupMessages({ library(data.table); library(rtracklayer); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[1] else "vestibular"
old_mods <- if (length(args) >= 2) args[2] else "alt_polyA,alt_TSS,expression,isoforms,splicing,stability"
old_mods <- if (old_mods == "none") character(0) else strsplit(old_mods, ",")[[1]]
out_file <- if (length(args) >= 3) args[3] else file.path("results", dataset, "output", "twas_hits_annot.tsv")

twas_hits <- fread(file.path("results", dataset, "output", "twas_hits.tsv"))
flip <- twas_hits$MODALITY %in% old_mods & twas_hits$MODEL %in% c("enet", "top1")
twas_hits$TWAS.Z[flip] <- -twas_hits$TWAS.Z[flip]
cat(sprintf("%d hits; TWAS.Z flipped for %d enet/top1 hits of old-weight modalities (%s)\n",
            nrow(twas_hits), sum(flip), paste(old_mods, collapse = ", ")))
twas_hits$ENSG <- sub("_.*$", "", twas_hits$ID)

gtf <- import("../../resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf")
anno_gtf <- mcols(gtf[gtf$type == "gene"]) %>%
    as.data.frame() %>%
    transmute(gene_id = sub("\\..*$", "", gene_id), gene_name = gene_name,
              gene_type = if ("gene_biotype" %in% names(.)) gene_biotype else NA_character_) %>%
    distinct(gene_id, .keep_all = TRUE)
twas_annot <- twas_hits %>%
    left_join(anno_gtf, by = c("ENSG" = "gene_id"), relationship = "many-to-many") %>%
    relocate(ENSG, gene_name, gene_type, .after = ID)
fwrite(twas_annot, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("hits by modality:\n"); print(table(twas_annot$MODALITY))
cat("written:", out_file, "\n")
