## One-time cache of the ZNF91 region (chr19:22.9-23.9 Mb) for fast figure re-renders:
## GWAS (balance_meta1, chr19 region) + nominal QTL pairs of the ZNF91 phenotypes (7 modalities; only needed columns).
## Output: data/cache_znf91_region.rds (summary statistics only).
suppressPackageStartupMessages({ library(arrow); library(data.table); library(dplyr) })
t0 <- Sys.time()
OUT <- "manuscript/figures/Fig5/data/cache_znf91_region.rds"
REG <- c(22.9e6, 23.9e6); GENE_ENSG <- "ENSG00000167232"
gw <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($3==19 && $4>=%d && $4<=%d)'",
                          "analyses/04_gwas_sumstats/results/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz", REG[1], REG[2]))
NOM <- "%s/output/vestibular/qtl/nominal/%s/vestibular_%s.cis_qtl_pairs.19.parquet"
mods <- c("isoforms", "expression", "alt_polyA", "alt_TSS", "RNA_editing", "splicing", "stability")
qtl <- lapply(mods, function(m) {
  root <- if (m == "RNA_editing") "analyses/05_pantry/pheast" else "analyses/05_pantry/pheast"   # RNA editing: host-gene expression regressed out (2026-10-14)
  d <- open_dataset(sprintf(NOM, root, m, m)) %>% select(phenotype_id, variant_id, pval_nominal) %>% collect()
  d$pos <- as.integer(sub("^[^_]+_([0-9]+)_.*$", "\\1", d$variant_id))
  d <- d[d$pos >= REG[1] & d$pos <= REG[2] & d$pval_nominal > 0, ]
  d$modality <- m
  d[grepl(paste0("^", GENE_ENSG), d$phenotype_id) | TRUE, ]          # keep all phenotypes in region (fallback choice)
})
names(qtl) <- mods
saveRDS(list(gwas = gw, qtl = qtl, region = REG), OUT)
cat(sprintf("cached: GWAS %d rows; QTL rows %s | %.1f s\n", nrow(gw),
            paste(sapply(qtl, nrow), collapse = "/"), as.numeric(difftime(Sys.time(), t0, units = "secs"))))


### END