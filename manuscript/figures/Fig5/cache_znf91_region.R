## ZNF91 region (chr19:22.9-23.9 Mb) extracted once for Fig 5 and Supp Fig 4 (first step; the plotting scripts read
## only this file): GWAS (balance_meta1) + nominal QTL pairs of all phenotypes in the region (7 modalities; only the
## needed columns). Output: data/cache_znf91_region.rds (summary statistics only). Run from the project root:
##   Rscript manuscript/figures/Fig5/cache_znf91_region.R
suppressPackageStartupMessages({ library(arrow); library(data.table); library(dplyr) })
t0 <- Sys.time()
OUT <- "manuscript/figures/Fig5/data/cache_znf91_region.rds"
REG <- c(22.9e6, 23.9e6)
gw <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($3==19 && $4>=%d && $4<=%d)'",
                          "analyses/04_gwas_sumstats/results/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz", REG[1], REG[2]))
NOM <- "%s/output/vestibular/qtl/nominal/%s/vestibular_%s.cis_qtl_pairs.19.parquet"
mods <- c("isoforms", "expression", "alt_polyA", "alt_TSS", "RNA_editing", "splicing", "stability")
qtl <- lapply(mods, function(m) {   # RNA editing phenotypes: host-gene expression regressed out
  d <- open_dataset(sprintf(NOM, "analyses/05_pantry/pheast", m, m)) %>% select(phenotype_id, variant_id, pval_nominal) %>% collect()
  d$pos <- as.integer(sub("^[^_]+_([0-9]+)_.*$", "\\1", d$variant_id))
  d <- d[d$pos >= REG[1] & d$pos <= REG[2] & d$pval_nominal > 0, ]
  d$modality <- m
  d
})
names(qtl) <- mods
saveRDS(list(gwas = gw, qtl = qtl, region = REG), OUT)
cat(sprintf("cached: GWAS %d rows; QTL rows %s | %.1f s\n", nrow(gw),
            paste(sapply(qtl, nrow), collapse = "/"), as.numeric(difftime(Sys.time(), t0, units = "secs"))))


### END