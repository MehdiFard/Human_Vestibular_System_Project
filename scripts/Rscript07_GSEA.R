setwd("/active/mohammadi_p/mehdi_projects/vestibular_mouse")
rm(list = ls())

library(data.table)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(enrichplot)
library(pathview)
library(GOSemSim)
library(OrganismDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Reading data
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
sac_vs_utr_w8 <- read.csv("differential_expression/DEA_result_utricle_vs_saccule_in_w8.csv")[, c(1,3,7)]
sac_vs_utr_w8_up <- sac_vs_utr_w8 %>% 
  filter(log2FoldChange > 2 & padj < 0.01) %>%
  pull(gene_id) %>% 
  sort(., decreasing = T)

gene_ids <- bitr(sac_vs_utr_w8_up, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)



# GO over-representation analysis
sac_vs_utr_w8_up_GOs <- enrichGO(gene         = sac_vs_utr_w8_up,
                                OrgDb         = org.Mm.eg.db,
                                keyType       = 'ENSEMBL',
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE,
                                pool          = TRUE)

head(sac_vs_utr_w8_up_GOs)

# Getting gene-GO map file
mmGO_BP <- godata('org.Mm.eg.db', ont="BP", keytype = "ENSEMBL")


# GO Gene Set Enrichment Analysis
sac_vs_utr_w8_up_GSEA <- gseGO(geneList     = sac_vs_utr_w8_up,
              OrgDb        = org.Mm.eg.db,
              keyType      = 'ENSEMBL',
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 3,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = T)
# pAdjustMethod = "none"


# KEGG pathway gene set enrichment analysis
search_kegg_organism('mmu', by='kegg_code')
# KEGG pathway analysis
sac_vs_utr_w8_up_KEGG <- gseKEGG(geneList     = sac_vs_utr_w8_up,
                         organism             = 'mmu',
                         keyType              = "kegg",
                         minGSSize            = 5,
                         pvalueCutoff         = 0.05,
                         verbose              = FALSE)


### Let's us gprofiler
# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
library(gprofiler2)
sac_vs_utr_w8_up_gprofiler <- gost(query = sac_vs_utr_w8_up, organism = "mmusculus",
                significant = T, correction_method = "bonferroni",
                user_threshold = 0.05)
# organism = "hsapiens"

# plot
gostplot(sac_vs_utr_w8_up_gprofiler, capped = FALSE, interactive = FALSE)



### Using enrichR
# https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html
# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
library(enrichR)

sac_vs_utr_w8 <- read.csv("differential_expression/DEA_result_saccule_vs_utricle_in_w22.csv")[, c(1,3,7)]
sac_vs_utr_w8_up <- sac_vs_utr_w8 %>% 
  filter(log2FoldChange > 2 & padj < 0.01) %>%
  pull(gene_id) %>% 
  sort(., decreasing = T)

gene_ids <- bitr(sac_vs_utr_w8_up, fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Mm.eg.db)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
## if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2019_Mouse")
if (websiteLive) {
  enriched <- enrichr(gene_ids$SYMBOL, dbs)
}

# Plot KEGG results
jpeg("plots/10_KEGG_utr_vs_sac_w8_ups.jpeg", width = 20, height = 20, units = "cm", res = 300,)
if (websiteLive) {
  plotEnrich(enriched[[4]], showTerms = 20, numChar = 60, y = "Count",
             orderBy = "P.value", title = "Utricle vs Saccule in W8 - Mm - UPs - KEGG")
}
dev.off()

# Plot GO BP results
jpeg("plots/11_GO_BP_utr_vs_sac_w8_ups.jpeg", width = 20, height = 20, units = "cm", res = 300,)
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 50, y = "Count",
             orderBy = "P.value", title = "Utricle vs Saccule in W8 - Mm - UPs - GO (BP)")
}
dev.off()

### GO and KEGG - Human data results
hs_utr_vs_sac <- read.csv("~/projects/vestibular_project_human/differential_expression/DEA_result_utricle_vs_saccule.csv")[, c(1,3,7)]

hs_utr_vs_sac_up_doen_gene <- hs_utr_vs_sac %>% 
  filter(log2FoldChange > 0.5 | log2FoldChange < -0.5) %>%
  filter(padj < 0.1) %>% 
  pull(gene_id) %>% 
  sort(., decreasing = T)

gene_ids_hs_utr_vs_sac <- bitr(hs_utr_vs_sac_up_doen_gene, fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)
# Enrichment
if (websiteLive) {
  enriched_hs <- enrichr(gene_ids_hs_utr_vs_sac$SYMBOL, dbs)
}

# Plot KEGG results
jpeg("plots/12_KEGG_utr_vs_sac_hs_all.jpeg", width = 20, height = 20, units = "cm", res = 300,)
if (websiteLive) {
  plotEnrich(enriched_hs[[4]], showTerms = 20, numChar = 60, y = "Count",
             orderBy = "P.value", title = "Utricle vs Saccule - Hs - KEGG")
}
dev.off()

# Plot GO BP results
jpeg("plots/13_GO_BP_utr_vs_sac_hs_all.jpeg", width = 20, height = 20, units = "cm", res = 300,)
if (websiteLive) {
  plotEnrich(enriched_hs[[3]], showTerms = 20, numChar = 60, y = "Count",
             orderBy = "P.value", title = "Utricle vs Saccule - Hs - GO (BP)")
}
dev.off()
