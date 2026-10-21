### Gene enrichment analysis for the vestibular high expressed genes across all GTEx tissues

rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
library(ggplot2)
library(stringr)


### Reading the genes
# Sig genes
genes <- fread("analyses/08_expression_gtex/results/vestibular_high_expressed_genes_accross_all_gtex.tsv", header = F)
genes_sig <- unique(genes$V1)

# Background genes
genes_background <- fread("analyses/08_expression_gtex/results/vestibular_and_gtex_genes_background.tsv", header = F)
genes_back <- unique(genes_background$V1)


# Gene id conversion
keytypes(org.Hs.eg.db)

genes_symbol <- bitr(
  genes$V1,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)


# GO Biological Process
go_bp  <- enrichGO(
  gene          = genes_sig,
  universe      = genes_back,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

BP <- as.data.frame(go_bp)
View(BP)

# Save output
fwrite(BP, "analyses/08_expression_gtex/results/enrichment_BP.tsv", row.names = F, quote = F, sep = "\t")

# Plotting
dotplot(go_bp, showCategory = 23)
barplot(go_bp, showCategory = 23)
heatplot(go_bp)

# By ggplot2
bp_plot <- BP %>%
  as.data.frame() %>%
  slice_min(order_by = p.adjust, n = 23) %>%
  mutate(
    # Description = str_wrap(Description, width = 40),
    Description = factor(Description, levels = rev(Description)),
    GeneRatio_num = Count / as.numeric(sub(".*/", "", GeneRatio))
  )
View(bp_plot)


ggplot(bp_plot, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_viridis_c(name = "-log10(adj. p)") +
  # scale_size_continuous(name = "Gene count") +
  # theme_classic() +
  theme_bw() +
  labs(
    x = "Gene ratio",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    # plot.margin = margin(5.5, 5.5, 5.5, 40),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    # legend.key.size = unit(0.4, "cm")
  ) +
  scale_size_continuous(
    name = "Gene count",
    breaks = c(5, 10, 15, 20, 25, 30),
    guide = guide_legend(keyheight = unit(0.4, "cm"),
                          keywidth  = unit(0.4, "cm"))
)

# Save
ggsave(
  "analyses/08_expression_gtex/figures/enrichment_GO_BP.pdf",
  width = 5,
  height = 5
)


# GO Molecular Function
go_mf  <- enrichGO(
  gene          = genes_sig,
  universe      = genes_back,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

MF <- as.data.frame(go_mf)
View(MF)


# GO Cellular Component
go_cc  <- enrichGO(
  gene          = genes_sig,
  universe      = genes_back,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

CC <- as.data.frame(go_cc)
View(CC)


### KEGG
search_kegg_organism('hsa', by='kegg_code')

# Convert gene ids
# Sig genes
genes_symbol_sig <- bitr(
  genes_sig,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
genes_symbol_sig_vector <- genes_symbol_sig$ENTREZID

# Background genes
genes_symbol_bg <- bitr(
  genes_back,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
genes_symbol_bg_vector <- genes_symbol_bg$ENTREZID

# Run kegg enrichment
ekegg <- enrichKEGG(
  gene          = genes_symbol_sig_vector,
  universe      = genes_symbol_bg_vector,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
)

kegg <- as.data.frame(ekegg)




### From cytoscape
gsea_cyto <- fread("analyses/08_expression_gtex/results/enrichment_stringdb_cytoscape.csv")
fwrite(gsea_cyto, "analyses/08_expression_gtex/results/enrichment_cytoscape.tsv", row.names = F, quote = F, sep = "\t")

# Filter
gsea_plot <- as.data.frame(gsea_cyto) %>%
  dplyr::rename(
    bg_genes = `background genes`,
    gene_count = genes,
    category_name = category,
    term = description,
    FDR = `FDR value`,
    pval = `p-value`
  ) %>%
  dplyr::filter(
    !category_name %in% c(
      "GO Biological Process",
      "COMPARTMENTS",
      "GO Cellular Component",
      "GO Molecular Function",
      "SMART Domains",
      "WikiPathways",
      "STRING Clusters"
    )
  ) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  # dplyr::slice_head(n = 15) %>%
  dplyr::mutate(
    term = factor(term, levels = rev(unique(term))),
    category_name = sub(" .*", "", category_name),
    term_category = paste0(term, " (", category_name, ")"),
    gene_ratio = gene_count / bg_genes
  )


# plot
ggplot(gsea_plot, aes(x = gene_ratio, y = term_category)) +
  geom_point(aes(size = gene_count, color = -log10(FDR))) +
  scale_color_viridis_c(
    name = "-log10(adj. p)",
    # option = "magma",
    guide = guide_colorbar(order = 2)
  ) +
  scale_size_continuous(
    name = "Gene count",
    breaks = c(5, 10, 15, 20, 25, 30),
    guide = guide_legend(
      order = 1,
      keyheight = unit(0.3, "cm"),
      keywidth  = unit(0.3, "cm")
    )
  ) +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

# Save
ggsave(
  "analyses/08_expression_gtex/figures/enrichment_Uniprot_and_other.pdf",
  width = 5,
  height = 5
)

### END