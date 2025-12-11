#!/usr/bin/env Rscript

rm(list = ls())
getwd()
setwd("Pantry/phenotyping/")

# load libraries
library(GenomicFeatures)
library(dplyr)
library(readr)
library(AnnotationDbi)
library(tximport)
library(DESeq2)
library(matrixStats)
library(BiocParallel)


# Build tx2gene from GTF
txdb <-makeTxDbFromGFF("input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# See what key types are available
keytypes(txdb)

# if TXNAME is listed in keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")  


# Matching transcripts and gene IDs
tx2gene_raw <- AnnotationDbi::select(txdb,
                        keys = k,
                        keytype = "TXNAME",
                        columns = "GENEID")

tx2gene <- tx2gene_raw %>%
    dplyr::select(tx = TXNAME, gene = GENEID) %>%
    distinct()

# Save it
# write_tsv(tx2gene, "input/ref/tx2gene.tsv")

### Prepare sample info
samples_saccule <- list.dirs("intermediate/saccule/expression/", full.names = F)[-1]
sample_utricle <- list.dirs("intermediate/utricle/expression/", full.names = F)[-1]

# Paired samples present in both tissues
common_samples <- intersect(samples_saccule, sample_utricle)
length(common_samples)  # should be 74 in your case

# Build file paths for kallisto abundance.tsv
files_saccule <- file.path("intermediate/saccule/expression", common_samples, "abundance.tsv")
files_utricle <- file.path("intermediate/utricle/expression", common_samples, "abundance.tsv")

# Make unique library names: sample_tissue
sample_names_saccule <- paste0(common_samples, "_saccule")
sample_names_utricle <- paste0(common_samples, "_utricle")

files <- c(files_saccule, files_utricle)
names(files) <- c(sample_names_saccule, sample_names_utricle)

# subject ID = original sample ID (the common_samples)
subject <- c(common_samples, common_samples)

tissue <- factor(
  c(rep("saccule", length(common_samples)),
    rep("utricle", length(common_samples))),
  levels = c("saccule", "utricle")
)

coldata <- data.frame(
    row.names = names(files),
    subject   = subject,
    tissue    = tissue
)
# coldata[1:6, ]


# Import kallisto quantification with tximport
txi <- tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM",  # or "no" if you prefer raw counts
    ignoreTxVersion = TRUE                    # important if IDs have version suffix
)

# Making DEseq2 object
dds <- DESeqDataSetFromTximport(
    txi     = txi,
    colData = coldata,
    design  = ~ subject + tissue
)

# Filter low-count genes (adjust thresholds if needed)
# raw counts
cts <- counts(dds)
# SD-based filtering
gene_sd    <- rowSds(cts)
sd_cutoff  <- quantile(gene_sd, 0.25)
keep_sd    <- gene_sd > sd_cutoff
# count-depth filtering
keep_counts <- rowSums(cts >= 10) >= 3
# combined filter
keep <- keep_sd & keep_counts
# apply to dds
dds <- dds[keep, ]

# Run DE analysis
register(MulticoreParam(28))
# register(SerialParam())
dea_obj_HumID_fitted <- DESeq(dds, parallel = TRUE)
saveRDS(dea_obj_HumID_fitted, "Pantry/phenotyping/intermediate/vestibular/dea_obj_HumID_fitted_afterDESeq.rds")


### END




# res <- results(dds, contrast = c("tissue", "utricle", "saccule"))
# summary(res)

# res_df <- as.data.frame(res)
# res_df$gene_id <- rownames(res_df)
# res_df <- res_df[order(res_df$padj), ]
# head(res_df)


# ### Volcano plot
# library(ggplot2)
# library(dplyr)

# # Make a plotting-friendly data frame
# volcano_df <- res_df %>%
#   mutate(
#     neglog10_padj = -log10(padj),
#     sig = case_when(
#       padj < 0.05 & log2FoldChange > 1  ~ "Up",
#       padj < 0.05 & log2FoldChange < -1 ~ "Down",
#       TRUE ~ "NS"
#     )
#   )

# # Basic volcano plot
# ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10_padj, color = sig)) +
#     geom_point(alpha = 0.7, size = 1.2) +
#     scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
#     theme_bw(base_size = 14) +
#     xlab("Log2 Fold Change") +
#     ylab("-Log10 adjusted p-value") +
#     ggtitle("Volcano plot") +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank()
#     )

# # Optional: Label the top genes
# library(ggrepel)

# top_labels <- volcano_df %>%
#     filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
#     slice_max(neglog10_padj, n = 20)

# ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10_padj, color = sig)) +
#     geom_point(alpha = 0.6, size = 1.1) +
#     geom_text_repel(data = top_labels,
#                     aes(label = gene_id),
#                     size = 3,
#                     max.overlaps = 50) +
#     scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
#     theme_bw(base_size = 14) +
#     xlab("Log2 Fold Change") +
#     ylab("-Log10 adjusted p-value") +
#     ggtitle("Volcano plot with top genes labeled") +
#     theme(plot.title = element_text(hjust = 0.5))




# # Optional: shrink LFCs (good for ranking/plots):
# # BiocManager::install("apeglm")
# library(apeglm)

# res_shrink <- lfcShrink(dds,
#                         coef = "tissue_utricle_vs_saccule",
#                         type = "apeglm")

# res_shrink_df <- as.data.frame(res_shrink)
# res_shrink_df$gene_id <- rownames(res_shrink_df)
# res_shrink_df <- res_shrink_df[order(res_shrink_df$padj), ]

# head(res_shrink_df)


# write.csv(res_shrink_df,
#           file = "DESeq2_saccule_vs_utricle_paired.csv",
#           row.names = FALSE)

# # Quick sanity-check plots
# # PCA
# vsd <- vst(dds, blind = FALSE)
# plotPCA(vsd, intgroup = "tissue")

# # MA plot
# plotMA(res_shrink, ylim = c(-5, 5))
