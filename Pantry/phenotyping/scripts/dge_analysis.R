#!/usr/bin/env Rscript

rm(list = ls())
getwd()

# Load packages
suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(data.table)
    library(stringr)
    library(GenomicFeatures)  # for GTF parsing and exon lengths
    library(GenomicRanges)
    library(edgeR)
    library(DESeq2)
    library(rtracklayer)
    library(jsonlite)
    library(matrixStats)
    library(sva)
    library(tidyr)
    library(readr)
})


# Making a vector contains protein coding and lncRNAs gene IDs to filter expression data
# Reading GTF file
gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
    dplyr::filter(type == "gene") %>% 
    dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
    dplyr::pull("gene_id")

# Read Saccule 
saccule_expr_count <- fread("Pantry/phenotyping/output/saccule/unnorm/expression_count.bed")[, -1:-3]
saccule_expr_count[1:4, 1:4]
sac_name <- names(saccule_expr_count) [-1]
names(saccule_expr_count)[-1] <- paste0(names(saccule_expr_count)[-1], "_sac") # Add tissue suffix


# Read Utricle
utricle_expr_count <- fread("Pantry/phenotyping/output/utricle/unnorm/expression_count.bed")[, -1:-3]
utr_name <- names(utricle_expr_count)[-1]
names(utricle_expr_count)[-1] <- paste0(names(utricle_expr_count)[-1], "_utr") # Add tissue suffix


# Merge saccule and utricle and keeping protien and lncRNA genes
sac_utr_expr_count <- saccule_expr_count %>%
    inner_join(utricle_expr_count, by = "phenotype_id") %>%
    dplyr::rename(gene_id = phenotype_id) %>%
    filter(gene_id %in% genes_protien_lncRNA)

# Read metadata
meta <- fread("src/metadata_rna.tsv")[, c("human_id", "tissue", "age", "sex", "batch")]
meta_sac_utr_unique <- unique(meta, by = "human_id")
meta_sac_utr_unique$human_id <- as.character(meta_sac_utr_unique$human_id)


# Counts matrix (genes x samples)
counts = sac_utr_expr_count
counts[1:3,1:3]
# stopifnot("gene_id" %in% names(counts))
# samples <- setdiff(names(counts), "gene_id")
# count_mtx <- as.matrix(counts[, ..samples]) # To matrix
# rownames(count_mtx) <- counts$gene_id
# count_mtx[1:5,1:5]


# Identify paired subjects
paired_subjects <- intersect(sac_name, utr_name)
length(paired_subjects)

# Keep only columns for paired subjects
counts = sac_utr_expr_count
samples <- setdiff(names(counts), "gene_id")

keep_cols <- samples[
    str_remove(samples, "_(sac|utr)$") %in% paired_subjects
]
length(keep_cols)   # should be 74*2 = 148 samples

# Subset counts matrix
count_mtx_paired <- as.matrix(counts[, ..keep_cols])
rownames(count_mtx_paired) <- counts$gene_id
count_mtx_paired <- round(count_mtx_paired)
count_mtx_paired[1:3,1:3]

# Build accurate sample_info (colData)
sample_info_extratc <- data.frame(
    sample = keep_cols,
    human_id = str_remove(keep_cols, "_(sac|utr)$"),
    tissue = ifelse(str_detect(keep_cols, "_sac$"), "saccule", "utricle"),
    stringsAsFactors = FALSE
)

sample_info <- sample_info_extratc %>%
    left_join(meta_sac_utr_unique[, -c("tissue")], by = "human_id") %>%
    mutate(
        tissue = factor(tissue, levels = c("saccule", "utricle")),
        sex = factor(sex),
        batch = factor(batch)
    )
rownames(sample_info) <- sample_info$sample

# Check alignment:
all(rownames(sample_info) == colnames(count_mtx_paired))


# Standard deviation for each gene across samples
gene_sd <- matrixStats::rowSds(count_mtx_paired)
sd_cutoff <- quantile(gene_sd, 0.25) # 25th percentile cutoff
# summary(gene_sd)
# hist(gene_sd)
sd_cutoff # example: maybe ~ 5–20 depending on counts
keep_sd <- gene_sd >= sd_cutoff
count_mtx_filtered <- count_mtx_paired[keep_sd, ]
dim(count_mtx_filtered)
count_mtx_filtered[1:5, 1:5]

# Build DESeq2 object with paired design
design = ~ human_id + tissue

sample_info$human_id <- as.factor(sample_info$human_id)
sample_info$tissue <- as.factor(sample_info$tissue)
sample_info$batch <- as.factor(sample_info$batch)
head(sample_info)

dds <- DESeqDataSetFromMatrix(
    countData = count_mtx_filtered,
    colData = sample_info,
    design = ~ batch + tissue
)
# batch + 
colData(dds)
counts(dds)[1:3]


# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
dds


library(BiocParallel)
register(MulticoreParam(32)) # or
# register(SerialParam())  # important
dea_obj_batch_fitted <- DESeq(dds, parallel = TRUE)
saveRDS(dea_obj_batch_fitted, "Pantry/phenotyping/intermediate/vestibular/dea_obj_batch_fitted_afterDESeq.rds")
dea_obj_batch_fitted = readRDS("Pantry/phenotyping/intermediate/vestibular/dea_obj_batch_fitted_afterDESeq.rds")

# Extract results: utricle vs saccule
resultsNames(dea_obj_batch_fitted)
res <- results(dea_obj_batch_fitted, contrast = c("tissue", "utricle", "saccule"))
View(res)
readr::write_tsv(res, "Pantry/phenotyping/intermediate/vestibular/deseq2_result.tsv")

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]
dim(res_df)
head(res_df)

results_utr_vs_sac_signif <- res_df %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "up", "down")) %>%
  # mutate(tissue = rep("saccule",dim(.)[1])) %>%
  mutate(comparision = rep("utr_vs_sac",dim(.)[1])) %>%
  dplyr::select(gene_id, comparision, log2FoldChange, padj, up_down)
dim(results_utr_vs_sac_signif)

# # Optional shrink:
# library(apeglm)
# res_shrink <- lfcShrink(
#   dds,
#   coef = "tissue_utricle_vs_saccule",
#   type = "apeglm"
# )


# Volcano plot (clean version)
library(ggplot2)
library(dplyr)

volcano_df <- res_df %>%
  mutate(
    neglog10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in utricle",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in saccule",
      TRUE ~ "Not significant"
    )
  )

ggplot(volcano_df, aes(log2FoldChange, neglog10_padj, color = sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c(
    "Up in utricle" = "red",
    "Up in saccule" = "blue",
    "Not significant" = "grey"
  )) +
  theme_bw() +
  xlab("log2 Fold Change (utricle vs saccule)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot: Paired DE between utricle and saccule")



### END

head(res_df)

# Visualizing - Volcano plots (https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html)
library(ggrepel)
library(ggplot2)

DegResUtrVsSac <- data.frame(gene_symbol = res_df$gene_id, log2FoldChange = res_df$log2FoldChange, pvalue_adj = res_df$padj)
# add a column of NAs
DegResUtrVsSac$Annotation <- "NO"
# if log2Foldchange > 1 and pvalue_adj < 0.05, set as "UP" 
DegResUtrVsSac$Annotation[DegResUtrVsSac$log2FoldChange > 1 & DegResUtrVsSac$pvalue_adj < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue_adj < 0.05, set as "DOWN"
DegResUtrVsSac$Annotation[DegResUtrVsSac$log2FoldChange < -1 & DegResUtrVsSac$pvalue_adj < 0.05] <- "DOWN"

# plot and color the points with "Annotation"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

DegResUtrVsSac_clean <- DegResUtrVsSac[!is.na(DegResUtrVsSac$pvalue_adj), ]

p <- ggplot(data=DegResUtrVsSac_clean) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue_adj), col=Annotation)) +
    theme_minimal() +
    labs(title = "Utricle vs Saccule - Human") +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") +
    theme_classic() +
    scale_colour_manual(values = mycolors); p

# vol_plots = grid.arrange(p3_saccule, p3_utricle, ncol = 2); vol_plots
ggsave(filename = "VolcanoPlot_utricle_vs_saccule_batch_fitted.jpeg", plot = p,
        path = "figures/", dpi = 300)
