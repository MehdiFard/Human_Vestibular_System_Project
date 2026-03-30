#!/usr/bin/env Rscript
### DGE analysis for  utricle vs saccule

rm(list = ls())

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

# Identify paired subjects
paired_subjects <- intersect(sac_name, utr_name)

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


# Make sure sample_info matches count matrix columns
# sample_info$human_id <- factor(sample_info$human_id)
sample_info$tissue   <- factor(sample_info$tissue)
sample_info$batch    <- factor(sample_info$batch)

# Optional: check order
stopifnot(all(colnames(count_mtx_paired) == rownames(sample_info)))

# Filter low-count genes before building DESeq2 object
# Keep genes with count >= 10 in at least 10 samples
keep <- rowSums(count_mtx_paired >= 10) >= 10
count_mtx_filtered <- count_mtx_paired[keep, ]

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mtx_filtered),
  colData   = sample_info,
  design    = ~ batch + tissue
)

### Run DESeq2 with parallel processing
library(BiocParallel)
register(MulticoreParam(16)) # or
# register(SerialParam())  # important
dea_obj_batch_fitted <- DESeq(dds, parallel = TRUE)
saveRDS(dea_obj_batch_fitted, "figures/data_tables/DEA_obj_DESeq_normalized.rds")
dea_obj_batch_fitted = readRDS("figures/data_tables/DEA_obj_DESeq_normalized.rds")

# Extract results: utricle vs saccule
resultsNames(dea_obj_batch_fitted)
res <- results(dea_obj_batch_fitted, contrast = c("tissue", "utricle", "saccule"))

# Save results to a data frame and write to TSV
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]
write_tsv(res_df[, c(7, 1:6)], "figures/data_tables/DEA_result_vestibular.tsv")

# Filter significant DEGs (padj < 0.05 and |log2FC| >= 1)
results_utr_vs_sac_signif <- res_df %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "up", "down")) %>%
  # mutate(tissue = rep("saccule",dim(.)[1])) %>%
  mutate(comparision = rep("utr_vs_sac",dim(.)[1])) %>%
  dplyr::select(gene_id, comparision, log2FoldChange, padj, up_down)


### Visualizing
# Volcano plots (https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html)

DegResUtrVsSac <- data.frame(gene_symbol = res_df$gene_id, log2FoldChange = res_df$log2FoldChange, pvalue_adj = res_df$padj)
# add a column of NAs
DegResUtrVsSac$Annotation <- "NO"
# if log2Foldchange > 2 and pvalue_adj < 0.05, set as "UP" 
DegResUtrVsSac$Annotation[DegResUtrVsSac$log2FoldChange > 2 & DegResUtrVsSac$pvalue_adj < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue_adj < 0.05, set as "DOWN"
DegResUtrVsSac$Annotation[DegResUtrVsSac$log2FoldChange < -2 & DegResUtrVsSac$pvalue_adj < 0.05] <- "DOWN"

# plot and color the points with "Annotation"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

# Filter out rows with NA in pvalue_adj
DegResUtrVsSac_clean <- DegResUtrVsSac[!is.na(DegResUtrVsSac$pvalue_adj), ]

# Volcano plot
p <- ggplot(data=DegResUtrVsSac_clean) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue_adj)), color = "deepskyblue4") + # , col=Annotation
    theme_minimal() +
    # labs(title = "Utricle vs Saccule") +
    xlab(expression(log[2]("FC"))) +
    ylab(expression(-log[10]*"pvalue")) +
    scale_x_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) +
    scale_y_continuous(limits = c(0, 4)) +
    geom_vline(xintercept = c(-2, 2), color = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dashed", linewidth = 0.4) +
    theme_classic() +
    # scale_colour_manual(values = mycolors) +
    theme(legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
    ); p

# Save the plot
ggsave(filename = "VolcanoPlot_utricle_vs_saccule_batch_corrected.pdf",
        plot = p,
        units = "mm",
        width = 80,
        height = 80,
        path = "figures/",
        dpi = 300)


### END

