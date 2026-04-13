### Here I plot the final figures should go into the paper
##--> PC plots for expression

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
})


## Load expression data
# Saccule
saccule_expr_count <- fread("Pantry/phenotyping/output/saccule/unnorm/expression_count.bed")[, -1:-3]
names(saccule_expr_count)[-1] <- paste0(names(saccule_expr_count)[-1], "_sac") # Add tissue suffix

# Utricle
utricle_expr_count <- fread("Pantry/phenotyping/output/utricle/unnorm/expression_count.bed")[, -1:-3]
names(utricle_expr_count)[-1] <- paste0(names(utricle_expr_count)[-1], "_utr") # Add tissue suffix

# Merge saccule and utricle
sac_utr_expr_count <- merge(saccule_expr_count, utricle_expr_count, by="phenotype_id")


# Read metadata
meta <- fread("src/metadata_rna.tsv")[, c("human_id", "tissue", "age", "sex", "batch")]
meta_unique <- unique(meta, by = "human_id")
meta_unique$human_id <- as.character(meta_unique$human_id)
head(meta_unique)


# ----------------------------
# User inputs
# ----------------------------
out_prefix       <- "rna_tissue"
gtf_file <- "Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf"

# ----------------------------
# 1) Load counts
# ----------------------------
counts <- sac_utr_expr_count
names(counts)[1] <- "gene_id"
stopifnot("gene_id" %in% names(counts))
samples <- setdiff(names(counts), "gene_id")


# ----------------------------
# 2) GENE LENGTHS FROM GTF (non-overlapping exon union per gene)
# ----------------------------
# Build TxDb from GTF
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
ex_by_gene <- exonsBy(txdb, by = "gene") # Extract exons grouped by gene (GRangesList)
ex_reduced <- reduce(ex_by_gene) # Reduce overlaps within each gene (union of exons)
lens <- sum(width(ex_reduced))  # Sum widths to get per-gene length (bp), named integer vector (names are gene IDs)

# Convert to data.table
gene_lengths <- data.table(gene_id = names(lens), length = as.integer(lens))
gene_lengths <- gene_lengths[!is.na(length) & length > 0]


# ----------------------------
# 3) Merge and make matrix
# ----------------------------
expr <- merge(counts, gene_lengths, by = "gene_id")
expr <- expr[!is.na(length) & length > 0]
mat <- as.matrix(expr[, ..samples])
rownames(mat) <- expr$gene_id

# ----------------------------
# 4) TPM
# ----------------------------
# RPK = counts / (length in kb)
rpk <- sweep(mat, 1, expr$length / 1000, "/")
# scaling factor per sample = sum RPK / 1e6
sf  <- colSums(rpk, na.rm = TRUE) / 1e6
tpm <- sweep(rpk, 2, sf, "/")
tpm[1:5, 1:5]

# ----------------------------
# 5) Filter lowly expressed genes
# ----------------------------
min_tpm <- 1
min_samples <- 0.2  # keep genes with TPM>1 in ≥20% of samples
keep <- rowSums(tpm > min_tpm) >= (min_samples * ncol(tpm))
tpm_filt <- tpm[keep, , drop = FALSE]
cat("Genes kept:", sum(keep), "/", nrow(tpm), "\n")

# ----------------------------
# 6) Log2 transform
# ----------------------------
tpm_log <- log2(tpm_filt + 1)

# ----------------------------
# 7) PCA (samples as rows)
# ----------------------------
pca <- prcomp(t(tpm_log), scale. = TRUE)
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)

# Adding tissue info
pc_df$Tissue <- ifelse(grepl("_sac$", pc_df$Sample), "Saccule",
                ifelse(grepl("_utr$", pc_df$Sample), "Utricle", NA))


## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
pc_df[, human_id := as.integer(sub("_(sac|utr)$", "", Sample))] # Extract numeric ID from Sample (remove "_sac" or "_utr")
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(Sample, human_id, Tissue, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)

# Explained variance for axis labels
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
# Quick check
table(pc_df$Tissue, useNA = "ifany")

# Plot by tissue
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TPM expression",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_count_tissue.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

# Plot by batch
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TPM expression",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_count_batch.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

# ----------------------------
# 8) Save outputs
# ----------------------------
fwrite(as.data.table(tpm, keep.rownames = "gene_id"),       paste0(out_prefix, "_TPM_all.tsv"),      sep = "\t")
fwrite(as.data.table(tpm_filt, keep.rownames = "gene_id"),  paste0(out_prefix, "_TPM_filtered.tsv"), sep = "\t")
fwrite(as.data.table(pc_df),                                 paste0(out_prefix, "_PCA_scores.tsv"),   sep = "\t")
fwrite(gene_lengths,                                         paste0(out_prefix, "_gene_lengths.tsv"), sep = "\t")


###------------------------------------ Normalization with TMM + logCPM and PCA plot
# Counts (genes x samples)
counts = sac_utr_expr_count
names(counts)[1] <- "gene_id"
stopifnot("gene_id" %in% names(counts))
samples <- setdiff(names(counts), "gene_id")
mat <- as.matrix(counts[, ..samples]) # To matrix
rownames(mat) <- counts$gene_id

# Create DGEList and apply TMM normalization
y <- DGEList(counts = mat)
y <- calcNormFactors(y, method = "TMM") # It computes normalization factors and stores them in y$samples$norm.factors, which are then used to adjust library sizes

# Convert to log2-CPM (counts per million)
logCPM <- cpm(y, log = TRUE, prior.count = 1)

## Filtering steps before PCA
# Remove genes with any NA or non-finite values
logCPM <- logCPM[apply(logCPM, 1, function(x) all(is.finite(x))), ]

# Remove genes with zero variance across samples
gene_var <- apply(logCPM, 1, var)
logCPM <- logCPM[gene_var > 0, ]

# Optional: remove genes with very low expression everywhere. (If you came directly from CPM/logCPM, low-expression genes can add noise)
mean_expr <- rowMeans(logCPM)
logCPM <- logCPM[mean_expr > 0.5, ]   # adjust threshold as needed

# (Optional but good) Remove outlier samples with zero variance or too many NAs
sample_var <- apply(logCPM, 2, var)
logCPM <- logCPM[, sample_var > 0, drop = FALSE]

# Check dimensions
cat("Genes remaining:", nrow(logCPM), "\nSamples remaining:", ncol(logCPM), "\n")


# 4. Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)


# Adding tissue info
pc_df$Tissue <- ifelse(grepl("_sac$", pc_df$Sample), "Saccule",
                ifelse(grepl("_utr$", pc_df$Sample), "Utricle", NA))

## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
pc_df[, human_id := as.character(sub("_(sac|utr)$", "", Sample))] # Extract numeric ID from Sample (remove "_sac" or "_utr")
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(Sample, human_id, Tissue, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)

# Explained variance for axis labels
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
# Quick check
table(pc_df$Tissue, useNA = "ifany")

# Plot by tissue
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM expression",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_tissue2.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

# Plot by batch
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM expression",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_batch2.png",
  width = 9,
  height = 7.5,
  dpi = 300
)


###-------------------------------------------- Normalizing Vestibular RNA expression data with TMM + logCPM and PCA

# vestibular
vestibular_expr_count <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression_count.bed")[, -1:-3]
vestibular_expr_count[1:5, 1:5]
counts <- as.matrix(vestibular_expr_count[, -1])
rownames(counts) <- vestibular_expr_count$phenotype_id
counts[1:5, 1:5]
dim(counts)


# Make sure counts are non-negative integers
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, , drop = FALSE]  # drop all-zero genes

# edgeR object (only for filtering; can still switch to DESeq2 later)
dge0 <- DGEList(counts = counts)

# Rule: CPM > 1 in at least 20% of samples (tweak as needed)
n <- ncol(counts)
min_samples <- max(10, ceiling(0.5 * n))  # at least 2 samples
keep <- rowSums(cpm(dge0) > 1) >= min_samples

counts_filt <- counts[keep, ]
cat("Kept", nrow(counts_filt), "genes out of", nrow(counts), "\n")

# TMM normalization
dge <- DGEList(counts = counts_filt)
dge <- calcNormFactors(dge, method = "TMM")  # composition bias correction

# Log2(TMM-normalized CPM); add a small prior.count to stabilize low counts
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# Filtering based on SD
library(matrixStats)

row_sd <- rowSds(logCPM)
thr <- quantile(row_sd, 0.50)
keep_sd <- row_sd >= thr

logCPM_top <- logCPM[keep_sd, , drop = FALSE]
cat("Kept", nrow(logCPM_top), "genes out of", nrow(logCPM), "by SD ≥ median\n")

# PCA (samples as rows, genes as columns)
pca <- prcomp(t(logCPM), scale. = TRUE)
pca$x[1:5, 1:5]
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]


## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)

# Explained variance for axis labels
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)


indx <- which(pc_df_merged$human_id == "31947386") # Find index of outlier sample to remove


# Plot by tissue
ggplot(pc_df_merged[-indx], aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM expression - vestibular",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_vestibular.png",
  width = 9,
  height = 7.5,
  dpi = 300
)






###-------------------------------> Vestibular - TPM

gtf_file <- "Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf"

# ----------------------------
# 1) Load counts
# ----------------------------
counts_df <- vestibular_expr_count
counts_df[1:5, 1:5]
names(counts_df)[1] <- "gene_id"
stopifnot("gene_id" %in% names(counts_df))
samples <- setdiff(names(counts_df), "gene_id")


# ----------------------------
# 2) GENE LENGTHS FROM GTF (non-overlapping exon union per gene)
# ----------------------------
# Build TxDb from GTF
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
ex_by_gene <- exonsBy(txdb, by = "gene") # Extract exons grouped by gene (GRangesList)
ex_reduced <- reduce(ex_by_gene) # Reduce overlaps within each gene (union of exons)
lens <- sum(width(ex_reduced))  # Sum widths to get per-gene length (bp), named integer vector (names are gene IDs)

# Convert to data.table
gene_lengths <- data.table(gene_id = names(lens), length = as.integer(lens))
gene_lengths <- gene_lengths[!is.na(length) & length > 0]
head(gene_lengths)

# ----------------------------
# 3) Merge and make matrix
# ----------------------------
expr <- merge(counts_df, gene_lengths, by = "gene_id")
expr <- expr[!is.na(length) & length > 0]
counts_mat <- as.matrix(expr[, ..samples])
rownames(counts_mat) <- expr$gene_id
counts_mat[1:5, 1:5]
dim(counts_mat)

# ----------------------------
# 4) TPM
# ----------------------------
# RPK = counts / (length in kb)
rpk <- sweep(counts_mat, 1, expr$length / 1000, "/")
# scaling factor per sample = sum RPK / 1e6
sf  <- colSums(rpk, na.rm = TRUE) / 1e6
tpm <- sweep(rpk, 2, sf, "/")
tpm[1:5, 1:5]
dim(tpm)

# ----------------------------
# 5) Filter lowly expressed genes
# ----------------------------
min_tpm <- 1
min_samples <- 0.2  # keep genes with TPM>1 in ≥20% of samples
keep <- rowSums(tpm > min_tpm) >= (min_samples * ncol(tpm))
tpm_filt <- tpm[keep, , drop = FALSE]
cat("Genes kept:", sum(keep), "/", nrow(tpm), "\n")

# ----------------------------
# 6) Log2 transform
# ----------------------------
tpm_log <- log2(tpm_filt + 1)

# ----------------------------
# 7) PCA (samples as rows)
# ----------------------------
pca <- prcomp(t(tpm_log), scale. = TRUE)


pca$x[1:5, 1:5]
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]


## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)

# Explained variance for axis labels
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
# Quick check
table(pc_df_merged$tissue, useNA = "ifany")

# Plot by tissue
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TPM expression - vestibular",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tpm_vestibular.png",
  width = 9,
  height = 7.5,
  dpi = 300
)





###-------------------------------> Vestibular - DESeq2 VST
counts_mat = counts_mat


# ---- Packages ----
library(DESeq2)
library(ggplot2)
library(dplyr)

# ---- 0) Start from your counts matrix ----
# counts: genes x samples, rownames = Ensembl IDs, colnames = sample IDs
stopifnot(is.matrix(counts_mat))
storage.mode(counts_mat) <- "integer"                 # DESeq2 expects integer counts
counts <- counts_mat[rowSums(counts_mat) > 0, , drop=FALSE]   # drop all-zero genes

# ---- 1) Filter low-expressed genes (count-based, before normalization) ----
n <- ncol(counts)
min_samples <- max(2, ceiling(0.2 * n))           # present in ≥20% samples (tweak as needed)
keep <- rowSums(counts >= 10) >= min_samples      # e.g., at least 10 counts in enough samples
counts_filt <- counts[keep, , drop=FALSE]

cat("Kept", nrow(counts_filt), "genes out of", nrow(counts), "\n")

# ---- 2) Build colData (sample metadata) ----
# If you already have a sample_info data.frame, use it here.
# Otherwise, create a minimal one and fill in later (tissue/source/batch if you have them).
sample_info <- data.frame(
  sample = colnames(counts_filt),
  row.names = colnames(counts_filt)
  # , tissue = ..., source = ..., batch = ...
)

# ---- 3) Create DESeq2 object & estimate size factors ----
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = meta_unique,
  design    = ~ 1   # design not needed for PCA; set ~ batch + tissue for DE later
)

dds <- estimateSizeFactors(dds)   # library-size normalization
# Optional: normalized (non-log) counts if you need them
norm_counts <- counts(dds, normalized = TRUE)
norm_counts[1:5, 1:5]


# ---- 4) Variance-stabilizing transform (best for PCA) ----
vst_mat <- assay(vst(dds, blind = TRUE))  # matrix: genes x samples
vst_mat[1:5, 1:5]

# Filtering based on SD
row_sd <- rowSds(vst_mat)
thr <- quantile(row_sd, 0.50)          # 50% quantile (median)
keep_sd <- row_sd >= thr               # keep top half by SD

vst_top <- vst_mat[keep_sd, , drop = FALSE]
cat("Kept", nrow(vst_top), "genes out of", nrow(vst_mat), "by SD ≥ median\n")

# ---- 5) PCA on samples ----
pca <- prcomp(t(vst_mat), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)



pca$x[1:5, 1:5]
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]


## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)


# Plot by tissue
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of DEseq2 vst expression - vestibular",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_deseq2_vst_vestibular.png",
  width = 9,
  height = 7.5,
  dpi = 300
)







########################### Remove batch effects with ComBat - batch correction

library(sva)

# Check that the matrix and metadata align
all(colnames(vst_mat) == meta_unique$human_id)  # should be TRUE
head(meta_unique)

# Just batch — no biological covariate
batch <- meta_unique$batch

# Model matrix with only intercept (no biological groups)
mod0 <- model.matrix(~ 1, data = meta_unique)

# Apply ComBat
vst_corrected <- ComBat(dat = vst_mat, batch = batch, mod = mod0, par.prior = TRUE)
vst_corrected[1:5, 1:5]
head(meta_unique)


# ---- 5) PCA on samples ----
pca <- prcomp(t(vst_corrected), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

pca$x[1:5, 1:5]
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]

## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, age, sex, batch)] # Inspect
dim(pc_df_merged)
pc_df_merged[1:5, 1:10]
indx <- which(pc_df_merged$human_id == "31947386") # Find index of outlier sample to remove

# Plot by tissue
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = batch, label = human_id)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of DEseq2 vst expression - vestibular - batch corrected",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_deseq2_vst_vestibular_batch_corrected.png",
  width = 9,
  height = 7.5,
  dpi = 300
)





###---------------------- Adding GTEx data
## Making a vector contains protein coding genes and lncRNAs
# Reading GTF file
library(rtracklayer)
gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf") 
# Getting a list of genes that are protein coding or lncRNA 
genes_protien_lncRNA <- gtf %>%
  dplyr::filter(type == "gene") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
  dplyr::pull("gene_id")


# Load GTEx colors from JSON
library(jsonlite)
# gtex_colors <- fromJSON("src/gtex_colors.json")
gtex_colors <- fromJSON("https://raw.githubusercontent.com/broadinstitute/gtex-viz/master/boxplot/dev/colors.json")

# Convert to data frame
gtex_df <- bind_rows(lapply(names(gtex_colors), function(tissue) {
  data.frame(
    tissue_name = tissue,
    tissue_abbrv = gtex_colors[[tissue]]$tissue_abbrv,
    tissue_color_hex = paste0("#", gtex_colors[[tissue]]$tissue_color_hex),
    stringsAsFactors = FALSE
  )
}))

## Clean tissue names
# Cells - Transformed fibroblasts == Cells - Cultured fibroblasts, so replacing it to match with tissue id in GTEx8
gtex_df <- gtex_df %>%
  mutate(
    tissue_name_clean = tissue_name %>%
      str_to_lower() %>%
      str_replace_all(c(" - " = "_", " \\(" = "_", "\\)" = "", " " = "_", "cells_transformed" = "cells_cultured"))
  )
dim(gtex_df)
head(gtex_df)

# Make a named vector (for ggplot)
gtex_colors <- setNames(gtex_df$tissue_color_hex, gtex_df$tissue_abbrv)
gtex_abb <- setNames(gtex_df$tissue_name_clean, gtex_df$tissue_abbrv)

# Getting list of GTEx tissue files
files_list <- list.files("src/GTEx/tissues_raw_count", pattern = "gene_reads")
tissue_names <- str_extract(files_list, "(?<=v8_)[^\\.]+")

# Check matching with gtex_abb from JSON
table(tissue_names %in% gtex_abb)

# Reading and merging all GTEx tissue count data
l1 = length(tissue_names)
for (i in 1:l1) {
  cat("Tissue:", tissue_names[i], "- loop No. =", i, "\n")
  
  gtex_tissue_count <- fread(paste0("src/GTEx/tissues_raw_count/", files_list[i]))[, -c(1,3)]
  colnames(gtex_tissue_count)[1] <- "gene_id"
  gtex_tissue_count[1:5, 1:5]
  
  l2 = dim(gtex_tissue_count)[2]
  abb <- names(which(gtex_abb == tissue_names[i]))
  annot_gtx_tissue <- data.frame(tissue = rep(abb, l2-1), row.names = colnames(gtex_tissue_count)[-1])
  head(annot_gtx_tissue)
  
  if(i == 1) {
    gtex_tissue_count_all = gtex_tissue_count
    annot_gtx_tissue_all = annot_gtx_tissue
  } else {
    gtex_tissue_count_all <- merge(gtex_tissue_count_all, gtex_tissue_count, by = "gene_id", sort = F)
    annot_gtx_tissue_all <- rbind(annot_gtx_tissue_all, annot_gtx_tissue)
  }
}
annot_gtx_tissue_all$batch <- "gtex" # Adding batch info
gtex_tissue_count_all$gene_id <- gsub("\\..*", "", gtex_tissue_count_all$gene_id) # Remove version from Ensembl IDs
dim(gtex_tissue_count_all)


## Reading vestibualar count data
vestibular_expr_count <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression_count.bed")[, -1:-3]
colnames(vestibular_expr_count)[1] <- "gene_id"
dim(vestibular_expr_count)


## Annotation
# Read metadata for batch info
order_vec <- as.vector(colnames(vestibular_expr_count)[-1])
meta_unique <- fread("src/metadata_rna.tsv")[, c("human_id", "batch")] %>%
  as_tibble() %>%
  distinct(human_id, .keep_all = TRUE) %>%
  mutate(
    human_id = as.character(human_id),
    tissue = "VESTIB",
    human_id = factor(human_id, levels = order_vec)
  ) %>%
  filter(!is.na(human_id)) %>%
  arrange(human_id) %>%
  mutate(human_id = as.character(human_id)) %>%
  column_to_rownames("human_id")  %>%
  dplyr::select(tissue, batch)

# Matching the order
table(rownames(meta_unique) == colnames(vestibular_expr_count)[-1]) # Check matching


# meta <- fread("src/metadata_rna.tsv")[, c("human_id", "batch")]
# meta_unique <- unique(meta, by = "human_id")
# meta_unique$human_id <- as.character(meta_unique$human_id)
# meta_unique$tissue <- rep("VESTIB", nrow(meta_unique))
# # Matching the order
# meta_unique <- meta_unique[match(colnames(vestibular_expr_count)[-1], meta_unique$human_id), ]


# Final annotation for gtex and vestibular data
annot_gtex_vest <- rbind(annot_gtx_tissue_all, meta_unique)
dim(annot_gtex_vest)
tail(annot_gtex_vest)
head(annot_gtex_vest)


## Merging GTEx and vestibular count data
# Function to collapse counts by gene_id (sum counts for same gene_id)
collapse_by_gene <- function(df) {
  df %>%
    mutate(gene_id = sub("\\..*$", "", gene_id)) %>%   # remove version suffixes; drop this line if not needed
    group_by(gene_id) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Merging by gene_id
count_mtx_gtex_vest <- gtex_tissue_count_all %>%
  filter(gene_id %in% genes_protien_lncRNA) %>%  # keep only protein coding and lncRNA genes
  collapse_by_gene() %>%
  inner_join(collapse_by_gene(vestibular_expr_count), by = "gene_id") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

dim(count_mtx_gtex_vest)
class(count_mtx_gtex_vest)
count_mtx_gtex_vest[1:5, 1:5]
table(rownames(count_mtx_gtex_vest) %in% genes_protien_lncRNA) # Should be TRUE
class(count_mtx_gtex_vest[,1])


# Make sure counts are integers
# storage.mode(count_mtx_gtex_vest) <- "integer"  # for DESeq2/edgeR


## Filtering lowly expressed genes
# n <- ncol(count_mtx_gtex_vest)
# min_samples <- max(10, ceiling(0.5 * n))  # at least 2 samples
# keep <- rowSums(count_mtx_gtex_vest >= 10) >= min_samples
# counts_filt <- count_mtx_gtex_vest[keep, , drop=FALSE]
# dim(counts_filt)



# Drop all-zero genes immediately
counts <- count_mtx_gtex_vest[rowSums(count_mtx_gtex_vest) > 0, , drop = FALSE]
dim(counts)
counts[1:5, 1:5]

library(edgeR)
# counts: genes x samples (integers)
stopifnot(is.matrix(counts))

# Build DGEList
dge <- DGEList(counts = counts)

# Compute TMM normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Get normalized log2 CPM matrix
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
dim(logCPM)
logCPM[1:5, 1:5]
save(logCPM, file = "figures/data_tables/gtex_vest_logCPM_normalized.RData")

# Filter by SD or select top genes before PCA
library(matrixStats)
row_sd <- rowSds(logCPM)
keep <- row_sd >= quantile(row_sd, 0.5)
logCPM_sd <- logCPM[keep, , drop = FALSE]
dim(logCPM_sd)

# PCA
pca <- prcomp(t(logCPM_sd), scale. = TRUE)
save(pca, file = "figures/data_tables/gtex_vest_logCPM_normalized_PCs.RData")
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)


pca$x[1:5, 1:5]
pc_df <- data.frame(human_id = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]
dim(pc_df)

dim(annot_gtex_vest)
annot_gtex_vest[1:5, ]
annot_gtex_vest[17480:17492, ]

table(rownames(annot_gtex_vest) %in% pc_df$human_id) # Should be TRUE


## Adding meta data, batch and tissue
pc_df_annot <- merge(
  pc_df,
  cbind(human_id = rownames(annot_gtex_vest), annot_gtex_vest),
  by = "human_id",
  sort = FALSE
)
dim(pc_df_annot)
pc_df_annot[1:5, 1:5]
pc_df_annot[17490:17492, c(1,16400:16407)]
pc_df_annot[1:5, c("human_id", "tissue", "batch")] 



# Add VESTIB color (choose any hex code you like)
vestib_color <- "#000000"  # black — common choice for new tissue
gtex_palette <- setNames(gtex_df$tissue_color_hex, gtex_df$tissue_abbrv)
gtex_palette <- c(gtex_palette, VESTIB = vestib_color)

# Adding tissue colors
pc_df_annot <- pc_df_annot %>%
  mutate(tissue_color_hex = gtex_palette[tissue])

# Check
unique(pc_df_annot$tissue_color_hex[pc_df_annot$tissue == "VESTIB"]) # should show "#000000"

# Inspect
pc_df_annot[1:5, c("human_id", "tissue", "batch", "tissue_color_hex")] 

# Plot by tissue
ggplot(pc_df_annot, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = gtex_palette, name = "Tissue") +
  labs(
    title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_gtex_vestibular_tissue.png",
  width = 18,
  height = 15,
  dpi = 300
)

# Plot by tissue
ggplot(pc_df_annot, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_gtex_vestibular_batch.png",
  width = 18,
  height = 15,
  dpi = 300
)



###########################################Batch correcting with ComBat-Seq - GTEx + Vestibular

# ---- 1) Filter low-expressed genes (count-based, before normalization) ----
n <- ncol(count_mtx_gtex_vest)
min_samples <- max(2, ceiling(0.2 * n))           # present in ≥20% samples (tweak as needed)
keep <- rowSums(counts >= 10) >= min_samples      # e.g., at least 10 counts in enough samples
counts_filt <- counts[keep, , drop=FALSE]
dim(counts_filt)
cat("Kept", nrow(counts_filt), "genes out of", nrow(counts), "\n")

# ---- 3) Create DESeq2 object & estimate size factors ----
is.matrix(counts_filt)
storage.mode(counts_filt) <- "integer"

library("DESeq2")
table(rownames(annot_gtex_vest) %in% colnames(counts_filt)) # Should be TRUE
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = annot_gtex_vest,
  design    = ~ 1   # design not needed for PCA; set ~ batch + tissue for DE later
)

dds <- estimateSizeFactors(dds)   # library-size normalization
# Optional: normalized (non-log) counts if you need them
norm_counts <- counts(dds, normalized = TRUE)
norm_counts[1:5, 1:5]


# ---- 4) Variance-stabilizing transform (best for PCA) ----
vst_mat <- assay(vst(dds, blind = TRUE))  # matrix: genes x samples
vst_mat[1:5, 1:5]

# Filtering based on SD
row_sd <- rowSds(vst_mat)
thr <- quantile(row_sd, 0.50)          # 50% quantile (median)
keep_sd <- row_sd >= thr               # keep top half by SD

vst_top <- vst_mat[keep_sd, , drop = FALSE]
cat("Kept", nrow(vst_top), "genes out of", nrow(vst_mat), "by SD ≥ median\n")


library(sva)

dim(annot_gtex_vest)
annot_gtex_vest[1:5, ]
annot_gtex_vest[17480:17492, ]


# Check that the matrix and metadata align
all(colnames(vst_mat) == rownames(annot_gtex_vest))  # should be TRUE
head(annot_gtex_vest)

# Ensure factors
annot_gtex_vest$batch  <- as.factor(annot_gtex_vest$batch)
annot_gtex_vest$tissue <- as.factor(annot_gtex_vest$tissue)

# Assign variables
batch <- annot_gtex_vest$batch
mod   <- model.matrix(~ 1, data = annot_gtex_vest)

# Run ComBat
vst_corrected <- ComBat(dat = vst_mat, batch = batch, mod = mod, par.prior = TRUE, ref.batch = "gtex")
vst_corrected[1:5, 1:5]
dim(vst_corrected)



# PCA
pca <- prcomp(t(vst_corrected), scale. = TRUE)
save(pca, file = "figures/data_tables/gtex_vest_DESeq2_vst_normalized_PCs.RData")
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)


pca$x[1:5, 1:5]
pc_df <- data.frame(human_id = rownames(pca$x), pca$x, check.names = FALSE)
pc_df[1:5, 1:5]
dim(pc_df)

dim(annot_gtex_vest)
annot_gtex_vest[1:5, ]
annot_gtex_vest[17480:17492, ]

table(rownames(annot_gtex_vest) %in% pc_df$human_id) # Should be TRUE


## Adding meta data, batch and tissue
pc_df_annot <- merge(
  pc_df,
  cbind(human_id = rownames(annot_gtex_vest), annot_gtex_vest),
  by = "human_id",
  sort = FALSE
)
dim(pc_df_annot)
pc_df_annot[1:5, 1:5]
pc_df_annot[17490:17492, c(1,16400:16407)]
pc_df_annot[1:5, c("human_id", "tissue", "batch")] 


# Add VESTIB color (choose any hex code you like)
vestib_color <- "#000000"  # black — common choice for new tissue
gtex_palette <- setNames(gtex_df$tissue_color_hex, gtex_df$tissue_abbrv)
gtex_palette <- c(gtex_palette, VESTIB = vestib_color)

# Adding tissue colors
pc_df_annot <- pc_df_annot %>%
  mutate(tissue_color_hex = gtex_palette[tissue])

# Check
unique(pc_df_annot$tissue_color_hex[pc_df_annot$tissue == "VESTIB"]) # should show "#000000"

# Inspect
pc_df_annot[1:5, c("human_id", "tissue", "batch", "tissue_color_hex")] 

# Plot by tissue
ggplot(pc_df_annot, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = gtex_palette, name = "Tissue") +
  labs(
    title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_gtex_vestibular_tissue_batch_corrected.png",
  width = 18,
  height = 15,
  dpi = 300
)

# Plot by tissue
ggplot(pc_df_annot, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_gtex_vestibular_batch.png",
  width = 18,
  height = 15,
  dpi = 300
)

###-----------------------------
getwd()
df = readRDS("Pantry/phenotyping/intermediate/vestibular/dea_obj_batch_fitted_afterDESeq.rds")
dim(df)
head(df)


library(dplyr)
library(ggplot2)
library(readr)

twas <- read_tsv("Pantry/pheast/twas_assoc/vestibular/twas_hits_annot.tsv")
head(twas)
View(twas)
names(twas)

# Basic derived columns
twas <- twas %>%
  mutate(
    neglog10_p = -log10(TWAS.P),
    mid_pos = (P0 + P1) / 2,              # midpoint for plotting
    is_sig = TWAS.P < 5e-8,               # example genome-wide gene-level threshold
    high_coloc = COLOC.PP4 >= 0.8
  )

twas$neglog10_p <- -log10(twas$TWAS.P)

library(ggrepel)  # install if needed

ggplot(twas, aes(x = TWAS.Z, y = neglog10_p)) +
  geom_point(aes(color = high_coloc), alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_color_manual(values = c("grey60", "red")) +
  labs(
    x = "TWAS Z-score",
    y = "-log10(TWAS p-value)",
    color = "PP4 ≥ 0.8?",
    title = "TWAS volcano plot"
  ) +
  theme_bw()
  

install.packages("ggrepel")
library(ggrepel)  # install if needed
  top_hits <- twas %>%
  filter(TWAS.P < 5e-8 | COLOC.PP4 >= 0.8)


library(ggplot2)

## Define thresholds
p_twas_thresh  <- 5e-8
coloc_thresh   <- 0.8

## Flag high-coloc genes (TRUE/FALSE)
twas$high_coloc <- twas$COLOC.PP4 >= coloc_thresh

## Select top hits to label (e.g. significant & high colocalization)
top_hits <- subset(twas, TWAS.P <= p_twas_thresh & high_coloc)

## Plot
ggplot(twas, aes(x = TWAS.Z, y = -log10(TWAS.P))) +
  geom_point(aes(color = high_coloc), alpha = 0.7) +
  geom_hline(yintercept = -log10(p_twas_thresh), linetype = "dashed") +
  geom_text(
    data = top_hits,
    aes(label = GENE.NAME),
    size = 3,
    vjust = -0.5
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey60", "TRUE" = "red"),
    name   = paste0("PP4 ≥ ", coloc_thresh)
  ) +
  theme_bw()



###------------------------------------
library(data.table)
cr <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/cross_modality.cis_qtl.txt.gz")
dim(cr)
head(cr)
View(cr[1:5, ])
table(cr$qval < 0.05)

cr_ind <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/cross_modality.cis_independent_qtl.txt.gz")
dim(cr_ind)
head(cr_ind)
View(cr_ind[1:20, ])
table(cr_ind$pval_beta < 0.05)
names(cr_ind)
fwrite(cr_ind, "cr_ind.tsv", row.names = F, col.names = T, sep = "\t")

# x is a character vector OR a column in a data.table
cr_ind[, c("modality","feature") := tstrsplit(phenotype_id, ":", fixed = TRUE, keep = c(1,2))]

# gene id = first ENSG... token in the feature part
cr_ind[, gene_id := sub("^(ENSG[0-9]+).*", "\\1", feature)]

unique(cr_ind$modality)
length(unique(cr_ind$gene_id))

###---------------- ZNF91 ecxpression vs age
expr <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression.bed")
class(expr)

dim(expr)
View(expr[1:5, 1:10])
names(expr)[1:5]


expr_znf91 <- expr %>%
  as.data.frame() %>%
  dplyr::filter(phenotype_id == "ENSG00000167232") %>%
  dplyr::select(-c(1:3)) %>%
  data.table::transpose(keep.names = "human_id")

expr_znf91 <- expr_znf91[-1, ]
colnames(expr_znf91) <- c("human_id", "ZNF91")
head(expr_znf91)


mata_date <- fread("src/metadata_geno.tsv")
mata_date$human_id <- as.character(mata_date$human_id)
dim(mata_date)
head(mata_date)
length(unique(mata_date$human_id))

expr_znf91_age <- merge(mata_date[, -1], expr_znf91, by = "human_id")
expr_znf91_age$ZNF91 <- as.numeric(expr_znf91_age$ZNF91)

head(expr_znf91_age)
expr_znf91_age <- expr_znf91_age[order(expr_znf91_age$ZNF91, decreasing = TRUE), ]

# Fit linear regression: ZNF91 ~ age
lm_znf91 <- lm(as.numeric(ZNF91) ~ age, data = expr_znf91_age)
lm_sum   <- summary(lm_znf91)
lm_beta  <- round(coef(lm_znf91)["age"], 4)
lm_pval  <- round(coef(lm_sum)["age", "Pr(>|t|)"], 4)
lm_r2    <- round(lm_sum$r.squared, 3)
lm_label <- paste0("beta = ", lm_beta, "\np = ", lm_pval, "\nR2 = ", lm_r2)

# > lm_label
# [1] "β = -0.2289\np = 0.2083\nR² = 0.015"

# [1] "β = -0.0991\np = 0.4853\nR² = 0.005"

# Plot ZNF91 expression vs age
ggplot(expr_znf91_age, aes(x = age, y = as.numeric(ZNF91))) + # [-1, ]
  geom_point(color = "deepskyblue4", size = 1, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = lm_label,
            hjust = 1.1, vjust = 1.5, size = 2.5, color = "black") +
  labs(
    title = "ZNF91 Expression vs Age in Vestibular Samples",
    x = "Age",
    y = "ZNF91 Expression (tpm)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    plot.title = element_blank(),
  )

ggsave(
  last_plot(),
  filename = "figures/znf91_expression_vs_age_vestibular.pdf",
  units = "mm",
  width = 80,
  height = 80,
  dpi = 300
)

# Plot ZNF91 expression vs sex
ggplot(expr_znf91_age, aes(x = sex, y = as.numeric(ZNF91))) +
  geom_point(color = "deepskyblue4", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "ZNF91 Expression vs Sex in Vestibular Samples",
    x = "Sex",
    y = "ZNF91 Expression (tpm)"
  ) +
  theme_classic()

# Plot ZNF91 expression vs sex - violin plot
ggplot(expr_znf91_age, aes(x = sex, y = as.numeric(ZNF91), fill = sex)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
  labs(
    title = "ZNF91 Expression vs Sex in Vestibular Samples",
    x = "Sex",
    y = "ZNF91 Expression (tpm)"
  ) +
  theme_classic() +
  theme(legend.position = "none") 

ggsave(
  last_plot(),
  filename = "figures/znf91_expression_vs_sex_vestibular.png",
  width = 9,
  height = 7.5,
  dpi = 300
)


###-------------------------------> looking at ZNF91 isoforms ratio for different genotypes
ENSG00000167232__ENST00000300619
ENSG00000167232__ENST00000593341

### Load the weights
load("Pantry/pheast_no_outlier/intermediate/vestibular/twas/isoforms/ENSG00000167232__ENST00000300619.wgt.RDat")
View(wgt.matrix)
View(snps)
View(cv.performance)
View(hsq)
View(hsq.pv)
View(N.tot)


### Load isoform expression data
isoform_expr <- fread("Pantry/phenotyping/output/vestibular/isoforms.bed.gz")

# Extract isoform expression for ENST00000300619
isoform_expr_cln <- isoform_expr %>%
  dplyr::select(-c(`#chr`, start, end)) %>%
  dplyr::filter(phenotype_id == "ENSG00000167232__ENST00000300619") %>%
  data.table::transpose(keep.names = "human_id") %>%
  { data.table::setnames(., as.character(unlist(.[1, ]))); .[-1, ] } %>%
  dplyr::rename(sample = phenotype_id)


View(isoform_expr_cln)
head(isoform_expr_cln)


### Extract genotype for SNP 19_23307401_A_G
system("plink2 \
  --bfile Pantry/pheast/input/vestibular/vestibular_genotype \
  --snp 19_23307401_A_G \
  --export vcf \
  --out  Pantry/pheast/input/vestibular/19_23307401_A_G_genotype"
  )


# Read vcf
vcf <- fread("Pantry/pheast/input/vestibular/19_23307401_A_G_genotype.vcf", skip = "#CHROM")
gt <- tstrsplit(vcf$FORMAT, ":", fixed = TRUE)

geno <- t(vcf[, -(1:9)])   # samples only
colnames(geno) <- vcf$ID

geno_df <- data.frame(
  sample = colnames(vcf)[10:ncol(vcf)],
  genotype = geno[,1]
)

ref <- vcf$REF
alt <- vcf$ALT

geno_df$geno_alleles <- dplyr::case_when(
  geno_df$genotype == "0/0" ~ paste0(ref, ref),
  geno_df$genotype == "0/1" ~ paste0(ref, alt),
  geno_df$genotype == "1/1" ~ paste0(alt, alt),
  TRUE ~ NA_character_
)
geno_df$sample <- gsub("^FTL_", "", geno_df$sample)
View(geno_df)

### Merging isoform expression with genotype
# Check matching sample IDs
table(isoform_expr_cln$sample %in% geno_df$sample) #

# Merge
isoform_expr_geno <- merge(geno_df, isoform_expr_cln, by = "sample")
class(isoform_expr_geno)
View(isoform_expr_geno)
head(isoform_expr_geno)


names(isoform_expr_geno)[4] <- "ENST00000300619"
names(isoform_expr_geno)


ggplot(isoform_expr_geno, aes(x = genotype, y = as.numeric(ENST00000300619), fill = genotype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
  labs(
    title = "ENST00000300619 Expression vs Genotype at 19_23307401_A_G",
    x = "genotype",
    y = "ENST00000300619 Expression (tpm)"
  ) +
  theme_classic() +
  theme(legend.position = "none") 

ggsave(
  last_plot(),
  filename = "figures/znf91_expression_vs_sex_vestibular.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

###----------------- Summarize twas results for all significant genes
library(readr)
twas <- fread("Pantry/pheast_no_outlier/twas/vestibular/output/twas_hits_annot.tsv")
head(twas)
View(twas[1:5, ])
length(unique(twas$ID))
length(unique(twas$ENSG))
