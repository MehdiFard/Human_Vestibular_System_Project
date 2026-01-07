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
})


###--------------------
# Plotting PCA of saccule and utricle expression data - TMM + logCPM
###--------------------

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

# Counts matrix (genes x samples)
counts = sac_utr_expr_count
stopifnot("gene_id" %in% names(counts))
samples <- setdiff(names(counts), "gene_id")
count_mtx <- as.matrix(counts[, ..samples]) # To matrix
rownames(count_mtx) <- counts$gene_id

# Create DGEList and apply TMM normalization
y <- DGEList(counts = count_mtx)
y <- calcNormFactors(y, method = "TMM") # It computes normalization factors and stores them in y$samples$norm.factors, which are then used to adjust library sizes

# Convert to log2-CPM (counts per million)
logCPM <- cpm(y, log = TRUE, prior.count = 1)

## Filtering steps before PCA
# Remove genes with any NA or non-finite values
logCPM <- logCPM[apply(logCPM, 1, function(x) all(is.finite(x))), ]

# Remove genes with zero variance across samples
gene_var <- apply(logCPM, 1, var)
logCPM <- logCPM[gene_var > 0, ]

# Remove genes with very low expression everywhere. (If you came directly from CPM/logCPM, low-expression genes can add noise)
mean_expr <- rowMeans(logCPM)
logCPM <- logCPM[mean_expr > 0.5, ]   # adjust threshold as needed

# (Optional but good) Remove outlier samples with zero variance or too many NAs
sample_var <- apply(logCPM, 2, var)
logCPM <- logCPM[, sample_var > 0, drop = FALSE]

# Check dimensions
cat("Genes remaining:", nrow(logCPM), "\nSamples remaining:", ncol(logCPM), "\n")

# Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2) # Explained variance for axis labels
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)


# Adding tissue info
pc_df$Tissue <- ifelse(grepl("_sac$", pc_df$Sample), "Saccule",
                ifelse(grepl("_utr$", pc_df$Sample), "Utricle", NA))

## Adding meta data
setDT(pc_df) # Ensure it's a data.tables
setDT(meta_sac_utr_unique) # Ensure it's a data.tables
pc_df[, human_id := as.character(sub("_(sac|utr)$", "", Sample))] # Extract numeric ID from Sample (remove "_sac" or "_utr")
pc_df_merged <- merge(pc_df, meta_sac_utr_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(Sample, human_id, Tissue, tissue, age, sex, batch)] # Inspect

# Quick check
table(pc_df$Tissue, useNA = "ifany")

# outlier sample index
indx <- which(pc_df_merged$human_id %in% c("31947386", "31907773", "32073088"))


# Plot by tissue
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression (Saccule and Utricle)",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_sac_utr_tissue.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

# Plot by batch
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = batch)) + # [-indx, ]
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression (Saccule and Utricle)",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_sac_utr_batch.png",
  width = 9,
  height = 7.5,
  dpi = 300
)


###-------------------------- Running PCA after within-individual mean-centering or paired normalization.

# Extract everything before the first underscore
col_ids <- sub("_.*$", "", colnames(logCPM))

# Common IDs
paired_ids <- intersect(sac_name, utr_name) 

# Boolean vector for the columns to keep
keep_cols <- col_ids %in% paired_ids
logCPM_paired <- logCPM[, keep_cols]

# logCPM is a matrix with genes as rows and samples as cols
logcpm_df <- logCPM_paired %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = c("human_id", "tissue"),
    names_sep = "_"
  )


# Calculating mean and centered values for each gene within each sample
logcpm_centered <- logcpm_df %>%
  group_by(gene_id, human_id) %>%
  mutate(
    mean_expr = mean(value, na.rm = TRUE),
    value_centered = value - mean_expr
  ) %>%
  ungroup()


# Convert back to wide matrix (genes × samples)
logCPM_centered <- logcpm_centered %>%
  unite("sample", human_id, tissue, sep = "_") %>%
  dplyr::select(gene_id, sample, value_centered) %>%
  tidyr::pivot_wider(names_from = sample, values_from = value_centered) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()


# Run PCA
pca <- prcomp(t(logCPM_centered), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2) # Explained variance for axis labels
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)


# Adding tissue info
pc_df$Tissue <- ifelse(grepl("_sac$", pc_df$Sample), "Saccule",
                ifelse(grepl("_utr$", pc_df$Sample), "Utricle", NA))

## Adding meta data
setDT(pc_df) # Ensure it's a data.tables
setDT(meta_sac_utr_unique) # Ensure it's a data.tables
pc_df[, human_id := as.character(sub("_(sac|utr)$", "", Sample))] # Extract numeric ID from Sample (remove "_sac" or "_utr")
pc_df_merged <- merge(pc_df, meta_sac_utr_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(Sample, human_id, Tissue, tissue, age, sex, batch)] # Inspect

# Quick check
table(pc_df$Tissue, useNA = "ifany")

# outlier sample index
indx <- which(pc_df_merged$human_id == "31947386")
indx <- which(pc_df_merged$human_id %in% c("31947386", "31907773", "32073088"))


# Plot by tissue
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression (Saccule and Utricle)",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_sac_utr_tissue_mean_centered.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

# Plot by batch
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression (Saccule and Utricle)",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_sac_utr_batch_mean_centered.png",
  width = 9,
  height = 7.5,
  dpi = 300
)

### END