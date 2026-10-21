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
gtf <- readGFF("resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
  dplyr::filter(type == "gene") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
  dplyr::pull("gene_id")

cfg <- yaml::read_yaml("config/config.yaml")   # + lab-specific values (config_lab.yaml, not in git)
if (file.exists("config/config_lab.yaml")) cfg <- modifyList(cfg, yaml::read_yaml("config/config_lab.yaml"))

# Read Saccule and removing outliers
saccule_expr_count <- fread("analyses/05_pantry/phenotyping/output/saccule/unnorm/expression_count.bed") %>%
  dplyr::select(-any_of(c("#chr", "start", "end", unlist(cfg$samples$outlier_samples))))

sac_name <- names(saccule_expr_count) [-1]
names(saccule_expr_count)[-1] <- paste0(names(saccule_expr_count)[-1], "_sac") # Add tissue suffix

# Read Utricle and removing outliers
utricle_expr_count <- fread("analyses/05_pantry/phenotyping/output/utricle/unnorm/expression_count.bed")%>%
  dplyr::select(-any_of(c("#chr", "start", "end", unlist(cfg$samples$outlier_samples))))

utr_name <- names(utricle_expr_count)[-1]
names(utricle_expr_count)[-1] <- paste0(names(utricle_expr_count)[-1], "_utr") # Add tissue suffix


# Merge saccule and utricle and keeping protien and lncRNA genes
sac_utr_expr_count <- saccule_expr_count %>%
  inner_join(utricle_expr_count, by = "phenotype_id") %>%
  dplyr::rename(gene_id = phenotype_id) %>%
  filter(gene_id %in% genes_protien_lncRNA)
dim(sac_utr_expr_count) # Check dimensions after filtering

# Read metadata
meta <- fread("analyses/03_sample_matching/results/sample_sheets/metadata_rna.tsv")[, c("human_id", "tissue", "age", "sex", "batch")]
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

# Keep genes with SD at or above the 75th percentile across samples
gene_sd <- apply(logCPM, 1, sd)
summary(gene_sd) # Check SD distribution
sd_cutoff <- quantile(gene_sd, 0.5, na.rm = TRUE)
logCPM_cln <- logCPM[gene_sd >= sd_cutoff, ]

dim(logCPM_cln) # Check dimensions after filtering

# (Optional but good) Remove outlier samples with zero variance or too many NAs
# sample_var <- apply(logCPM_cln, 2, var)
# logCPM_cln <- logCPM_cln[, sample_var > 0, drop = FALSE]

# Check dimensions
cat("Genes remaining:", nrow(logCPM_cln), "\nSamples remaining:", ncol(logCPM_cln), "\n")

logCPM_cln[1:5, 1:5] # Inspect the logCPM matrix after filtering

# Run PCA
pca <- prcomp(t(logCPM_cln), scale. = TRUE)
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

# # outlier sample index


# Plot by tissue
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = Tissue)) +
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
  filename = "manuscript/figures/SuppFig1/panels/pca_expression_tmm_logcpm_sac_utr_tissue.png",
  units = "mm",
  width = 60,
  height = 60,
  dpi = 300
)



# Plot by tissue
p_tissue <- ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("Saccule" = "#1b9e77", "Utricle" = "#d95f02")) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.2, alpha = 1))) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = grid::unit(2, "mm"),
    legend.key.height = grid::unit(2, "mm"),
    legend.key.width = grid::unit(2, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm"),
    legend.box.spacing = grid::unit(0.2, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  ); p_tissue


ggsave(
  plot = p_tissue,
  filename = "manuscript/figures/SuppFig1/panels/pca_expression_tmm_logcpm_sac_utr_tissue.pdf",
  units = "mm",
  width = 80,
  height = 80,
  dpi = 300
)


# Plot by batch
ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) + # [-indx, ]
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
  filename = "manuscript/figures/SuppFig1/panels/pca_expression_tmm_logcpm_sac_utr_batch.png",
  width = 6,
  height = 6,
  dpi = 300
)


###-------------------------- Running PCA after within-individual mean-centering or paired normalization.

# Extract everything before the first underscore
col_ids <- sub("_.*$", "", colnames(logCPM_cln))

# Common IDs
paired_ids <- intersect(sac_name, utr_name) 

# Boolean vector for the columns to keep
keep_cols <- col_ids %in% paired_ids
logCPM_paired <- logCPM_cln[, keep_cols]

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
indx <- which(pc_df_merged$human_id %in% unlist(cfg$samples$outlier_samples))


# Plot by tissue
p_tissue_mean_centered <- ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("Saccule" = "#1b9e77", "Utricle" = "#d95f02")) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Tissue"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.2, alpha = 1))) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = grid::unit(2, "mm"),
    legend.key.height = grid::unit(2, "mm"),
    legend.key.width = grid::unit(2, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm"),
    legend.box.spacing = grid::unit(0.2, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  ); p_tissue_mean_centered

# Save plot
ggsave(
  plot = p_tissue_mean_centered,
  filename = "manuscript/figures/SuppFig1/panels/pca_sac_vs_utr_tmm_logcpm_centered.pdf",
  units = "mm",
  width = 80,
  height = 80,
  dpi = 300
)

# Plot by batch
p_batch_mean_centered <- ggplot(pc_df_merged, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 1.8, alpha = 0.9) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ); p_batch_mean_centered

# Save plot
ggsave(
  plot = p_batch_mean_centered,
  filename = "manuscript/figures/SuppFig1/panels/pca_expression_tmm_logcpm_sac_utr_batch_mean_centered.pdf",
  units = "mm",
  width = 50,
  height = 50,
  dpi = 300
)

### END