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
  library(rtracklayer)
  library(jsonlite)
  library(matrixStats)
  library(sva)
})

###--------------------
# Plotting PCA of GTEx + Vestibular expression data
###--------------------

## Making a vector contains protein coding and lncRNAs gene IDs to filter expression data
# Reading GTF file
gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
  dplyr::filter(type == "gene") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
  dplyr::pull("gene_id")

# Load GTEx colors from JSON
gtex_colors <- fromJSON("https://raw.githubusercontent.com/broadinstitute/gtex-viz/master/boxplot/dev/colors.json")
# gtex_colors <- fromJSON("src/gtex_colors.json")

# Convert to data frame
gtex_annot <- bind_rows(lapply(names(gtex_colors), function(tissue) {
  data.frame(
    tissue_name = tissue,
    tissue_abbrv = gtex_colors[[tissue]]$tissue_abbrv,
    tissue_color_hex = paste0("#", gtex_colors[[tissue]]$tissue_color_hex),
    stringsAsFactors = FALSE
  )
}))

# Clean tissue names
gtex_annot <- gtex_annot %>%
  mutate(
    tissue_name_clean = tissue_name %>%
      str_to_lower() %>%
      str_replace_all(c(" - " = "_", " \\(" = "_", "\\)" = "", " " = "_", "cells_transformed" = "cells_cultured")) # Cells - Transformed fibroblasts == Cells - Cultured fibroblasts, so replacing it to match with tissue id in GTEx v8
  )

# Making a named vector for tissue abbreviation lookup
gtex_abb <- setNames(gtex_annot$tissue_name_clean, gtex_annot$tissue_abbrv)

# Need to remove some tissues with low sample size
# BLDDER ("bladder"), CVXECT ("cervix_ectocervix"), CVSEND ("cervix_endocervix"), FLLPNT ("fallopian_tube"), KDNMDL ("kidney_medulla"), TESTIS ("testis")
indx = which(names(gtex_abb) %in% c("TESTIS", "BLDDER", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL"))
gtex_abb <- gtex_abb[-indx]



## Reading GTEx tissue raw count data
files_list <- list.files("src/GTEx/tissues_raw_count", pattern = "gene_reads")
# Removing files with low sample size tissues
patterns <- c("testis", "bladder", "cervix_ectocervix", "cervix_endocervix", "fallopian_tube",  "kidney_medulla")
files_list_filtered <- files_list[!grepl(paste(patterns, collapse = "|"), files_list)]

# Clean tissue names
tissue_names <- str_extract(files_list_filtered, "(?<=v8_)[^\\.]+")

# Check matching with gtex_abb from JSON
table(tissue_names %in% gtex_abb)
table(gtex_abb %in% tissue_names)
gtex_abb[!gtex_abb %in% tissue_names]


# Reading and merging all GTEx tissue count data
# l1 = length(tissue_names)
for (i in seq_along(tissue_names)) {
  tissue = tissue_names[i]
  cat("Tissue:", tissue, "- loop No. =", i, "\n")
  
  # gtex_tissue_count <- fread(paste0("src/GTEx/tissues_raw_count/", files_list_filtered[i]))[, -c(1,3)]
  # gtex_tissue_count[1:5, 1:5]
  # colnames(gtex_tissue_count)[1] <- "gene_id"


  gtex_tissue_count <- fread(paste0("src/GTEx/tissues_raw_count/", files_list_filtered[i])) %>%
    dplyr::select(-id, -Description) %>%
    dplyr::rename(gene_id = Name) # %>%
    # dplyr::mutate(gene_id = gsub("\\..*", "", gene_id)) %>%
    # dplyr::filter(gene_id %in% genes_protien_lncRNA)
  
  l2 = dim(gtex_tissue_count)[2]
  abb <- names(which(gtex_abb == tissue))
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

dim(gtex_tissue_count_all) # Check dimensions
gtex_tissue_count_all[1:5, 1:5] # Check data
table(annot_gtx_tissue_all$tissue) # Check tissue annotation
table(annot_gtx_tissue_all$tissue)[table(annot_gtx_tissue_all$tissue) < 80]

## Reading vestibualar count data
# vestibular_expr_count <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression_count.bed")[, -1:-3]
# colnames(vestibular_expr_count)[1] <- "gene_id"

vestibular_expr_count <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression_count.bed") %>%
    dplyr::select(-`#chr`, -start, -end, -"31947386", -"31907773", -"32073088") %>% # Removing outlier samples based on PCA and metadata
    dplyr::rename(gene_id = phenotype_id) %>%
    filter(gene_id %in% genes_protien_lncRNA)
dim(vestibular_expr_count) # Check dimensions

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

# Check matching
table(rownames(meta_unique) == colnames(vestibular_expr_count)[-1]) 

# Making final annotation for gtex and vestibular data
annot_gtex_vest <- rbind(annot_gtx_tissue_all, meta_unique)
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

# Check
table(rownames(count_mtx_gtex_vest) %in% genes_protien_lncRNA) # Should be TRUE

## Filtering low expression genes and Normalization
# First round of filtering: Drop all-zero genes immediately
counts <- count_mtx_gtex_vest[rowSums(count_mtx_gtex_vest) > 0, , drop = FALSE]

# Normalization using TMM from edgeR (counts: genes x samples (integers))
stopifnot(is.matrix(counts)) # Ensure counts is a matrix
dge <- DGEList(counts = counts) # Build DGEList object
dge <- calcNormFactors(dge, method = "TMM") # Compute TMM normalization factors
logCPM <- cpm(dge, log = TRUE, prior.count = 1) # Get normalized log2 CPM matrix
save(logCPM, file = "figures/data_tables/gtex_vest_logCPM_normalized_already_filtered.RData") # Save normalized logCPM matrix
# save(logCPM, file = "figures/data_tables/gtex_vest_logCPM_normalized.RData") # Save normalized logCPM matrix


# Second round of filtering: Filter by SD or select top genes before PCA (using matrixStats library)
row_sd <- rowSds(logCPM)
keep <- row_sd >= quantile(row_sd, 0.5)
logCPM_sd <- logCPM[keep, , drop = FALSE]

## PCA
pca <- prcomp(t(logCPM_sd), scale. = TRUE)
save(pca, file = "figures/data_tables/gtex_vest_logCPM_normalized_already_filtered_PCs.RData")
# save(pca, file = "figures/data_tables/gtex_vest_logCPM_normalized_PCs.RData")
# load("figures/data_tables/gtex_vest_logCPM_normalized_PCs.RData")
load("figures/data_tables/gtex_vest_logCPM_normalized_already_filtered_PCs.RData")
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

# Preparing PC data frame for plotting
pc_df <- data.frame(human_id = rownames(pca$x), pca$x, check.names = FALSE)


## Adding meta data, batch and tissue
pc_df_annot <- merge(
  pc_df,
  cbind(human_id = rownames(annot_gtex_vest), annot_gtex_vest),
  by = "human_id",
  sort = FALSE
)


# Add VESTIB color
vestib_color <- "deepskyblue4"  # black — common choice for new tissue
gtex_palette <- setNames(gtex_annot$tissue_color_hex, gtex_annot$tissue_abbrv)
gtex_palette <- c(gtex_palette, VESTIB = vestib_color)

# Adding tissue colors
pc_df_annot <- pc_df_annot %>%
  mutate(tissue_color_hex = gtex_palette[tissue])


# Check
unique(pc_df_annot$tissue_color_hex[pc_df_annot$tissue == "VESTIB"]) # should show "#000000"
pc_df_annot[1:5, c("human_id", "tissue", "batch", "tissue_color_hex")] 

# View(pc_df_annot[pc_df_annot$tissue == "VESTIB", c("human_id", "tissue", "batch", "tissue_color_hex")] )

# indx <- which(pc_df_annot$human_id %in% c("31947386", "31907773", "32073088"))

pc_df_annot[1:5, 1:5]
# Plot by tissue
ggplot(pc_df_annot, aes(PC1, PC2, color = tissue)) + # pc_df_annot[-indx, ]
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = gtex_palette, name = "Tissue") +
  labs(
    title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) + 
    theme_classic(base_size = 14)


pc_plot_df <- pc_df_annot[-indx, ] |>
  dplyr::filter(tissue != "TESTIS")

pc_plot_df = pc_df_annot


library(ggplot2)
library(dplyr)
library(ggrepel)

label_tissues <- c("VESTIB", "BRNCHA", "PTTARY")

label_df <- pc_plot_df %>%
  filter(tissue %in% label_tissues) %>%
  group_by(tissue) %>%
  summarize(
    PC1 = median(PC1, na.rm = TRUE),
    PC2 = median(PC2, na.rm = TRUE),
    .groups = "drop"
  )


pc_plot_df <- pc_df_annot |> # pc_df_annot[-indx, ]
  dplyr::filter(tissue != "TESTIS") |>
  dplyr::mutate(
    tissue = case_when(
      grepl("^BRN", tissue) ~ "BRN_TISSUES",
      TRUE ~ tissue
    )
  )
label_tissues <- c("VESTIB", "BRN_TISSUES", "PTTARY")

label_df <- pc_plot_df %>%
  filter(tissue %in% label_tissues) %>%
  group_by(tissue) %>%
  summarize(
    PC1 = median(PC1, na.rm = TRUE),
    PC2 = median(PC2, na.rm = TRUE),
    .groups = "drop"
  )

gtex_palette["BRN_TISSUES"] <- "#EEEE00"   # example purple for brain



ggplot(pc_plot_df, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = gtex_palette, name = "Tissue") +
  geom_label_repel(
    data = label_df,
    aes(label = tissue),
    color = "black",
    fill = "white",
    label.size = 0.2,
    size = 3,
    min.segment.length = 0,
    segment.color = "black",
    box.padding = 0.5,      # space around label
    point.padding = 0.4,    # space around point
    force = 2               # increase push strength
  ) +
  labs(
    # title = "PCA of logCPM expression (GTEx + Vestibular)",
    x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
      )




# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_gtex_vestibular_tissue_Final.png",
  width = 10,
  height = 10,
  dpi = 300
)


####






# + 
#   theme(
#     legend.position = "none"
#   )
# # or
# +
#   theme(
#   legend.position = "right",
#   legend.key.size = unit(0.3, "cm"),
#   legend.text = element_text(size = 6)
# )








# ### To add sample Ids
# # Subset your vestibular samples
# vestib_subset <- subset(pc_df_annot, tissue == "VESTIB")

# ggplot(pc_df_annot[-indx, ], aes(PC1, PC2, color = tissue)) +
#   geom_point(size = 3, alpha = 0.9) +

#   # Draw VESTIB samples again on top (larger, outlined)
#   geom_point(
#     data = vestib_subset,
#     aes(PC1, PC2),
#     size = 4,
#     shape = 21,
#     fill = "white",
#     stroke = 1.2,
#     color = "black"
#   ) +

#   # Add human_id text — only for VESTIB samples
#   geom_text(
#     data = vestib_subset,
#     aes(label = human_id),
#     size = 3,
#     vjust = -2,
#     # fontface = "bold"
#   ) +

#   scale_color_manual(values = gtex_palette, name = "Tissue") +
#   labs(
#     title = "PCA of logCPM expression (GTEx + Vestibular)",
#     x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
#     y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
#   ) +
#   theme_classic(base_size = 14)



###--------------------
# Plotting PCA Vestibular expression data - TMM normalized logCPM
###--------------------
# Reading vestibular count data
vestibular_expr_count <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression_count.bed")[, -1:-3]
counts <- as.matrix(vestibular_expr_count[, -1])
rownames(counts) <- vestibular_expr_count$phenotype_id

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

# Filtering based on SD (matrixStats)
row_sd <- rowSds(logCPM)
thr <- quantile(row_sd, 0.50)
keep_sd <- row_sd >= thr
logCPM_top <- logCPM[keep_sd, , drop = FALSE]
cat("Kept", nrow(logCPM_top), "genes out of", nrow(logCPM), "by SD ≥ median\n")

# PCA (samples as rows, genes as columns)
pca <- prcomp(t(logCPM), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2) # Explained variance for axis labels
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)

## Adding meta data
setDT(pc_df) # Ensure both are data.tables
meta_unique$human_id <- rownames(meta_unique)
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, batch)] # Inspect

# Find index of outlier sample to remove
indx <- which(pc_df_merged$human_id == "31947386") 

# Plot by Batch
ggplot(pc_df_merged[-indx], aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression - vestibular",
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
  height = 9,
  dpi = 300
)

###--------------------
# Plotting PCA Vestibular expression data - TMM normalized logCPM + Batch correction with ComBat
###--------------------

# Check that the matrix and metadata align
all(colnames(logCPM) == meta_unique$human_id)  # should be TRUE

# Just batch — no biological covariate
batch <- meta_unique$batch

# Model matrix with only intercept (no biological groups)
mod0 <- model.matrix(~ 1, data = meta_unique)

# Apply ComBat
logcpm_corrected <- ComBat(dat = logCPM, batch = batch, mod = mod0, par.prior = TRUE)

## PCA
pca <- prcomp(t(logcpm_corrected), scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
pc_df <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE)

## Adding meta data
setDT(pc_df) # Ensure both are data.tables
setDT(meta_unique) # Ensure both are data.tables
colnames(pc_df)[1] <- "human_id"
pc_df_merged <- merge(pc_df, meta_unique, by = "human_id", all.x = TRUE) # Merge metadata by human_id
pc_df_merged[1:5, .(human_id, tissue, batch)] # Inspect
indx <- which(pc_df_merged$human_id == "31947386") # Find index of outlier sample to remove

# Plot by tissue
ggplot(pc_df_merged[-indx, ], aes(x = PC1, y = PC2, color = batch, label = human_id)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of TMM + logCPM expression - vestibular - batch corrected",
    x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
    color = "Batch"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  last_plot(),
  filename = "figures/pca_expression_tmm_logcpm_vestibular_batch_corrected.png",
  width = 9,
  height = 9,
  dpi = 300
)


### END
