#!bin/bash Rscript

###---> PCA on expression result from GTEx and Vestibular tissue

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
    library(readr)
})

### Assembling Kallisto read counts into gene expression table for vestibular and GTEx tissues
# First making transcript to gene table
library(tximport)
gtf_file <- "Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf"
txdb <- makeTxDbFromGFF(file = gtf_file)
# saveDb(x=txdb, file = "gencode.v32.annotation.TxDb")
keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, columns = "GENEID", keytype = "TXNAME")
tx2gene <- unique(tx2gene[, c("TXNAME", "GENEID")])
colnames(tx2gene) <- c("tx_id", "gene_id")

# Saving it
write_tsv(tx2gene, "src/tx2gene.tsv")

# Assembling vestibular samples expression 
base_dir <- "Pantry/phenotyping/intermediate/vestibular/expression"

# One folder per sample, each containing abundance.h5
files <- list.files(base_dir,
                    pattern = "abundance\\.h5$",
                    recursive = TRUE,
                    full.names = TRUE)

# Sample id = parent folder name
names(files) <- basename(dirname(files))

# Optional: make sure names are unique
stopifnot(!anyDuplicated(names(files)))

# Getting the counts
txi_vest <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene,          # from your TxDb code
                countsFromAbundance = "no")  # keeps est_counts

# Expression matrix for vestibular samples (raw counts)
gene_counts_vestib <- txi_vest$counts   # gene-level matrix (genes x samples)
fwrite(as.data.frame(gene_counts_vestib) %>% rownames_to_column("gene_id"), "src/gene_counts_vestib.tsv", sep = "\t")
# write_tsv(as.data.frame(gene_counts_vestib) %>% rownames_to_column("gene_id"), "src/gene_counts_vestib.tsv.gz")

### Preping GTEx raw counts from Dan's pantry for comparison
gtf_file <- "Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf"
base_dir <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/phenos"
out_dir  <- "src/GTEx/tissues_raw_count_from_Dan_pantry"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Discover tissues
tissue_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
tissues <- setdiff(basename(tissue_dirs), c("ref"))

# Process each tissue
for (tissue in tissues) {
    tissue_expr_dir <- file.path(base_dir, tissue, "intermediate", "expression")
    if (!dir.exists(tissue_expr_dir)) next

    # Find all abundance.tsv files under sample_id folders
    files <- list.files(
        tissue_expr_dir,
        pattern = "^abundance\\.h5$",
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(files) == 0) {
        message("No abundance.tsv found for tissue: ", tissue)
        next
    }

    # sample id = parent folder name
    names(files) <- basename(dirname(files))

    # Safety: ensure unique sample ids
    if (anyDuplicated(names(files))) {
        dup <- names(files)[duplicated(names(files))]
        stop("Duplicate sample IDs found in tissue ", tissue, ": ", paste(unique(dup), collapse = ", "))
    }

    # Import and summarize to gene-level raw counts (sum est_counts)
    txi <- tximport(
        files,
        type = "kallisto",
        tx2gene = tx2gene,
        countsFromAbundance = "no"
    )

    counts <- txi$counts  # gene x sample
    out_dt <- data.table(gene_id = rownames(counts), counts)
    # Optional: remove rownames attribute
    rownames(out_dt) <- NULL

    # Save (tsv). Name includes tissue.
    out_file <- file.path(out_dir, paste0("gene_counts_", tissue, ".tsv.gz"))
    fwrite(out_dt, out_file, sep = "\t")

    message("Saved: ", out_file, " (genes=", nrow(out_dt), ", samples=", ncol(out_dt)-1, ")")
}


###--------------------- Now running PCA
# Using gtf file to make a vector contains protein coding and lncRNAs gene IDs to filter expression data
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

# Need to remove some tissues with low sample size (> 20): CVXECT ("cervix_ectocervix"), CVSEND ("cervix_endocervix"), FLLPNT ("fallopian_tube"), KDNMDL ("kidney_medulla"), TESTIS ("testis")
indx = which(names(gtex_abb) %in% c("TESTIS", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL"))
gtex_abb <- gtex_abb[-indx]


## Reading GTEx tissue raw count data
files_list <- list.files("src/GTEx/tissues_raw_count_from_Dan_pantry/", pattern = ".tsv.gz$")
# Removing tissues with low sample size (> 20)
patterns <- c("TESTIS", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL")
# patterns <- c("testis", "bladder", "cervix_ectocervix", "cervix_endocervix", "fallopian_tube",  "kidney_medulla")
files_list_filtered <- files_list[!grepl(paste(patterns, collapse = "|"), files_list)]

# Clean tissue names from file names
tissue_names <- sub("^gene_counts_(.*)\\.tsv\\.gz$", "\\1", files_list_filtered)

# Check matching with gtex_abb from JSON
table(names(gtex_abb) %in% tissue_names)
names(gtex_abb)[!names(gtex_abb) %in% tissue_names]
gtex_abb[!names(gtex_abb) %in% tissue_names]


# Reading and merging all GTEx tissue count data
path_to_gtex_expr <- "src/GTEx/tissues_raw_count_from_Dan_pantry/"

# Loop through each tissue file, read the count data, and merge into a single data frame
for (i in seq_along(tissue_names)) {
    tissue = tissue_names[i]
    cat("Tissue:", tissue, "- loop No. =", i, "\n")
    
    # reading GTEx tissue count data
    gtex_tissue_count <- fread(paste0(path_to_gtex_expr, "gene_counts_", tissue, ".tsv.gz"))
    
    # Annotation for GTEx tissue samples
    l = dim(gtex_tissue_count)[2]
    annot_gtx_tissue <- data.frame(tissue = rep(tissue, l-1), row.names = colnames(gtex_tissue_count)[-1])
    
    # Save it
    if(i == 1) {
        gtex_tissue_count_all = gtex_tissue_count
        annot_gtx_tissue_all = annot_gtx_tissue
    } else {
        gtex_tissue_count_all <- merge(gtex_tissue_count_all, gtex_tissue_count, by = "gene_id", sort = F)
        annot_gtx_tissue_all <- rbind(annot_gtx_tissue_all, annot_gtx_tissue)
    }
}
annot_gtx_tissue_all$batch <- "gtex" # Adding batch info

gtex_tissue_count_all[1:3, 1:3]
head(annot_gtx_tissue_all)

## Reading vestibualar count data
count_vestib = fread("src/gene_counts_vestib.tsv", header = TRUE)
vestibular_expr_count <- fread("src/gene_counts_vestib.tsv", header = TRUE) %>%
    dplyr::filter(gene_id %in% genes_protien_lncRNA)

## Annotation - read metadata for batch info
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

dim(annot_gtex_vest)
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
    # collapse_by_gene() %>%
    inner_join(collapse_by_gene(vestibular_expr_count), by = "gene_id") %>%
    column_to_rownames("gene_id") %>%
    as.matrix()

# Check
dim(count_mtx_gtex_vest)
count_mtx_gtex_vest[1:3, 1:3]
table(rownames(count_mtx_gtex_vest) %in% genes_protien_lncRNA) # Should be TRUE

## Filtering low expression genes and Normalization
# First round of filtering: Drop all-zero genes immediately
counts <- count_mtx_gtex_vest[rowSums(count_mtx_gtex_vest) > 0, , drop = FALSE]

# Normalization using TMM from edgeR (counts: genes x samples (integers))
stopifnot(is.matrix(counts)) # Ensure counts is a matrix
dge <- DGEList(counts = counts) # Build DGEList object
dge <- calcNormFactors(dge, method = "TMM") # Compute TMM normalization factors
logCPM <- cpm(dge, log = TRUE, prior.count = 1) # Get normalized log2 CPM matrix
save(logCPM, file = "figures/data_tables/gtex_vest_logCPM_normalized_vestib_gtex_same_pipeline.RData") # Save normalized logCPM matrix
logCPM <- load("figures/data_tables/gtex_vest_logCPM_normalized_vestib_gtex_same_pipeline.RData") # Load it back when needed

# Second round of filtering: Filter by SD or select top genes before PCA (using matrixStats library)
row_sd <- rowSds(logCPM)
keep <- row_sd >= quantile(row_sd, 0.75); table(keep)
logCPM_sd <- logCPM[keep, , drop = FALSE]
dim(logCPM_sd)


## PCA
pca <- prcomp(t(logCPM_sd), scale. = TRUE)
# save(pca, file = "figures/data_tables/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData")
load("figures/data_tables/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData")
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
unique(pc_df_annot$tissue_color_hex[pc_df_annot$tissue == "VESTIB"]) # should show "deepskyblue4"
pc_df_annot[1:5, c("human_id", "tissue", "batch", "tissue_color_hex")] 

# Let's keep out the outliers
indx <- which(pc_df_annot$human_id %in% c("31907773", "32073088")) # "31947386"

# Plot by tissue
ggplot(pc_df_annot[-indx, ], aes(PC1, PC2, color = tissue)) + # pc_df_annot[-indx, ]
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = gtex_palette, name = "Tissue") +
    labs(
        title = "PCA of logCPM expression (GTEx + Vestibular)",
        x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
        y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)")
    ) + 
    theme_classic(base_size = 14)


# A better plot
library(grid)

# Making df
pc_plot_df <- pc_df_annot |>
    # dplyr::filter(tissue != "TESTIS") |>
    dplyr::mutate(
        tissue = case_when(
        grepl("^BRN", tissue) ~ "BRN_TISSUES",
        TRUE ~ tissue
        )
    )

# To lable tissues
label_tissues <- c("VESTIB", "BRN_TISSUES", "PTTARY", "WHLBLD", "LCL", "SPLEEN", "GEUVADIS",
                    "PANCREAS","ADRENAL", "THYROID")

# Calculate median PC1 and PC2 for each tissue to position labels, and add nudges to reduce overlap
label_df <- pc_plot_df %>%
    filter(tissue %in% label_tissues) %>%
    group_by(tissue) %>%
    summarize(
        PC1 = median(PC1, na.rm = TRUE),
        PC2 = median(PC2, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # optional: manual nudges to reduce overlap (tune as needed)
    mutate(
        nudge_x = case_when(
        tissue == "VESTIB" ~ 0.5,
        tissue == "BRN_TISSUES" ~ -0.5,
        TRUE ~ 0
        ),
        nudge_y = case_when(
        tissue == "WHLBLD" ~ 0.5,
        tissue == "LCL" ~ -0.5,
        TRUE ~ 0
        )
    )


# Add a single color for all brain tissues
gtex_palette["BRN_TISSUES"] <- "#EEEE00"

# Plot
ggplot(pc_plot_df[-indx, ], aes(PC1, PC2, color = tissue)) +
    geom_point(size = 1, alpha = 0.9) +
    scale_color_manual(values = gtex_palette, name = "Tissue") +
    geom_label(
        data = label_df,
        aes(x = PC1 + nudge_x, y = PC2 + nudge_y, label = tissue),
        inherit.aes = FALSE,
        color = "black",
        fill = "white",
        label.size = 0.08,
        size = 1.5,
        label.padding = unit(0.3, "lines")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, col="grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, col="grey50") +
    labs(
        x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    )



###--------------------------- PCA on each tissue using centroid method
# Align annotation to matrix columns
stopifnot(all(colnames(count_mtx_gtex_vest) %in% rownames(annot_gtex_vest)))
annot <- annot_gtex_vest[colnames(count_mtx_gtex_vest), , drop = FALSE]

# edgeR object from kallisto est_counts
dge <- DGEList(counts = count_mtx_gtex_vest)
group <- factor(annot[, "tissue"])

# Filter lowly expressed genes (edgeR-recommended)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization (works fine for fractional counts)
dge <- calcNormFactors(dge, method = "TMM")

# logCPM
logcpm <- cpm(dge, log = TRUE, prior.count = 1)  # genes x samples

# drop zero-variance genes (rare after filtering)
logcpm <- logcpm[rowVars(logcpm) > 0, , drop = FALSE]

# Tissue centroids: genes x tissues
tissues <- unique(annot[, "tissue"])

centroids <- vapply(
    tissues,
    FUN = function(tt) {
        cols <- which(annot[, "tissue"] == tt)
        rowMeans(logcpm[, cols, drop = FALSE])
    },
    FUN.VALUE = numeric(nrow(logcpm))
)
rownames(centroids) <- rownames(logcpm)
colnames(centroids) <- tissues

# Filter genes by SD across tissue centroids (keep top 25% most variable)
row_sd <- rowSds(centroids)
summary(row_sd)
keep <- row_sd >= quantile(row_sd, 0.75); table(keep)
centroids_sd <- centroids[keep, , drop = FALSE]
dim(centroids_sd)


# PCA (tissues as observations)
pca <- prcomp(t(centroids_sd), center = TRUE, scale. = FALSE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- data.frame(
    tissue = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
)

# If each tissue belongs to a single batch, you can attach it like this:
batch_by_tissue <- as.data.frame(annot) %>% distinct(tissue, batch)
pca_df <- pca_df %>% left_join(batch_by_tissue, by = "tissue")

# 5) Plot
ggplot(pca_df, aes(PC1, PC2, color = tissue)) +
    geom_point(size = 1, alpha = 0.9) +
    scale_color_manual(values = gtex_palette, name = "Tissue") +
    # label only selected tissues
    geom_text(
        data = subset(pca_df, tissue %in% tissues_to_label),
        aes(label = tissue),
        size = 1.8,
        vjust = -1,
        color = "black"
    ) +

    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, col="grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, col="grey50") +
    labs(
        x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    )


    # Save plot
    ggsave(
    last_plot(),
    filename = "figures/pca_expression_gtex_vestibular_tissue_Final_cln_vestib_gtex_same_pipeline_centroid.pdf",
    width = 5,
    height = 5,
    dpi = 300
)


### END