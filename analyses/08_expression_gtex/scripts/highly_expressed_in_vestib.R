### Comparing gene expression in Vestibular tissue with GTEx and finding genes that are highly expressed in vestibular tissue
rm(list=ls())

# Load libraries
library(data.table)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(rtracklayer)
library(jsonlite)


# define pseudocount
pseudocount = 1

# Getting protein coding and lncRNAs gene IDs
gtf <- readGFF("resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
    dplyr::filter(type == "gene") %>% 
    dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
    dplyr::pull("gene_id")
cfg <- yaml::read_yaml("config/config.yaml")   # + lab-specific values (config_lab.yaml, not in git)
if (file.exists("config/config_lab.yaml")) cfg <- modifyList(cfg, yaml::read_yaml("config/config_lab.yaml"))

# Vestibular tissue, calculating median gene expression
expr_vest_tpm <- fread("analyses/05_pantry/phenotyping/output/vestibular/unnorm/expression.bed") %>%
    dplyr::select(-`#chr`, -start, -end, -any_of(unlist(cfg$samples$outlier_samples))) %>%
    dplyr::rename(gene_id = phenotype_id) %>%
    filter(gene_id %in% genes_protien_lncRNA)


# Convert numeric columns to matrix and compute row medians and filtering out genes with median TPM <= 0.1
median_expr_vest <- data.table(
    gene_id = expr_vest_tpm$gene_id,
    VESTIB = rowMedians(as.matrix(expr_vest_tpm[, -1]))
) # [VESTIB > 0.1]


# GTEx tissues median gene expression
path_to_gtex_expr <- paste0(cfg$refs$gtex_pantry_phenos, "/")
# Getting tissue names excluding "ref" and low sample size
tissue_names <- setdiff(list.dirs(path_to_gtex_expr, recursive = FALSE, full.names = FALSE), c("ref", "TESTIS", "BLDDER", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL"))


# A df for storing the loop output
median_expr_tbl = median_expr_vest

# Looping across files
for (i in seq_along(tissue_names)) {
    tissue = tissue_names[i]
    cat("Looking at:", tissue, "Tissue", "\n")
    expr_gtex <- fread(paste0(path_to_gtex_expr, tissue,"/output/unnorm/expression.bed")) %>%
        dplyr::select(-`#chr`, -start, -end) %>%
        dplyr::rename(gene_id = phenotype_id) %>%
        filter(gene_id %in% genes_protien_lncRNA)

    # Convert numeric columns to matrix and compute row medians
    median_expr_gtex <- data.table(
        gene_id = expr_gtex$gene_id,
        gtex_tissue = rowMedians(as.matrix(expr_gtex[, -1]))
        ) # [gtex_tissue > 0.1]

    # Compare median expression between vestibular and GTEx tissue
    col = paste0(tissue)
    median_expr_tbl <- median_expr_tbl %>%
        inner_join(median_expr_gtex, by = "gene_id") %>%
        # remove genes basically off in both tissues
        filter(VESTIB >= 0.1 | gtex_tissue >= 0.1) %>%
        dplyr::rename(!!col := gtex_tissue)

}


# Filtering - Keeping the genes if the tpm expression is higher than 0.1 in at least 10 tissues
median_expr_tbl_filt <- median_expr_tbl %>%
    mutate(n_tissues_expr = rowSums(across(-gene_id, ~ . > 0.1))) %>%
    filter(n_tissues_expr >= 10) %>%
    select(-n_tissues_expr)


# Tissue names from expr table
tissue_cols <- setdiff(names(median_expr_tbl), c("gene_id", "VESTIB"))


# Then compute max tissue (with log2 transform)
dt_max <- median_expr_tbl_filt[, {

    # log2 transform all tissues
    vals_log <- log2(c(VESTIB, unlist(.SD, use.names = FALSE)) + 1)
    
    vestib_log <- vals_log[1]
    others_log <- vals_log[-1]
    
    j <- which.max(others_log)
    
    max_other_log <- others_log[j]
    
    diff <- vestib_log - max_other_log   # log2 difference
    lfc  <- diff                         # already log2 fold-change
    
    .(
        VESTIB = vestib_log,
        max_other = max_other_log,
        max_tissue = tissue_cols[j],
        diff = diff,
        log2FC = lfc
        )

}, by = gene_id, .SDcols = tissue_cols]


# Adding t-test result
ttest <- fread("analyses/08_expression_gtex/results/vestib_vs_gtex_t_tests.tsv")[, c(1,8)]
head(ttest)
table(ttest$padj_two_sided < 0.05)
dt_max_ttest <- dt_max %>%
        left_join(ttest, by = "gene_id")


# remove genes where max_other is still zero
dt_max_ttest_cln <- dt_max_ttest[max_other > 0]


### Plot expression of a gene across all tissues
# Getting GTEx tissue color
gtex_colors <- fromJSON("https://raw.githubusercontent.com/broadinstitute/gtex-viz/master/boxplot/dev/colors.json")

# Convert to data frame
gtex_annot <- bind_rows(lapply(names(gtex_colors), function(tissue) {
    data.frame(
        tissue_name = tissue,
        tissue_abbrv = gtex_colors[[tissue]]$tissue_abbrv,
        tissue_color_hex = paste0("#", gtex_colors[[tissue]]$tissue_color_hex),
        stringsAsFactors = FALSE
    )
}))

# Making a vector of tissues and colors
gtex_palette <- setNames(gtex_annot$tissue_color_hex, gtex_annot$tissue_abbrv)
vestib_color <- "deepskyblue4"  # black — common choice for new tissue
gtex_palette <- c( VESTIB = vestib_color, gtex_palette)

# Adding tissue colors
dt_max_ttest_cln[, tissue_color := gtex_palette[max_tissue]]


# Keeping significant genes
dt_plot <- dt_max_ttest_cln %>%
    dplyr::filter(padj_two_sided < 0.01)

# Plot
ggplot(dt_plot, aes(x = max_other, y = VESTIB, color = -log10(padj_two_sided))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 0.5) +
    # scale_color_viridis_c(name = "-log10(p)") +
    scale_color_viridis_c(name = expression(-log[10](p))) +
    guides(color = guide_colorbar(barheight = grid::unit(8, "mm"),
                                    barwidth = grid::unit(2, "mm"),
                                    ticks = FALSE,
                                    title.position = "top")) +
    labs(
        x = expression("Highest median expression across GTEx tissues ("*log[2]*"tpm)"),
        y = expression("Median expression in VESTIB tissue ("*log[2]*"tpm)")
    ) +
    theme_classic() + # base_size = 6
    theme(
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = grid::unit(5, "mm"),
        legend.key.width = grid::unit(2.2, "mm"),
        legend.spacing.y = grid::unit(0.6, "mm"),
        legend.box.spacing = grid::unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = c(0.99, 0.4),
        legend.justification = c(1, 0.5)
    )


ggsave(
    last_plot(),
    filename = "analyses/08_expression_gtex/figures/vestibular_vs_gtex_max_tissue_expression_new.pdf",
    units = "mm",
    width = 90,
    height = 90,
    dpi = 300
)



###------ Getting genes with highest expression in VESTIBULAR compared to GTEx tissues to run enrichment analysis 
lfc_cut <- 2      # 4-fold
dt_max_filt <- dt_plot[log2FC > lfc_cut]
dim(dt_max_filt)
head(dt_max_filt)

# Write highly expressed genes in vestibular
write.table(dt_max_filt[, "gene_id"], "analyses/08_expression_gtex/results/vestibular_high_expressed_genes_accross_all_gtex.tsv", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)

# Write background genes
write.table(dt_plot[, "gene_id"], "analyses/08_expression_gtex/results/vestibular_and_gtex_genes_background.tsv", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)


### END