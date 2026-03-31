### Comparing gene expression in Vestibular tissue with GTEx tissue to see if there is a tissue specific gene
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
gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
    dplyr::filter(type == "gene") %>% 
    dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
    dplyr::pull("gene_id")

# Vestibular tissue, calculating median gene expression
expr_vest_tpm <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression.bed") %>%
    dplyr::select(-`#chr`, -start, -end, -"31947386", -"31907773", -"32073088") %>%
    dplyr::rename(gene_id = phenotype_id) %>%
    filter(gene_id %in% genes_protien_lncRNA)
dim(expr_vest_tpm)

# Convert numeric columns to matrix and compute row medians
median_expr_vest <- data.table(
    gene_id = expr_vest_tpm$gene_id,
    median_expr_vest = rowMedians(as.matrix(expr_vest_tpm[, -1]))
)

# GTEx tissues median gene expression
path_to_gtex_expr <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/phenos/"
# Getting tissue names excluding "ref" and "TESTIS"
tissue_names <- setdiff(list.dirs(path_to_gtex_expr, recursive = FALSE, full.names = FALSE), c("ref", "TESTIS"))

# Looping over vestibular tissues
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
        median_expr_gtex = rowMedians(as.matrix(expr_gtex[, -1])))

    # Compare median expression between vestibular and GTEx tissue
    col = paste0("median_expr_", tissue)
    expr_compare <- median_expr_vest %>%
        inner_join(median_expr_gtex, by = "gene_id") %>%
        # remove genes basically off in both tissues
        filter(median_expr_vest >= 0.1 | median_expr_gtex >= 0.1) %>%
        mutate(
        fold_change = (median_expr_vest + pseudocount) / (median_expr_gtex + pseudocount),
        log2_fc = log2((median_expr_vest + pseudocount) / (median_expr_gtex + pseudocount))) %>%
        filter(log2_fc > 1 | log2_fc < -1)  %>%
        dplyr::rename(!!col := median_expr_gtex)

    # Save comparison table
    fwrite(expr_compare,
        file = paste0("figures/expr_comparison/VESTIB_vs_", tissue, "_expr_comparison.tsv"),
        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    # Plot
    expr_compare_plot <- expr_compare %>%
    mutate(
        log2_vest = log2(median_expr_vest + pseudocount),
        log2_gtex = log2(.data[[col]] + pseudocount),
        high_fc = ifelse(abs(log2_fc) >= 2, "high_fc", "other")
    )

    # ggplot2
    xlab = bquote(.(tissue) ~ log[2]~"(median TPM + pseudocount)")
    ylab = bquote("VESTIB" ~ log[2]~"(median TPM + pseudocount)")
    p = ggplot(expr_compare_plot,
            aes(x = log2_gtex, y = log2_vest)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_point(aes(color = high_fc), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("other" = "grey70", "high_fc" = "red")) +
        labs(
            x = paste0(tissue, " log2(median TPM + pseudocount)"),
            y = "VESTIB log2(median TPM + pseudocount)",
            x = xlab,
            y = ylab,
            color = "|log2 FC| ≥ 2",
            title = paste0("VESTIB vs ", tissue, " Median Gene Expression")
        ) +
        theme_classic(); p

    # Let's save plot
    # dir.create("figures/expr_comparison", showWarnings = FALSE)
    ggsave(p, file = paste0("figures/expr_comparison/VESTIB_vs_", tissue, "_expr_comparison.png"),
        width = 9, height = 8, dpi = 300)

}


### Candidate tissue specific genes.
# ENSG00000274012  chr 14:49862551_49862849, Gene = RN7SL2, Type = scRNA, RNA component of signal recognition particle 7SL2
# ENSG00000283293, chr 6:52995620-52995951,  Gene = RN7SK,  Type = snRNA, RNA component of 7SK nuclear ribonucleoprotein
# ENSG00000276168, chr 14:49586580_49586878, Gene = RN7SL1, Type = scRNA, RNA component of signal recognition particle 7SL1
# ENSG00000277209, chr 14:20343071_20343411, Gene = RPPH1, Type = ncRNA,  ribonuclease P RNA component H1


### END

# # Adding gene names
# label_genes <- c(
#     "ENSG00000120885",
#     "ENSG00000100473"
# )

# p = ggplot(expr_compare_plot,
#         aes(x = log2_gtex, y = log2_vest)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_point(aes(color = high_fc), alpha = 0.6, size = 1.5) +
#     scale_color_manual(values = c("other" = "grey70", "high_fc" = "red")) +
#     geom_text(
#         data = expr_compare_plot %>% filter(gene_id %in% label_genes),
#         aes(label = gene_id),
#         size = 3,
#         vjust = -0.7
#     ) +
#     labs(
#         x = bquote(.(tissue_names[i])~log[2]~"(median TPM + pseudocount)"),
#         y = bquote("VESTIB"~log[2]~"(median TPM + pseudocount)"),
#         color = "|log2 FC| ≥ 2",
#         title = paste0("VESTIB vs ", tissue_names[i], " Median Gene Expression")
#     ) +
#     theme_classic(); p

# # Plot
# ggsave(p, file = paste0("figures/expr_comparison/VESTIB_vs_", tissue, "_expr_comparison_2.png"),
#         width = 9, height = 8, dpi = 300)


# # Coloring based on log2_fc
# xlab = bquote(.(tissue) ~ log[2]~"(median TPM + pseudocount)")
# ylab = bquote("VESTIB" ~ log[2]~"(median TPM + pseudocount)")
# p = ggplot(expr_compare_plot,
#         aes(x = log2_gtex, y = log2_vest, color = log2_fc)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     geom_point(alpha = 0.7, size = 1.5) +
#     scale_color_gradient2(
#         low = "blue",
#         mid = "white",
#         high = "red",
#         midpoint = 0,
#         name = "log2 FC"
#     ) +
#     labs(
#         x = xlab,
#         y = ylab,
#         title = paste0("VESTIB vs ", tissue, " Median Gene Expression")
#     ) +
#     theme_classic(); p



###---------------------------------------------- To run t-tests
# Getting protein coding and lncRNAs gene IDs
gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")

# Getting a list of genes that are protein coding or lncRNA
genes_protien_lncRNA <- gtf %>%
    dplyr::filter(type == "gene") %>% 
    dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
    dplyr::pull("gene_id")

# Vestibular tissue, calculating median gene expression
expr_vest_tpm <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression.bed") %>%
    dplyr::select(-`#chr`, -start, -end) %>%
    dplyr::rename(gene_id = phenotype_id) %>%
    filter(gene_id %in% genes_protien_lncRNA)

median_expr_vest <- data.table(
    gene_id = expr_vest_tpm$gene_id,
    VESTIB = log2(
        rowMedians(as.matrix(expr_vest_tpm[, -1])) + pseudocount
    )
)

# Reading GTEx tissue expression
for (i in seq_along(tissue_names)) {

    tissue = tissue_names[i]
    cat("Looking at:", tissue, "Tissue", "\n")
    expr_gtex <- fread(paste0(path_to_gtex_expr, tissue,"/output/unnorm/expression.bed")) %>%
        dplyr::select(-`#chr`, -start, -end) %>%
        dplyr::rename(gene_id = phenotype_id) %>%
        filter(gene_id %in% genes_protien_lncRNA)

    # Convert numeric columns to matrix and compute row medians and log2 transform
    median_expr_gtex <- data.table(
        gene_id = expr_gtex$gene_id,
        median_expr_gtex = log2(
            rowMedians(as.matrix(expr_gtex[, -1])) + pseudocount
        ))

    # Compare median expression between vestibular and GTEx tissue
    col = tissue
    expr_compare <- median_expr_vest %>%
        # left_join(median_expr_gtex, by = "gene_id") %>%
        inner_join(median_expr_gtex, by = "gene_id") %>%
        dplyr::rename(!!col := median_expr_gtex)

    if(i == 1){
        expr_tabal_log2_median <- expr_compare
    } else {
        expr_tabal_log2_median <- expr_tabal_log2_median %>%
        inner_join(expr_compare[, -2], by = "gene_id")
    }
}
dim(expr_tabal_log2_median)
expr_tabal_log2_median[1:5, 1:5]

## Running t-test
dt <- copy(expr_tabal_log2_median)

# Removing genes with zero expression
dt <- dt[rowSums(dt[, -1, with = FALSE] > 0) > 0]

# Filtering genes only expressed in vestibilur tissue
dt[, all_other_zero :=
        rowSums(.SD[, setdiff(names(.SD), "VESTIB"), with = FALSE] > 0) == 0,
    .SDcols = setdiff(names(dt), "gene_id")]

# put them in another df
df_all_other_zero <- dt[all_other_zero == TRUE]

# Save
fwrite(df_all_other_zero,
    file = "figures/data_tables/vestib_vs_gtex_only_expressed_in_vestib_log2.tsv",
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# keeping the other genes and removing "all_other_zero" column
dt_filtered <- dt[all_other_zero == FALSE, -"all_other_zero"]
dim(dt_filtered)
View(dt_filtered[1:5, 1:5])

# Getting tissue names
tissue_cols <- setdiff(names(dt), "gene_id")

# run t-test
res <- dt_filtered[, {
    vestib <- VESTIB
    others <- unlist(.SD[, setdiff(names(.SD), "VESTIB"), with = FALSE],
                    use.names = FALSE)
    others <- others[is.finite(others)]
    # one-sided: VESTIB higher
    p_more <- t.test(others, mu = vestib, alternative = "less")$p.value
    # one-sided: VESTIB lower
    p_less <- t.test(others, mu = vestib, alternative = "greater")$p.value
    # two-sided (optional)
    p_two  <- t.test(others, mu = vestib, alternative = "two.sided")$p.value
    .(
        VESTIB = vestib,
        mean_other = median(others),
        diff = vestib - median(others),
        p_vestib_higher = p_more,
        p_vestib_lower  = p_less,
        p_two_sided     = p_two
    )
}, by = gene_id, .SDcols = setdiff(names(dt_filtered), "gene_id")]

# Mulitiple testing correction
res[, padj_two_sided  := p.adjust(p_two_sided, method = "BH")]

# Order
res[order(padj_two_sided)][1:20]

# Save comparison table
fwrite(res,
    file = "figures/data_tables/vestib_vs_gtex_t_tests.tsv",
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
head(res)

# Plotting t-test results
thr = 0.01/nrow(res)
res$signif <- res$padj_two_sided < thr
ggplot(res, aes(x = mean_other, y = VESTIB, color = signif)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
        values = c("FALSE" = "grey70", "TRUE" = "red"),
        name = "padj < 2.8e-7"
    ) +
    theme_classic() +
    labs(
        x = "Median Expression Across GTEx Tissues (TPM)",
        y = "Median Expression in Vestibular Tissue (TPM)",
        title = "Vestibular vs GTEx Median Gene Expression"
    )


# Another coloring scheme
res$class <- "NS"
res$class[res$padj_two_sided < thr & res$diff > 0] <- "Vestib > GTEx"
res$class[res$padj_two_sided < thr & res$diff < 0] <- "Vestib < GTEx"

p <- ggplot(res, aes(mean_other, VESTIB, color = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.75, size = 1.5) +
    scale_color_manual(
        values = c(
        "Vestib > GTEx" = "chocolate",
        "Vestib < GTEx" = "cadetblue4",
        "NS" = "grey75"
        )
    ) +
    theme_classic() +
    labs(
        x = "Median Expression (log2 TPM) Across GTEx Tissues",
        y = "Median Expression (log2 TPM) in Vestibular Tissue (TPM)",
        color = "P-adj < 2.8e-7",
        title = "Vestibular vs GTEx Median Gene Expression"
    )
p

p <- p +
    theme(
        legend.position = c(0.9, 0.5)   # (x, y) in npc coordinates (0–1)
    )

p <- p +
    theme(
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold")
    )
p

# Save it
ggsave(p, file = "figures/VESTIB_vs_all_GTEx_expr_comparison.png",
        width = 6, height = 6, dpi = 300)



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

# Function for ploting
plot_gene_bar <- function(gene, dt_filtered, gtex_palette) {

    x <- dt_filtered[gene_id == gene]

    long <- melt(
        x,
        id.vars = "gene_id",
        variable.name = "tissue",
        value.name = "expr"
    )

    ggplot(long, aes(x = reorder(tissue, -expr), y = expr, fill = tissue)) +
        geom_col() +
        scale_fill_manual(values = gtex_palette, drop = FALSE) +
        theme_classic() +
        theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
        ) +
        labs(
        title = gene,
        x = NULL,
        y = "log2(median TPM)"
    )
}


# Plot
p = plot_gene_bar("ENSG00000167232", dt_filtered, gtex_palette); p

# Function for ploting with gene symbol
plot_gene_bar <- function(gene, gene_symbol, dt_filtered, gtex_palette) {

    x <- dt_filtered[gene_id == gene]

    long <- melt(
        x,
        id.vars = "gene_id",
        variable.name = "tissue",
        value.name = "expr"
    )

    ggplot(long, aes(x = reorder(tissue, -expr), y = expr, fill = tissue)) +
        geom_col() +
        scale_fill_manual(values = gtex_palette, drop = FALSE) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "none"
        ) +
        labs(
            title = paste0(gene, " (", gene_symbol, ")"),
            x = NULL,
            y = "log2(median TPM)"
        ) +
        theme(
            axis.title = element_text(size = 6),
            axis.text = element_text(size = 4),
            plot.title = element_text(size = 6)
        )
}


p <- plot_gene_bar(
    "ENSG00000167232",
    "ZNF91",
    dt_filtered,
    gtex_palette
)

# Save the plot
ggsave(
    filename = "figures/expr_across_gtex_ENSG00000167232.pdf",
    plot = p,
    units = "mm",
    width = 90,
    height = 80,
    dpi = 300
)


# ### Plot many genes at the same time
# top_genes <- res[order(padj_two_sided)][1:12, gene_id]
# long_top <- melt(
#     dt_filtered[gene_id %in% top_genes],
#     id.vars = "gene_id",
#     variable.name = "tissue",
#     value.name = "expr"
# )
# # Plot
# ggplot(long_top, aes(x = tissue, y = expr, fill = (tissue == "VESTIB"))) +
#     geom_col(show.legend = FALSE) +
#     facet_wrap(~ gene_id, scales = "free_y") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#     labs(x = NULL, y = "log2(median TPM + pseudocount)")

### END