### Comparing gene expression in Vestibular tissue with GTEx tissue to see if there is a tissue specific gene
rm(list=ls())

# Load libraries
library(data.table)
library(dplyr)
library(matrixStats)
library(ggplot2)


# Vestibular tissue, calculating median gene expression
expr_vest_tpm <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression.bed") %>%
    select(-`#chr`, -start, -end) %>%
    dplyr::rename(gene_id = phenotype_id)


# Convert numeric columns to matrix and compute row medians
median_expr_vest <- data.table(
    gene_id = expr_vest_tpm$gene_id,
    median_expr_vest = rowMedians(as.matrix(expr_vest_tpm[, -1]))
)

# GTEx tissues median gene expression
path_to_gtex_expr <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/phenos/"
tissue_names <- setdiff(list.dirs(path_to_gtex_expr, recursive = FALSE, full.names = FALSE), "ref")

# Looping over vestibular tissues
for (i in seq_along(tissue_names)) {
    tissue = tissue_names[i]
    cat("Looking at:", tissue, "Tissue", "\n")
    expr_gtex <- fread(paste0(path_to_gtex_expr, tissue,"/output/unnorm/expression.bed")) %>%
        select(-`#chr`, -start, -end) %>%
        dplyr::rename(gene_id = phenotype_id)

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
        fold_change = (median_expr_vest + 1) / (median_expr_gtex + 1),
        log2_fc = log2((median_expr_vest + 1) / (median_expr_gtex + 1))) %>%
        filter(log2_fc > 1 | log2_fc < -1)  %>%
        dplyr::rename(!!col := median_expr_gtex)

    # Save comparison table
    fwrite(expr_compare,
        file = paste0("figures/expr_comparison/VESTIB_vs_", tissue, "_expr_comparison.tsv"),
        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    # Plot
    expr_compare_plot <- expr_compare %>%
    mutate(
        log2_vest = log2(median_expr_vest + 1),
        log2_gtex = log2(.data[[col]] + 1),
        high_fc = ifelse(abs(log2_fc) >= 2, "high_fc", "other")
    )

    # ggplot2
    xlab = bquote(.(tissue) ~ log[2]~"(median TPM + 1)")
    ylab = bquote("VESTIB" ~ log[2]~"(median TPM + 1)")
    p = ggplot(expr_compare_plot,
            aes(x = log2_gtex, y = log2_vest)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_point(aes(color = high_fc), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("other" = "grey70", "high_fc" = "red")) +
        labs(
            x = paste0(tissue, " log2(median TPM + 1"),
            y = "VESTIB log2(median TPM + 1)",
            x = xlab,
            y = ylab,
            color = "|log2 FC| ≥ 2",
            title = paste0("VESTIB vs ", tissue, " Median Gene Expression")
        ) +
        theme_classic(); p

    # Let's save it
    # dir.create("figures/expr_comparison", showWarnings = FALSE)
    ggsave(p, file = paste0("figures/expr_comparison/VESTIB_vs_", tissue, "_expr_comparison_2.png"),
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
#     "ENSG00000274012",
#     "ENSG00000283293",
#     "ENSG00000276168",
#     "ENSG00000277209"
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
#         x = bquote(.(tissue_names[i])~log[2]~"(median TPM + 1)"),
#         y = bquote("VESTIB"~log[2]~"(median TPM + 1)"),
#         color = "|log2 FC| ≥ 2",
#         title = paste0("VESTIB vs ", tissue_names[i], " Median Gene Expression")
#     ) +
#     theme_classic(); p


# # Coloring based on log2_fc
# xlab = bquote(.(tissue) ~ log[2]~"(median TPM + 1)")
# ylab = bquote("VESTIB" ~ log[2]~"(median TPM + 1)")
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
