## Supp Fig 3: multimodal xQTL mapping in the otolith organs.
##   a  protein-coding genes by number of phenotypes per gene, per modality   (Pantry paper figure1.R, panel C)
##   b  genes by total number of xQTLs                                         (figure2abc.R, panel B)
##   c  modality of the xQTLs by rank within the gene                          (figure2abc.R, panel C)
##   d  xQTL position along the gene, one row per modality, with +-50 kb TSS / TES insets   (figure2d.R)
##   e  functional annotation enrichment of the xVariants (mean +- SD over tissues) and frequency   (figure2e.R)
## Inputs: analyses/05_pantry/results/xqtl_summary (phenotypes, xQTLs, xQTL positions, enrichment).
## Writes panels/ (one PDF per panel), source_data/ (data behind each panel) and SuppFig3.pdf / SuppFig3.png.
## Run from the project root:  Rscript manuscript/figures/SuppFig3/suppfig3.R
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
X <- "analyses/05_pantry/results/xqtl_summary"
O <- "manuscript/figures/SuppFig3"
for (d in c("panels", "source_data")) dir.create(file.path(O, d), recursive = TRUE, showWarnings = FALSE)

modalities <- c(expression = "Expression", isoforms = "Isoform ratio", splicing = "Intron excision",
                alt_TSS = "Alt. TSS", alt_polyA = "Alt. polyA", stability = "RNA stability", RNA_editing = "RNA editing")
modality_colors <- c(Expression = "#e41a1c", `Isoform ratio` = "#377eb8", `Intron excision` = "#4daf4a",
                     `Alt. TSS` = "#984ea3", `Alt. polyA` = "#ff7f00", `RNA stability` = "#fdc11c", `RNA editing` = "#984210")
save_panel <- function(p, stem, w, h) ggsave(file.path(O, "panels", paste0(stem, ".pdf")), p, width = w, height = h)
save_data <- function(d, stem) write_tsv(d, file.path(O, "source_data", paste0(stem, ".tsv")))

## ---- a: phenotypes per protein-coding gene ------------------------------------------------------------------------
data_a <- read_tsv(file.path(X, "vestibular.phenotypes.tsv.gz"), col_types = "cccc") |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(gene_biotype == "protein_coding") |>
    dplyr::count(modality, gene_id, name = "n_phenos") |>
    mutate(n_phenos = pmin(n_phenos, 10), modality = fct_rev(modality)) |>
    group_by(modality, n_phenos) |>
    summarise(thousands = n() / 1000, .groups = "drop")
p_a <- ggplot(data_a, aes(x = thousands, y = modality, fill = n_phenos)) +
    geom_col(width = 0.8) +
    scale_fill_viridis_c(breaks = c(1, 4, 7, 10), labels = c(1, 4, 7, "10+")) +
    xlab("Protein-coding genes (×1000)") +
    ylab("Modality") +
    labs(fill = "Phenotypes\nper gene") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.title.y = element_text(margin = margin(r = 10)),   # keep "Modality" clear of the tick labels
          legend.key.size = unit(9, "pt"),
          legend.text = element_text(size = 8),
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 8))
save_panel(p_a, "SuppFig3a_phenotypes_per_gene", 4.16, 2.08)
save_data(data_a, "SuppFig3a_phenotypes_per_gene")

## ---- b: xQTLs per gene ------------------------------------------------------------------------------------------
qtls <- read_tsv(file.path(X, "vestibular.comb.qtls.tsv"), col_types = "cicccd") |>
    mutate(modality = factor(modalities[modality], levels = modalities),
           rank = fct_other(as.character(rank), keep = as.character(1:9), other_level = "10+"))
data_b <- qtls |>
    dplyr::count(gene_id) |>
    mutate(n_qtls = pmin(n, 5) |> as.character() |> str_replace("5", "5+")) |>
    dplyr::count(n_qtls, name = "genes")
p_b <- ggplot(data_b, aes(x = n_qtls, y = genes)) +
    geom_col(position = "dodge", width = 0.5, color = "black", fill = NA, linewidth = 0.5) +
    geom_text(aes(label = genes), vjust = -0.4, size = 2.2) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    xlab("Total xQTLs per gene") +
    ylab("No. of Genes") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
save_panel(p_b, "SuppFig3b_xqtls_per_gene", 2.2, 2.3)
save_data(data_b, "SuppFig3b_xqtls_per_gene")

## ---- c: modality by xQTL rank -----------------------------------------------------------------------------------
rank_counts <- qtls |> dplyr::count(rank) |> slice(1:5)
p_c <- qtls |>
    filter(rank %in% 1:5) |>
    ggplot(aes(x = rank, fill = modality)) +
    geom_bar(position = "fill", width = 0.8) +
    geom_text(aes(x = rank, y = 1.27, label = n, fill = NULL), data = rank_counts, hjust = 1, angle = 90, size = 3) +
    annotate("text", x = -0.2, y = 1.2, label = "n =", size = 3) +
    scale_fill_manual(values = modality_colors) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(1, 5), ylim = c(0, 1), clip = "off") +
    labs(fill = "Modality") +
    xlab("Rank of xQTL within gene") +
    ylab("Proportion of xQTLs") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.key.size = unit(10, "pt"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.margin = NULL,
          plot.margin = margin(30, 5.5, 5.5, 5.5))
save_panel(p_c, "SuppFig3c_xqtl_rank_by_modality", 3.5, 2.3)
save_data(qtls |> dplyr::count(rank, modality) |> mutate(proportion = n / sum(n), .by = rank),
          "SuppFig3c_xqtl_rank_by_modality")

## ---- d: xQTL position along the gene (100 bins over [-1, 2]) with +-50 kb insets at the TSS and TES --------------
pos <- read_tsv(file.path(X, "vestibular.comb.qtls.rel_pos.tsv"), col_types = "ccddd") |>
    mutate(modality = factor(modalities[modality], levels = modalities))
main_row <- function(m, last) {
    d <- pos |> filter(modality == m, rel_pos_gene >= -1, rel_pos_gene <= 2)
    p <- ggplot(d, aes(rel_pos_gene)) +
        geom_histogram(binwidth = 0.03, boundary = -1, fill = modality_colors[[m]]) +
        geom_vline(xintercept = c(0, 1), alpha = 0.5) +
        annotate("text", x = -0.95, y = Inf, label = m, hjust = 0, vjust = 1.2, size = 2.5) +
        scale_x_continuous(limits = c(-1, 2), expand = c(0, 0), breaks = c(0, 1), labels = c("Gene start", "Gene end")) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 2), expand = expansion(mult = c(0, 0.05))) +
        theme_classic() + theme(axis.text = element_text(colour = "black"), plot.margin = margin(1, 4, 1, 2)) +
        labs(x = if (last) "xVariant position normalized to xGene length" else NULL, y = NULL)
    if (!last) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_line())
    p
}
inset <- function(m, var, mid) {   # 50 bins over +-50 kb, no y axis
    d <- pos |> filter(modality == m, abs(.data[[var]]) <= 5e4)
    ggplot(d, aes(.data[[var]])) +
        geom_histogram(bins = 50, fill = modality_colors[[m]]) +
        scale_x_continuous(limits = c(-5e4, 5e4), expand = c(0, 0), breaks = c(-5e4, 0, 5e4),
                           labels = c("-50 Kb", mid, "+50 Kb")) +
        theme_classic() +
        theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_text(colour = "black", size = 5.5), axis.title = element_blank(),
              axis.line.x = element_line(linewidth = 0.3), axis.ticks.x = element_line(linewidth = 0.3),
              plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0))
}
## inset placement within each row (fractions of the row): left inset before "Gene start", right after "Gene end"
IN_W <- 0.14; IN_H <- 0.46; IN_Y <- 0.44; TSS_X <- 0.11; TES_X <- 0.73
mods <- unname(modalities)
main <- lapply(seq_along(mods), function(i) main_row(mods[i], i == length(mods)))
save_data(map2_dfr(mods, main, \(m, p) ggplot_build(p)$data[[1]] |>
                       transmute(modality = m, xmin = sprintf("%g", xmin), xmax = sprintf("%g", xmax), n_qtls = count)),
          "SuppFig3d_xqtl_position_along_gene")
rows <- cowplot::align_plots(plotlist = main, align = "v", axis = "lr")
rows <- lapply(seq_along(mods), function(i) {
    h_axis <- if (i == length(mods)) 0.32 else 0          # bottom row also holds the x-axis labels
    ggdraw(rows[[i]]) +
        draw_plot(inset(mods[i], "rel_pos_TSS", "Start"), x = TSS_X, y = h_axis + IN_Y * (1 - h_axis),
                  width = IN_W, height = IN_H * (1 - h_axis)) +
        draw_plot(inset(mods[i], "rel_pos_TES", "End"), x = TES_X, y = h_axis + IN_Y * (1 - h_axis),
                  width = IN_W, height = IN_H * (1 - h_axis))
})
d_stack <- plot_grid(plotlist = rows, ncol = 1, rel_heights = c(rep(1, length(mods) - 1), 1.45))
p_d <- plot_grid(ggdraw() + draw_label("No. xQTLs", angle = 90, size = 11), d_stack, nrow = 1, rel_widths = c(0.035, 1))
save_panel(p_d, "SuppFig3d_xqtl_position_along_gene", 6.6, 6.8)

## ---- e: functional annotation enrichment --------------------------------------------------------------------------
combine_categories <- function(df, grouping) {
    df |>
        mutate(category = if_else(category %in% names(grouping), grouping[category], category)) |>
        summarise(count_qtls = sum(count_qtls), total_qtls = unique(total_qtls),
                  count_bg = sum(count_bg), total_bg = unique(total_bg),
                  .by = c(tissue, modality, category)) |>
        mutate(frac_qtls = count_qtls / total_qtls, frac_bg = count_bg / total_bg,
               log2_enrich = log2(frac_qtls / frac_bg))
}
add_pseudocounts <- function(df, amount) {
    df |>
        dplyr::rename(count_qtls_raw = count_qtls, total_qtls_raw = total_qtls,
                      frac_qtls_raw = frac_qtls, log2_enrich_raw = log2_enrich) |>
        mutate(count_qtls = count_qtls_raw + amount,
               total_qtls = total_qtls_raw + amount / frac_bg,
               frac_qtls = count_qtls / total_qtls,
               log2_enrich = log2(frac_qtls / frac_bg))
}
categories <- c(
    enhancer_d = "Enhancer", promoter_d = "Promoter", open_chromatin_region_d = "Open chromatin",
    promoter_flanking_region_d = "Promoter-flanking", CTCF_binding_site_d = "CTCF binding site",
    TF_binding_site_d = "TF binding site", `3_prime_UTR_variant_d` = "3' UTR", `5_prime_UTR_variant_d` = "5' UTR",
    frameshift_variant_d = "Frameshift", intron_variant_d = "Intron", missense_variant_d = "Missense",
    non_coding_transcript_exon_variant_d = "NC transcript", splice_acceptor_variant_d = "Splice acceptor",
    splice_donor_variant_d = "Splice donor", splice_region_variant_d = "Splice region", stop_gained_d = "Stop gained",
    synonymous_variant_d = "Synonymous", splicing = "Splicing", truncating = "Truncating")

enrich_all <- read_tsv(file.path(X, "enrich.vestibular.comb.tsv"), col_types = "ccciiiiddd") |>
    combine_categories(c(frameshift_variant_d = "truncating", stop_gained_d = "truncating",
                         splice_acceptor_variant_d = "splicing", splice_donor_variant_d = "splicing",
                         splice_region_variant_d = "splicing")) |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(count_qtls >= 2) |>
    add_pseudocounts(0.5)
enrich <- enrich_all |>
    summarise(log2_enrich_mean = mean(log2_enrich), log2_enrich_sd = sd(log2_enrich),
              frac_qtls_mean = mean(frac_qtls), .by = c(modality, category)) |>
    mutate(category = factor(categories[category]) |> fct_reorder(log2_enrich_mean, var))
stripes <- tibble(y = seq(1, length(levels(enrich$category)), by = 2) - 0.5)

p_e1 <- enrich |>
    mutate(category = as.integer(category) + (7/20) - (1/10) * as.integer(modality)) |>
    ggplot(aes(x = log2_enrich_mean, xmin = log2_enrich_mean - log2_enrich_sd, xmax = log2_enrich_mean + log2_enrich_sd,
               y = category, color = modality, shape = modality)) +
    geom_rect(aes(ymin = y, ymax = y + 1, xmin = -Inf, xmax = Inf, color = NULL, shape = NULL, x = NULL, y = NULL),
              fill = "#eeeeee", data = stripes, show.legend = FALSE) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_linerange(linewidth = 0.6, show.legend = FALSE) +
    geom_point(size = 1.5, stroke = 0.75) +
    scale_color_manual(values = modality_colors) +
    scale_shape_manual(values = c(16, 17, 15, 4, 8, 5, 3)) +   # 7 modalities
    expand_limits(x = 5.9) +
    scale_y_continuous(breaks = 1:length(levels(enrich$category)), labels = levels(enrich$category), expand = c(0, 0)) +
    coord_cartesian(ylim = c(1 - 0.6, length(levels(enrich$category)) + 0.6)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.background = element_rect(color = "black", linewidth = 0.25),
          legend.justification = c(1, 0),
          legend.key.height = unit(7, "pt"),
          legend.key.width = unit(7, "pt"),
          legend.margin = margin(2, 3, 2, 3),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.05),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          plot.margin = unit(c(5.5, 1, 5.5, 5.5), "pt")) +
    xlab(expression(log[2]*" fold enrichment in xVariants    ")) +
    ylab(NULL) +
    labs(color = "Modality", shape = "Modality")
p_e2 <- enrich |>
    mutate(modality = fct_rev(modality)) |>
    ggplot(aes(y = category, x = frac_qtls_mean, fill = modality)) +
    geom_col(position = "dodge", width = 0.7, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6)) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.text = element_text(color = "black"), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), legend.position = "none", plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")) +
    ylab(NULL) +
    xlab("Freq. in xVariants           ")   # label moved left since the axis is truncated in the figure
p_e <- plot_grid(p_e1, p_e2, align = "h", axis = "tb", rel_widths = c(0.5, 0.4))
save_panel(p_e, "SuppFig3e_functional_enrichment", 6.5, 6.5)
save_data(enrich_all |>
              dplyr::select(category, modality, tissue, log2_enrichment = log2_enrich, freq_in_xvariants = frac_qtls) |>
              mutate(category = categories[category],
                     log2_enrichment = sprintf("%g", log2_enrichment),
                     freq_in_xvariants = sprintf("%g", freq_in_xvariants)),
          "SuppFig3e_functional_enrichment")

## ---- Supp Fig 3: a, b, c on top; d and e below (canvas at the panels' native sizes; scale to the page width) ------
top <- plot_grid(p_a, p_b, p_c + theme(plot.margin = margin(40, 30, 5.5, 25)),   # room for "n =", rank counts, legend
                 nrow = 1, labels = c("a", "b", "c"), label_size = 14, rel_widths = c(4.3, 2.4, 4.0))
bottom <- plot_grid(p_d, p_e, nrow = 1, labels = c("d", "e"), label_size = 14, rel_widths = c(6.6, 6.2))
fig <- plot_grid(top, bottom, ncol = 1, rel_heights = c(2.5, 6.8))
W <- 13.2; H <- 9.3
ggsave(file.path(O, "SuppFig3.pdf"), fig, width = W, height = H, device = cairo_pdf)
ggsave(file.path(O, "SuppFig3.png"), fig, width = W, height = H, dpi = 300, bg = "white", device = ragg::agg_png)
cat("[done]\n")


### END