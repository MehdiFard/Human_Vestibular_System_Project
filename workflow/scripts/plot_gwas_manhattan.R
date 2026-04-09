library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(tomics)

if (is.null(snakemake@params[["width"]])) {
  width <- 8
} else {
  width <- snakemake@params[["width"]]
}

if (is.null(snakemake@params[["height"]])) {
  height <- 2
} else {
  height <- snakemake@params[["height"]]
}

if (is.null(snakemake@params$ylim)) {
  ylim <- c(1, 1e-10)
} else {
  ylim <- snakemake@params$ylim
}

theme_set(theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 6, angle = 90, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_blank()
  ))

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]

if (!is.null(snakemake@params$p_col)) {
  p_col <- snakemake@params$p_col
} else {
  p_col <- snakemake@config[["col_names"]][["p_value"]]
}

dat <- fread(snakemake@input$gwas)

procd_sumstats <- process_sumstats_for_manhattan(dat, stat_cols = p_col, chr_col = chr_col, bp_col = bp_col)

if(!is.null(snakemake@input$lead_snps)) {
  lead_snps <- fread(snakemake@input$lead_snps)

  setnames(lead_snps, c(snakemake@config[["col_names"]][["chromosome"]], snakemake@config[["col_names"]][["base_pair_location"]], 'genes'), c('chr', 'bp', 'gene'))

  ggsave(
    draw_manhattan(procd_sumstats, lead_snps = lead_snps, stat_col = p_col, y_limits = c(1, 1e-20), repel_args = list(size = 4, angle = 0, colour = "black", min.segment.length = 0)) +
    geom_hline(yintercept = 5e-8, linetype = "dashed", color = "black"),
    file = snakemake@output[[1]],
    width = width,
    height = height
  )
} else {
  ggsave(
    draw_manhattan(procd_sumstats, stat_col = p_col, y_limits = c(1, 1e-20), repel_args = list(size = 4, angle = 0, colour = "black", min.segment.length = 0)) +
    geom_hline(yintercept = 5e-8, linetype = "dashed", color = "black"),
    file = snakemake@output[[1]],
    width = width,
    height = height
  )
}
