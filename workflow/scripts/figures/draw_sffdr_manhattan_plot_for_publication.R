library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(tomics)

width <- 6
height <- 4
ylim <- c(1, 1e-28)

theme_set(theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 6, angle = 90, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ))

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]

p_col <- snakemake@params$p_col

if(p_col == 'fp') {
  y_label <- "-log<sub>10</sub>(functional p-value)"
} else {
  y_label <- "-log<sub>10</sub>(p-value)"
}

dat <- fread(snakemake@input$gwas)

procd_sumstats <- process_sumstats_for_manhattan(dat, stat_cols = p_col, chr_col = chr_col, bp_col = bp_col)

lead_snps <- fread(snakemake@input$lead_snps)

setnames(lead_snps, c(snakemake@config[["col_names"]][["chromosome"]], snakemake@config[["col_names"]][["base_pair_location"]], 'genes'), c('chr', 'bp', 'gene'))

lead_snps[gene == 'ZNF675', gene := 'ZNF91']

ggsave(
  draw_manhattan(procd_sumstats, lead_snps = lead_snps, stat_col = p_col, y_limits = ylim, y_axis_breaks = 10^(-seq(0, 28, by = 2)), palette = c("gray60", "lightpink3"), repel_args = list(size = 4, angle = 0, colour = "black", min.segment.length = 0, nudge_y = 1)) +
  ylab(y_label) +
  xlab("Chromosome") +
  geom_hline(yintercept = 5e-8, linetype = "dashed", color = "black"),
  file = snakemake@output[[1]],
  width = width,
  height = height
)
