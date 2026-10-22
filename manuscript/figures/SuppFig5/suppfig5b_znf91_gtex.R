## Supplementary Figure 5b: ZNF91 expression across GTEx tissues and the otolith organs (VESTIB).
## Value per tissue = log2(median TPM + 1) over its samples, as in analyses/08_expression_gtex/scripts/expression_comparsion.R
## (all vestibular samples; GTEx tissues run through Pantry, TESTIS left out). Only ZNF91 is read, so this runs in seconds.
## Output: panels/expr_across_gtex_ENSG00000167232.{pdf,png}, source_data/expr_across_gtex_ENSG00000167232.tsv.
## Run from the project root (env r_plots):  Rscript manuscript/figures/SuppFig5/suppfig5b_znf91_gtex.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(jsonlite) })
cfg <- yaml::read_yaml("config/config.yaml")   # + lab-specific values (config_lab.yaml, not in git)
if (file.exists("config/config_lab.yaml")) cfg <- modifyList(cfg, yaml::read_yaml("config/config_lab.yaml"))
O <- "manuscript/figures/SuppFig5"
for (d in c("panels", "source_data")) dir.create(file.path(O, d), recursive = TRUE, showWarnings = FALSE)
GENE <- "ENSG00000167232"; SYMBOL <- "ZNF91"

## log2(median TPM + 1) of one gene in one Pantry expression table (reads only the header and that gene's row)
med_log2 <- function(f) {
  d <- fread(cmd = sprintf("grep -P '^#chr|\\t%s\\t' %s", GENE, shQuote(f)))
  stopifnot(nrow(d) == 1)
  log2(median(as.numeric(unlist(d[, -(1:4)]))) + 1)
}
gtex_root <- cfg$refs$gtex_pantry_phenos
tissues <- setdiff(list.dirs(gtex_root, recursive = FALSE, full.names = FALSE), c("ref", "TESTIS"))
tissues <- tissues[file.exists(file.path(gtex_root, tissues, "output/unnorm/expression.bed"))]
dt <- rbind(data.table(tissue = "VESTIB", expr = med_log2("analyses/05_pantry/phenotyping/output/vestibular/unnorm/expression.bed")),
            rbindlist(lapply(tissues, function(t) data.table(tissue = t, expr = med_log2(file.path(gtex_root, t, "output/unnorm/expression.bed"))))))
dt <- dt[order(-expr)][, tissue := factor(tissue, levels = tissue)]
fwrite(dt, file.path(O, "source_data", "expr_across_gtex_ENSG00000167232.tsv"), sep = "\t")

## GTEx tissue colours (resources/gtex_colors.json), otolith organs in deep sky blue
gc <- fromJSON("resources/gtex_colors.json")
pal <- c(VESTIB = "deepskyblue4", setNames(paste0("#", sapply(gc, `[[`, "tissue_color_hex")), sapply(gc, `[[`, "tissue_abbrv")))

p <- ggplot(dt, aes(tissue, expr, fill = tissue)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sprintf("%s (%s)", GENE, SYMBOL), x = NULL, y = "log2(median TPM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5), axis.title = element_text(size = 6), plot.title = element_text(size = 6))
for (dv in c("pdf", "png"))
  ggsave(file.path(O, "panels", paste0("expr_across_gtex_ENSG00000167232.", dv)), p, units = "mm", width = 150, height = 80,
         dpi = 300, bg = "white", device = if (dv == "pdf") cairo_pdf else ragg::agg_png)
cat(sprintf("%d tissues + VESTIB; VESTIB rank %d\n", length(tissues), which(levels(dt$tissue) == "VESTIB")))
cat("[done]\n")
