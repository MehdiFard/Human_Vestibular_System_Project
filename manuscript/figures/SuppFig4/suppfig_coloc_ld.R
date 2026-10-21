## Supp Fig 4: regional plots of the ZNF91 locus (a) + colocalization posteriors of the dizziness xTWAS hits (b) +
## cross-ancestry LD of the lead variants (c). Outputs:
##   panels/SuppFig4bc_coloc_LD.{pdf,png}   b + c, 7.2 x 2.7 in
##   SuppFig4_5.5x5.5in.{pdf,png}           a + b + c, laid out at 7.2 x 7.2 in and scaled to 5.5 in with ghostscript
## Inputs: TWAS hits (balance_meta1, analyses/06_twas), data/ld_cross_ancestry_3snps.tsv (ld_cross_ancestry.py),
## panels/regional/Supp_Fig4.rds (../Fig5/fig5a_suppfig4_regional.R). Run from the project root (env locuszoomr).
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(cowplot); library(data.table); library(dplyr); library(tidyr) })   # ggrepel: draws the SNP tags of the saved Supp Fig 4
OUT <- "manuscript/figures/SuppFig4"; PAN <- "manuscript/figures/SuppFig4/panels"
FONT <- "Nimbus Sans"; W <- 7.2
th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), legend.key.size = unit(6, "pt"), plot.margin = margin(4, 6, 2, 4))

## ---- a: coloc posterior probabilities (H0-H4) for every dizziness xTWAS hit
tw <- fread("analyses/06_twas/results/vestibular/output/twas_hits_annot.tsv")[TRAIT == "balance_meta1"]
tw[, gene := fifelse(is.na(gene_name) | gene_name == "", ENSG, gene_name)]
tw[, lab := fcase(MODALITY == "isoforms",    sprintf("Isoform %s (%s)", sub(".*__", "", ID), gene),
                  MODALITY == "RNA_editing", sprintf("Editing chr19:%s (%s)", format(as.integer(sub(".*_", "", ID)), big.mark = ","), gene),
                  default = sprintf("Expression (%s)", gene))]
pp <- tw %>% select(lab, starts_with("COLOC.PP")) %>%
  pivot_longer(starts_with("COLOC.PP"), names_to = "H", values_to = "PP") %>%
  mutate(H = factor(sub("COLOC.PP", "H", H), levels = paste0("H", 0:4)),     # H4 stacked from 0 (left)
         lab = factor(lab, levels = tw[order(COLOC.PP4)]$lab))
hlab <- c(H4 = "H4: shared variant", H3 = "H3: distinct variants", H2 = "H2: QTL only", H1 = "H1: GWAS only", H0 = "H0: neither")
pA <- ggplot(pp, aes(PP, lab, fill = H)) + geom_col(width = 0.72) +
  geom_vline(xintercept = 0.8, linetype = "dashed", linewidth = 0.25) +
  geom_text(data = tw[, .(lab = factor(lab, levels = levels(pp$lab)), PP4 = COLOC.PP4)],
            aes(x = 1.015, y = lab, label = sprintf("%.2f", PP4)), inherit.aes = FALSE, hjust = 0, size = 5.5 / .pt, family = FONT) +
  annotate("text", x = 1.015, y = nrow(tw) + 0.85, label = "PP4", hjust = 0, size = 5.5 / .pt, family = FONT, fontface = "bold") +
  scale_fill_manual(values = c(H4 = "#b2182b", H3 = "#ef8a62", H2 = "#d1e5f0", H1 = "#67a9cf", H0 = "grey85"),
                    breaks = paste0("H", 4:0), labels = hlab, name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.8, 1), labels = c("0", "0.25", "0.5", "0.8", "1"),
                     limits = c(0, 1.1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Colocalization posterior probability", y = NULL) + th +
  theme(axis.text.y = element_text(size = 5.6), legend.position = "bottom", legend.text = element_text(size = 5.5),
        legend.key.size = unit(5.5, "pt"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

## ---- b: cross-ancestry LD of the three lead variants (computed from 1000 Genomes, ld_cross_ancestry.py)
ld <- fread("manuscript/figures/SuppFig4/data/ld_cross_ancestry_3snps.tsv")
pair_lab <- c("rs695017|rs607168" = "rs695017 – rs607168 (GWAS – editing QTL)",
              "rs695017|rs167154" = "rs695017 – rs167154 (GWAS – isoform QTL)",
              "rs167154|rs607168" = "rs167154 – rs607168 (isoform – editing QTL)")
ld[, pair := factor(pair_lab[paste(snp1, snp2, sep = "|")], levels = pair_lab)]
pops <- c("EUR", "AMR", "EAS", "SAS", "AFR")
nlab <- ld[, .(n = n[1]), by = pop][match(pops, pop), sprintf("%s\n(n = %d)", pop, n)]
ld[, pop := factor(pop, levels = pops)]
pB <- ggplot(ld, aes(pop, r2, colour = pair, group = pair)) + geom_line(linewidth = 0.4) + geom_point(size = 1.4) +
  scale_colour_manual(values = setNames(c("#984210", "#377eb8", "grey40"), pair_lab), name = NULL) +
  scale_x_discrete(labels = nlab) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = "1000 Genomes superpopulation", y = expression("LD " * r^2)) + th +
  theme(axis.text.x = element_text(size = 5.6, lineheight = 0.9), legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(size = 5.5), legend.key.height = unit(6, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"))

## ---- standalone and merged with Supp Fig 4
LBL <- list(label_size = 10, label_fontfamily = FONT)
row_bc <- function(labs) do.call(plot_grid, c(list(pA, pB, nrow = 1, rel_widths = c(1.45, 1), labels = labs, align = "h", axis = "b"), LBL))
sv <- function(p, stem, w, h, dir = OUT) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
}
sv(row_bc(c("a", "b")), "SuppFig4bc_coloc_LD", W, 2.7, dir = PAN)   # b + c as a standalone panel
s4 <- readRDS("manuscript/figures/SuppFig4/panels/regional/Supp_Fig4.rds")
## square version for Word (5.5 x 5.5 in): laid out at 7.2 x 7.2 in (shorter regional tracks), scaled to 5.5 in with
## ghostscript (vector) + 600 dpi PNG with the DPI tag
merged_sq <- plot_grid(do.call(plot_grid, c(list(s4, labels = "a"), LBL)), row_bc(c("b", "c")), ncol = 1, rel_heights = c(4.75, 2.45))
src <- tempfile(fileext = ".pdf")          # intermediate
ggsave(src, merged_sq, width = 7.2, height = 7.2, device = cairo_pdf)
dst <- file.path(OUT, "SuppFig4_5.5x5.5in")
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=pdfwrite", "-dFIXEDMEDIA", "-dPDFFitPage", "-dAutoRotatePages=/None",
                "-dEmbedAllFonts=true", "-dDEVICEWIDTHPOINTS=396", "-dDEVICEHEIGHTPOINTS=396", paste0("-sOutputFile=", dst, ".pdf"), src))
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=png16m", "-r600", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                paste0("-sOutputFile=", dst, ".png"), paste0(dst, ".pdf")))
print(tw[order(-COLOC.PP4), .(lab, TWAS.P, PP3 = COLOC.PP3, PP4 = COLOC.PP4)])
cat("[done]\n")
