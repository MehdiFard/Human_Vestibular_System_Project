## Fig 5c: LocusCompare, GWAS vs QTL -log10 P at the ZNF91 locus (3.6 x 1.8 in, the panel c slot) -> panels/Fig5c.{rds,pdf,png}.
## Also builds the alternative views considered for this slot (TWAS vs PP4, coloc PP0-PP4 bars, cross-ancestry LD;
## the last two are in Supp Fig 4 via suppfig_coloc_ld.R).
## Inputs: TWAS hits (analyses/06_twas), data/cache_znf91_region.rds, data/ld_rs695017_1kgEUR.tsv.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(cowplot); library(data.table); library(dplyr); library(tidyr) })
F5 <- "manuscript/figures/Fig5/data"; OUT <- "manuscript/figures/Fig5/panels"; dir.create(OUT, showWarnings = FALSE)
FONT <- "Nimbus Sans"; WC <- 3.6; HC <- 1.8
th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), legend.key.size = unit(6, "pt"), plot.margin = margin(4, 6, 2, 3))
modcol <- c(expression = "#e41a1c", isoforms = "#377eb8", RNA_editing = "#984210")
modlab <- c(expression = "Expression", isoforms = "Isoform ratio", RNA_editing = "RNA editing")

tw <- fread("analyses/06_twas/results/vestibular/output/twas_hits_annot.tsv")[TRAIT == "balance_meta1"]
tw[, gene := fifelse(is.na(gene_name) | gene_name == "", ENSG, gene_name)]
tw[, feat := fcase(MODALITY == "isoforms", sub(".*__", "", ID),
                   MODALITY == "RNA_editing", paste0("site ", format(as.integer(sub(".*_", "", ID)), big.mark = ",")),
                   default = "expression")]
tw[, lab := paste(gene, feat)]

## ---- A: improved TWAS significance vs colocalization
pA <- ggplot(tw, aes(-log10(TWAS.P), COLOC.PP4)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.02, fill = "grey92") +
  annotate("text", x = 8.2, y = 0.99, label = "PP4 > 0.8", size = 5.5 / .pt, family = FONT, hjust = 0, vjust = 1, colour = "grey40") +
  geom_point(aes(fill = MODALITY, shape = MODALITY), size = 1.8, stroke = 0.25, colour = "black") +
  geom_text_repel(aes(label = lab), size = 5 / .pt, family = FONT, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.15, max.overlaps = Inf, seed = 2) +
  scale_fill_manual(values = modcol, labels = modlab, name = NULL) +
  scale_shape_manual(values = c(expression = 22, isoforms = 21, RNA_editing = 24), labels = modlab, name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.8, 1), expand = c(0, 0)) +
  labs(x = expression(-log[10] * italic(P)[TWAS]), y = "Colocalization PP4") + th +
  theme(legend.position = "right")

## ---- B: coloc posterior probabilities per hit
pp <- tw %>% select(lab, MODALITY, TWAS.P, starts_with("COLOC.PP")) %>%
  pivot_longer(starts_with("COLOC.PP"), names_to = "H", values_to = "PP") %>%
  mutate(H = factor(sub("COLOC.PP", "H", H), levels = paste0("H", 0:4)),     # H4 stacked from 0 (left)
         lab = factor(lab, levels = tw[order(COLOC.PP4)]$lab))
hlab <- c(H4 = "H4 shared variant", H3 = "H3 distinct variants", H2 = "H2 QTL only", H1 = "H1 GWAS only", H0 = "H0 neither")
pB <- ggplot(pp, aes(PP, lab, fill = H)) + geom_col(width = 0.75) +
  geom_vline(xintercept = 0.8, linetype = "dashed", linewidth = 0.25) +
  geom_text(data = tw[, .(lab = factor(lab, levels = levels(pp$lab)), PP4 = COLOC.PP4)],
            aes(x = 1.01, y = lab, label = sprintf("%.2f", PP4)), inherit.aes = FALSE,
            hjust = 0, size = 5.5 / .pt, family = FONT) +
  scale_fill_manual(values = c(H4 = "#b2182b", H3 = "#ef8a62", H2 = "#d1e5f0", H1 = "#67a9cf", H0 = "grey85"),
                    breaks = paste0("H", 4:0), labels = hlab, name = "Coloc hypothesis") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 0.8, 1), limits = c(0, 1.09)) +
  annotate("text", x = 1.01, y = nrow(tw) + 0.75, label = "PP4", hjust = 0, size = 5.5 / .pt, family = FONT, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(x = "Posterior probability", y = NULL) + th + theme(axis.text.y = element_text(size = 5))

## ---- C: LocusCompare (GWAS vs QTL), RNA editing site and isoform
cache <- readRDS(file.path(F5, "cache_znf91_region.rds"))
ld <- fread(file.path(F5, "ld_rs695017_1kgEUR.tsv"))[, .(pos, r2)][, .SD[1], by = pos]
gw <- as.data.table(cache$gwas)[, .(SNP, rsID, pos = BP, gwas = -log10(P))]
## labpos: fixed label positions (data units; label's right edge at lx) so labels never cover variants
lc <- function(m, ph, title, labpos) {
  q <- as.data.table(cache$qtl[[m]])[phenotype_id == ph, .(variant_id, qtl = -log10(pval_nominal))]
  d <- merge(gw, q, by.x = "SNP", by.y = "variant_id")
  d <- merge(d, ld, by = "pos", all.x = TRUE)
  d[, ldc := cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"))]
  lab <- merge(d[rsID %in% labpos$rsID, .(rsID, qtl, gwas)], as.data.table(labpos), by = "rsID")
  ggplot(d[order(r2, na.last = FALSE)], aes(qtl, gwas, fill = ldc)) +
    geom_point(shape = 21, size = 1.1, stroke = 0.15, colour = "grey20") +
    geom_segment(data = lab, aes(x = lx, y = ly, xend = qtl, yend = gwas), inherit.aes = FALSE,
                 linewidth = 0.2, colour = "grey30") +
    geom_label(data = lab, aes(x = lx, y = ly, label = rsID), inherit.aes = FALSE, hjust = 1, size = 5 / .pt,
               family = FONT, fill = scales::alpha("white", 0.85), label.size = 0, label.padding = unit(0.6, "pt")) +
    scale_y_continuous(limits = c(0, 28.5), expand = expansion(mult = c(0.02, 0))) +
    scale_fill_manual(values = c("0-0.2" = "#4062d8", "0.2-0.4" = "#00d0e8", "0.4-0.6" = "#1fbf1f", "0.6-0.8" = "#ffa500",
                                 "0.8-1" = "#e5131f"), na.value = "grey80", name = expression(r^2), drop = FALSE) +
    labs(x = bquote(-log[10] * italic(P) ~ .(title)), y = expression(-log[10] * italic(P)[GWAS])) + th
}
## legend inside the editing plot (its bottom right is empty), so both plots use the full width
## editing: rs695017 above the cluster; rs167154 / rs607168 below it (empty area under the red points)
lab_edit <- data.frame(rsID = c("rs695017", "rs167154", "rs607168"), lx = c(9.6, 7.2, 11.6), ly = c(27.8, 20.8, 20.8))
## isoform: rs695017 on top, rs607168 below it, rs167154 below-left of its point
lab_iso  <- data.frame(rsID = c("rs695017", "rs607168", "rs167154"), lx = c(5.2, 4.7, 6.1), ly = c(27.8, 26.0, 20.6))
pC <- plot_grid(lc("RNA_editing", "ENSG00000167232__19_23338819", "editing QTL", lab_edit) +
                  guides(fill = guide_legend(override.aes = list(size = 1.2), ncol = 2, byrow = FALSE)) +
                  theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
                        legend.title = element_text(size = 5.5), legend.text = element_text(size = 5),
                        legend.key.size = unit(5, "pt"), legend.key.width = unit(4, "pt"),
                        legend.spacing.x = unit(1, "pt"), legend.margin = margin(1, 1, 1, 1)),
                lc("isoforms", "ENSG00000167232__ENST00000300619", "isoform QTL", lab_iso) +
                  theme(axis.title.y = element_blank(), legend.position = "none"),
                nrow = 1, rel_widths = c(1, 1))

## ---- D: cross-ancestry LD of the three lead variants (1000 Genomes superpopulations; computed earlier from
## resources/ldref_1kg_hg38 with the IGSR population table)
ldx <- data.table(pop = rep(c("EUR", "AMR", "EAS", "SAS", "AFR"), each = 3),
                  pair = rep(c("rs695017 – rs607168\n(GWAS – editing QTL)", "rs695017 – rs167154\n(GWAS – isoform QTL)",
                               "rs167154 – rs607168\n(isoform – editing)"), 5),
                  r2 = c(1.00, 0.96, 0.96,  0.98, 0.91, 0.90,  0.97, 0.93, 0.90,  0.94, 0.82, 0.77,  0.97, 0.72, 0.71))
ldx[, pop := factor(pop, levels = c("EUR", "AMR", "EAS", "SAS", "AFR"))]
pD <- ggplot(ldx, aes(pop, r2, colour = pair, group = pair)) + geom_line(linewidth = 0.4) + geom_point(size = 1.4) +
  scale_colour_manual(values = c("rs695017 – rs607168\n(GWAS – editing QTL)" = "#984210",
                                 "rs695017 – rs167154\n(GWAS – isoform QTL)" = "#377eb8",
                                 "rs167154 – rs607168\n(isoform – editing)" = "grey40"),
                      breaks = unique(ldx$pair), name = NULL) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = "1000 Genomes superpopulation", y = expression(LD ~ r^2)) + th +
  theme(legend.text = element_text(size = 5.5), legend.key.height = unit(14, "pt"))

sv <- function(p, n) { ggsave(file.path(OUT, paste0(n, ".pdf")), p, width = WC, height = HC, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(n, ".png")), p, width = WC, height = HC, dpi = 400, bg = "white") }
saveRDS(pC, file.path(OUT, "Fig5c.rds"))   # for the vector composite (fig5_assemble.R)
sv(pC, "Fig5c")   # panel c = LocusCompare; the alternative views above are not saved (B and D are in Supp Fig 4)
cat("[done]\n")
