## ====================================================================
## Fig 5 and Supp Fig 4 (ZNF91 locus), two-column layout, 7.2 in (18.3 cm) wide.
##   Fig 5a:     GWAS, isoform ratio, RNA editing and expression QTL tracks + genes (chr19:23.0-23.8 Mb)
##   Supp Fig 4a: left = GWAS, isoform, RNA editing, expression; right = Alt polyA, Alt TSS, splicing, stability;
##                gene track under each column
## QTL track per modality: the ZNF91 phenotype with the smallest p (else the top phenotype in the region).
## Inputs: data/cache_znf91_region.rds (cache_znf91_region.R), data/ld_rs695017_1kgEUR.tsv (ld_rs695017_1kg_eur.py),
## TWAS hits (analyses/06_twas). Outputs: panels/Fig5a.{rds,pdf,png}, ../SuppFig4/panels/Supp_Fig4*.
## Run from the project root on a compute node:
##   Rscript manuscript/figures/Fig5/fig5a_suppfig4_regional.R
## ====================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(ggplot2); library(ggrepel); library(cowplot)
  library(locuszoomr); library(ensembldb); library(AnnotationHub); library(rtracklayer)
  library(dplyr); library(tidyr); library(purrr)
})
pdf(NULL)   # background device for grob/label measurement (otherwise R writes Rplots.pdf); outputs use their own devices

## Gene track: all genes, unnamed Ensembl genes unlabelled, label widths measured at the real panel width;
## thin axis lines (0.25 pt) on every track.
OUT_F5 <- "manuscript/figures/Fig5/panels"; OUT_S4 <- "manuscript/figures/SuppFig4/panels"
dir.create(OUT_S4, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_F5, recursive = TRUE, showWarnings = FALSE)
FONT   <- "Nimbus Sans"            # Helvetica-metric font available on the cluster
W      <- 7.2                      # full page width (in)

CHR <- 19L; GENE <- "ZNF91"; GENE_ENSG <- "ENSG00000167232"
REGION   <- c(23.0e6, 23.8e6)
QTL_YLIM <- c(0, 15)
LD_FILE  <- "manuscript/figures/Fig5/data/ld_rs695017_1kgEUR.tsv"
TWAS     <- "analyses/06_twas/results/vestibular/output/twas_hits_annot.tsv"

## tagged variants (line to ID); nudge = horizontal label offset in Mb
SNPS <- tibble::tribble(
  ~rsID,       ~pos,      ~panel,        ~nudge, ~nudge_y,
  "rs695017",  23421111L, "gwas",         0.12,  -1,
  "rs167154",  23307401L, "isoforms",    -0.07,   0,
  "rs607168",  23459336L, "RNA_editing",  0.07,   0)

NOM <- "%s/output/vestibular/qtl/nominal/%s/vestibular_%s.cis_qtl_pairs.19.parquet"
mods <- tibble::tribble(
  ~mod,          ~label,               ~color,
  "isoforms",    "Isoform ratio QTLs", "#377eb8",
  "expression",  "Expression QTLs",    "#e41a1c",
  "alt_polyA",   "Alt polyA QTLs",     "#ff7f00",
  "alt_TSS",     "Alt TSS QTLs",       "#984ea3",
  "RNA_editing", "RNA editing QTLs",   "#984210",
  "splicing",    "Splicing QTLs",      "#4daf4a",
  "stability",   "Stability QTLs",     "#fdc11c") %>%
  mutate(parquet = sprintf(NOM, ifelse(mod == "RNA_editing", "analyses/05_pantry/pheast",
                                       "analyses/05_pantry/pheast"), mod, mod))

## ---------------------------------------------------------------- data
## region data from the cache (cache_znf91_region.R; chr19:22.9-23.9 Mb) -> renders in seconds
cache  <- readRDS("manuscript/figures/Fig5/data/cache_znf91_region.rds")
gwas19 <- as.data.table(cache$gwas)
gwas_reg <- gwas19 %>% transmute(SNP, rsID, chr = as.integer(CHR), pos = as.integer(BP), p = as.numeric(P)) %>%
  filter(pos >= REGION[1], pos <= REGION[2], !is.na(p), p > 0)
ensDb <- AnnotationHub(localHub = TRUE)[["AH119325"]]

top_pheno <- function(mod, label) {
  q <- cache$qtl[[mod]] %>% transmute(phenotype_id, pval_nominal, vpos = pos) %>%
    filter(vpos >= REGION[1], vpos <= REGION[2], pval_nominal > 0)
  g <- q %>% filter(grepl(paste0("^", GENE_ENSG), phenotype_id)); u <- if (nrow(g)) g else q
  top <- u %>% slice_min(pval_nominal, n = 1, with_ties = FALSE) %>% pull(phenotype_id)
  cat(sprintf("  [%s] %s\n", label, top))
  u %>% filter(phenotype_id == top) %>% transmute(pos = vpos, p = pval_nominal, modality = label, phenotype_id)
}
qtl <- purrr::pmap_dfr(list(mods$mod, mods$label), top_pheno)

title_for <- function(label, ph) {
  if (label %in% c("Isoform ratio QTLs", "Alt TSS QTLs")) return(sprintf("%s: %s", label, regmatches(ph, regexpr("ENST[0-9]+", ph))))
  if (label == "Splicing QTLs") { s <- strsplit(sub(".*__", "", ph), "_")[[1]]
    return(sprintf("%s: chr%s:%.3f-%.3f Mb", label, s[1], as.integer(s[2]) / 1e6, as.integer(s[3]) / 1e6)) }
  if (label == "RNA editing QTLs") { s <- strsplit(sub(".*__", "", ph), "_")[[1]]
    return(sprintf("%s: ZNF91 site chr%s:%s", label, s[1], format(as.integer(s[2]), big.mark = ","))) }
  sprintf("%s: ZNF91", label)
}
titles <- qtl %>% distinct(modality, phenotype_id) %>% mutate(t = map2_chr(modality, phenotype_id, title_for))
titles <- setNames(titles$t, titles$modality)

## ---------------------------------------------------------------- track panels
th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        plot.margin = margin(3, 8, 4, 3))   # bottom: room for the "0" tick label
xlim_mb <- REGION / 1e6
track <- function(df, title, color, ylim, lab = NULL) {
  p <- ggplot(df, aes(pos / 1e6, -log10(p))) +
    geom_point(size = 0.6, alpha = 0.8, colour = color, stroke = 0) +
    annotate("text", x = xlim_mb[1] + 0.01, y = ylim[2], label = title, hjust = 0, vjust = 1,
             size = 6.5 / .pt, family = FONT) +
    scale_x_continuous(limits = xlim_mb, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0), breaks = scales::pretty_breaks(3), oob = scales::oob_squish) +
    labs(x = NULL, y = expression(-log[10] * italic(P))) + th +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (!is.null(lab) && nrow(lab))
    p <- p + geom_text_repel(data = lab, aes(pos / 1e6, -log10(p), label = rsID), inherit.aes = FALSE,
                             size = 6.5 / .pt, family = FONT, nudge_x = lab$nudge, nudge_y = lab$nudge_y,
                             segment.size = 0.25, segment.colour = "grey30", min.segment.length = 0,
                             box.padding = 0.2, max.overlaps = Inf, seed = 1)
  p
}
lab_for <- function(df, panel, by = "pos") {
  s <- SNPS %>% filter(panel == !!panel)
  if (by == "rsID") df %>% inner_join(s %>% select(rsID, nudge, nudge_y), by = "rsID") %>% distinct(rsID, pos, p, nudge, nudge_y)
  else df %>% inner_join(s %>% select(rsID, pos, nudge, nudge_y), by = "pos") %>% distinct(rsID, pos, p, nudge, nudge_y)
}
gwas_ymax <- ceiling(max(-log10(gwas_reg$p)) / 5) * 5
p_gwas <- track(gwas_reg, "Dizziness GWAS meta-analysis (all ancestries)", "black", c(0, gwas_ymax), lab_for(gwas_reg, "gwas", "rsID"))
qtl_track <- function(m) { i <- match(m, mods$mod); df <- qtl %>% filter(modality == mods$label[i])
  track(df, titles[[mods$label[i]]], mods$color[i], QTL_YLIM, lab_for(df, m)) }
## named protein-coding and lncRNA genes only (drops unnamed ENSG entries and pseudogenes; keeps e.g. LINC01224)
named_genes <- function(loc) {
  tx <- as.data.frame(loc$TX)
  unique(tx$gene_name[tx$gene_biotype %in% c("protein_coding", "lncRNA") & !is.na(tx$gene_name) &
                        tx$gene_name != "" & !grepl("^ENSG", tx$gene_name)])
}
genes <- function(highlight) {
  loc <- locus(data = gwas_reg %>% transmute(rsid = SNP, chr, pos, p), seqname = CHR, xrange = REGION, ens_db = ensDb)
  a <- list(loc, filter_gene_name = named_genes(loc), cex.text = 0.45, maxrows = 4)
  if (highlight) a <- c(a, list(highlight = GENE, highlight_col = "red"))
  do.call(gg_genetracks, a) + coord_cartesian(xlim = xlim_mb, expand = FALSE, clip = "on") +
    xlab(sprintf("Chromosome %d (Mb)", CHR)) +
    theme(text = element_text(family = FONT), axis.text.x = element_text(size = 6, colour = "black"),
          axis.title.x = element_text(size = 7), plot.margin = margin(1, 8, 2, 3),
          axis.line.x = element_line(linewidth = 0.25, colour = "black"),     # same weight as the other tracks
          axis.ticks.x = element_line(linewidth = 0.25, colour = "black"))
}
column <- function(panels, heights) plot_grid(plotlist = cowplot::align_plots(plotlist = panels, align = "v", axis = "lr"),
                                              ncol = 1, rel_heights = heights)

## ---------------------------------------------------------------- Fig 5b (LD from 1000 Genomes EUR)
ld <- fread(LD_FILE)
g5b <- gwas19 %>% transmute(SNP, rsID, chr = as.integer(CHR), pos = as.integer(BP), A1, A2, p = as.numeric(P)) %>%
  filter(!is.na(p), p > 0)
g5b <- g5b %>% left_join(ld %>% select(pos, a1, a2, r2), by = "pos", relationship = "many-to-many") %>%
  filter(is.na(a1) | (A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)) %>%
  group_by(SNP) %>% slice(1) %>% ungroup() %>% select(-a1, -a2) %>% as.data.frame()
loc5b <- locus(data = g5b, gene = GENE, flank = 2.5e5, ens_db = ensDb, labs = "rsID", index_snp = "rs695017")
loc5b$data$ld <- loc5b$data$r2
loc5b <- link_recomb(loc5b, genome = "hg38", recomb = import.bw("resources/ucsc_recomb/hg38/recomb1000GAvg.bw"))
cat(sprintf("Fig 5b: %d SNPs, %d with r2 (%.0f%%)\n", nrow(loc5b$data), sum(!is.na(loc5b$data$ld)), 100 * mean(!is.na(loc5b$data$ld))))
print(loc5b$data %>% filter(rsID %in% SNPS$rsID) %>% select(rsID, pos, p, ld))
xr <- loc5b$xrange / 1e6
lab5b <- loc5b$data %>% filter(rsID %in% SNPS$rsID) %>%
  mutate(x = pos / 1e6, y = -log10(p),
         nx = c(rs167154 = -0.06, rs695017 = 0, rs607168 = 0.06)[rsID],
         ny = c(rs167154 = 2.5, rs695017 = 3.5, rs607168 = 2.5)[rsID])
p5b_sc <- gg_scatter(loc5b, labels = NULL, recomb_col = "#2b6cb0", legend_pos = "topleft", size = 1.1,
                     ylim = c(0, 31), cex.axis = 0.75, cex.lab = 0.85) +
  geom_text_repel(data = lab5b, aes(x, y, label = rsID), inherit.aes = FALSE, size = 6.5 / .pt, family = FONT,
                  nudge_x = lab5b$nx, nudge_y = lab5b$ny, segment.size = 0.25, min.segment.length = 0,
                  box.padding = 0.15, max.overlaps = Inf, seed = 1) +
  coord_cartesian(xlim = xr, expand = FALSE, clip = "off") +
  theme(text = element_text(family = FONT, size = 7), axis.text = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        legend.key.size = unit(7, "pt"), axis.text.x = element_blank(), axis.title.x = element_blank(),
        plot.margin = margin(3, 3, 1, 3))
## all genes: shows how gene-dense this high-LD block is.
## Unnamed Ensembl genes are drawn without a label (name set to " "); all 38 genes need 12 rows (overlapping
## gene bodies). locuszoomr measures label widths on the open device, so open one of the panel's real width
## (otherwise R's default 7-in device makes labels look half as wide and they overlap).
unnamed <- is.na(loc5b$TX$gene_name) | loc5b$TX$gene_name == "" | grepl("^ENSG", loc5b$TX$gene_name)
loc5b$TX$gene_name[unnamed] <- " "
pdf(NULL, width = 3.3, height = 2)
p5b_genes <- gg_genetracks(loc5b, cex.text = 0.38, maxrows = 12,
                           highlight = GENE, highlight_col = "red") +
  coord_cartesian(xlim = xr, expand = FALSE, clip = "off") + xlab(sprintf("Chromosome %d (Mb)", CHR)) +
  theme(text = element_text(family = FONT), axis.text.x = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 7), plot.margin = margin(1, 3, 2, 3))
invisible(dev.off())
p5b <- plot_grid(plotlist = cowplot::align_plots(p5b_sc, p5b_genes, align = "v", axis = "lr"), ncol = 1, rel_heights = c(2.1, 1.7))

## ---------------------------------------------------------------- Fig 5c (TWAS vs coloc, dizziness)
modcol <- c(expression = "#e41a1c", isoforms = "#377eb8", RNA_editing = "#984210", splicing = "#4daf4a",
            alt_TSS = "#984ea3", alt_polyA = "#ff7f00", stability = "#fdc11c")
modlab <- c(expression = "Expression", isoforms = "Isoform ratio", RNA_editing = "RNA editing", splicing = "Intron excision",
            alt_TSS = "Alt. TSS", alt_polyA = "Alt. polyA", stability = "RNA stability")
tw <- fread(TWAS) %>% filter(TRAIT == "balance_meta1") %>%
  mutate(lab = case_when(MODALITY == "isoforms" & grepl("ENST00000300619", ID) ~ "ZNF91 isoform ENST00000300619",
                         MODALITY == "RNA_editing" & grepl("19_23338819", ID) ~ "ZNF91 RNA editing, chr19:23,338,819",
                         TRUE ~ NA_character_),
         gene = ifelse(is.na(gene_name) | gene_name == "", ENSG, gene_name))
p5c <- ggplot(tw, aes(-log10(TWAS.P), COLOC.PP4, colour = MODALITY)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.25, colour = "grey50") +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_text_repel(aes(label = lab), size = 6 / .pt, family = FONT, colour = "black", min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.3, nudge_x = -7, nudge_y = 0.1, direction = "y",
                  hjust = 1, na.rm = TRUE, seed = 1) +
  scale_colour_manual(values = modcol, labels = modlab, name = NULL) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 0.8, 1)) +
  labs(x = expression(-log[10] * italic(P)[TWAS]), y = "Colocalization PP4") + th +
  theme(legend.position = c(0.02, 0.45), legend.justification = c(0, 0.5), legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 6), legend.background = element_blank(), plot.margin = margin(6, 6, 2, 3))

## ---------------------------------------------------------------- assemble Fig 5
f5_tracks <- list(p_gwas, qtl_track("isoforms"), qtl_track("RNA_editing"), qtl_track("expression"), genes(TRUE))
col_a <- column(f5_tracks, c(1.25, 1, 1, 1, 0.95))
save_fig <- function(p, stem, w, h, dir) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
}
H5 <- 5.6
right <- plot_grid(p5b, p5c, ncol = 1, rel_heights = c(3.8, 1.8), labels = c("b", "c"),
                   label_size = 10, label_fontfamily = FONT)
fig5 <- plot_grid(col_a, right, nrow = 1, rel_widths = c(1, 1), labels = c("a", ""), label_size = 10, label_fontfamily = FONT)
# (earlier all-ggplot Fig 5 composites removed; the paper figure is made by fig5_assemble.R)
fig5_noc <- plot_grid(col_a, plot_grid(p5b, NULL, ncol = 1, rel_heights = c(3.8, 1.8)), nrow = 1,
                      labels = c("a", "b"), label_size = 10, label_fontfamily = FONT)

saveRDS(col_a, file.path(OUT_F5, "Fig5a.rds"))   # for the vector composite (fig5_assemble.R)
save_fig(col_a, "Fig5a", W / 2, H5, OUT_F5)

## ---------------------------------------------------------------- Supp Fig 4 (two columns)
left  <- column(list(p_gwas, qtl_track("isoforms"), qtl_track("RNA_editing"), qtl_track("expression"), genes(FALSE)),
                c(1.25, 1, 1, 1, 0.95))
rightc <- column(list(qtl_track("alt_polyA"), qtl_track("alt_TSS"), qtl_track("splicing"), qtl_track("stability"), genes(FALSE)),
                 c(1.25, 1, 1, 1, 0.95))
s4 <- plot_grid(left, rightc, nrow = 1)
saveRDS(s4, file.path(OUT_S4, "Supp_Fig4.rds"))   # for the merged Supp Fig 4 + coloc/LD (suppfig_coloc_ld.R)
save_fig(s4, "Supp_Fig4_regional_plots", W, H5, OUT_S4)
cat("[done]\n")


### END