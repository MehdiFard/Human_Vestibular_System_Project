## ====================================================================
## Supplementary Figure 3 — regional QTL plots around ZNF91 (chr19)
##
## Two fixes vs scripts/locuszoom_plot.R:
##   (a) Don't slice_min(p) per SNP — that biases the p-value
##       distribution and lifts -log10(p) off 0. Instead, pick ONE
##       phenotype per modality (prefer ZNF91 if available so the cis
##       window covers the locus; fall back to the strongest signal in
##       the region otherwise) and plot ALL of its SNPs.
##   (b) Lock x-axis across panels, fixed ylim = (0, 8) for QTL
##       panels, titles inside the panel area (top-left, indented),
##       enough vertical spacing between panels so y-axis labels at
##       the panel boundaries do not collide.
##
## Output:
##   figures/paper/Supp_Fig3_regional_plots.{pdf,png}        composite
##   figures/paper/Supp_Fig3_panel_<modality>.{pdf,png}      individual
## ====================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(locuszoomr)
  library(ensembldb)
  library(AnnotationHub)
  ## tidyverse pieces LAST so dplyr/tidyr/purrr verbs win
  library(dplyr)
  library(tidyr)
  library(purrr)
})

setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry")

OUT_DIR <- "figures/paper"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ----- Region around ZNF91 ----------------------------------------
## ZNF91 (GRCh38): chr19:23,393,838-23,432,224. Use +/- 1 Mb flank
## (matches `flank = 1e6` in scripts/locuszoom_plot.R).
CHR        <- 19L
GENE       <- "ZNF91"
GENE_ENSG  <- "ENSG00000167232"
REGION     <- c(22.4e6, 24.4e6)
QTL_YLIM   <- c(0, 8)

## SNP labels (from twas colocalization, see locuszoom_plot.R lines 35-43)
SNPS <- tibble::tribble(
  ~rsID,        ~pos,
  "rs695017",   23421111L,    ## gwas coloc'd snp
  "rs167154",   23307401L     ## isoform qtl coloc'd snp
)

## ----- GWAS meta-analysis -----------------------------------------
gwas_meta <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")

gwas_reg <- gwas_meta %>%
  rename(chr = CHR, pos = BP, p = P) %>%
  mutate(chr = suppressWarnings(as.integer(chr)),
         pos = as.integer(pos),
         p   = as.numeric(p)) %>%
  filter(!is.na(chr), chr == CHR,
         !is.na(pos), pos >= REGION[1], pos <= REGION[2],
         !is.na(p), p > 0)

## ----- AnnotationHub (gene track) ---------------------------------
ah <- AnnotationHub()
ensDb_v113 <- ah[["AH119325"]]

## ----- Modality config --------------------------------------------
mods <- tibble::tribble(
  ~mod,          ~label,             ~color,
  "isoforms",    "Isoform QTLs",     "#377eb8",
  "expression",  "Expression QTLs",  "#e41a1c",
  "alt_polyA",   "Alt_polyA QTLs",   "#ff7f00",
  "alt_TSS",     "Alt_TSS QTLs",     "#984ea3",
  "RNA_editing", "RNA editing QTLs", "#984210",
  "splicing",    "Splicing QTLs",    "#4daf4a",
  "stability",   "Stability QTLs",   "#fdc11c"
) %>%
  mutate(parquet = sprintf(
    "Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/%s/vestibular_%s.cis_qtl_pairs.19.parquet",
    mod, mod))

## ----- Top phenotype per modality ---------------------------------
## Prefer phenotypes for ZNF91 (ENSG00000167232) — their cis windows
## are centered on ZNF91 and cover the full plot region. If a modality
## has no ZNF91 phenotype in the region, fall back to the top signal.
top_pheno_in_region <- function(parquet_path, label,
                                gene_ensg = GENE_ENSG) {
  qtl <- read_parquet(parquet_path) %>%
    left_join(gwas_meta %>% select(variant_id = SNP, rsID), by = "variant_id") %>%
    separate(col = variant_id,
             into = c("vchr", "vpos", "A1", "A2"),
             sep = "_", remove = FALSE) %>%
    mutate(vchr = suppressWarnings(as.integer(vchr)),
           vpos = as.integer(vpos)) %>%
    filter(!is.na(vchr), vchr == CHR,
           !is.na(vpos), vpos >= REGION[1], vpos <= REGION[2],
           !is.na(pval_nominal), pval_nominal > 0)

  if (nrow(qtl) == 0L) {
    warning(sprintf("[%s] no SNPs in region", label))
    return(tibble(pos = integer(), p = double(), rsID = character(),
                  modality = character(), phenotype_id = character()))
  }

  qtl_gene <- qtl %>% filter(grepl(paste0("^", gene_ensg), phenotype_id))
  if (nrow(qtl_gene) > 0L) {
    qtl_use <- qtl_gene
    src <- "ZNF91"
  } else {
    qtl_use <- qtl
    src <- "top-in-region"
  }

  top <- qtl_use %>%
    slice_min(pval_nominal, n = 1, with_ties = FALSE) %>%
    pull(phenotype_id)
  cat(sprintf("  [%s] phenotype (%s): %s\n", label, src, top))

  qtl_use %>%
    filter(phenotype_id == top) %>%
    transmute(pos = vpos, p = pval_nominal, rsID = rsID,
              modality = label, phenotype_id = phenotype_id)
}

cat(sprintf("[..] selecting phenotype per modality in chr%d:%.1f-%.1f Mb\n",
            CHR, REGION[1] / 1e6, REGION[2] / 1e6))

modality_data <- purrr::pmap_dfr(
  list(mods$parquet, mods$label),
  function(pq, lab) top_pheno_in_region(pq, lab)
) %>%
  mutate(modality = factor(modality, levels = mods$label))

## ----- Gene symbol lookup for each chosen phenotype ---------------
get_symbol <- function(ensg) {
  if (is.na(ensg) || ensg == "") return(NA_character_)
  tryCatch({
    g <- ensembldb::genes(ensDb_v113, filter = ~ gene_id == ensg)
    sym <- if (length(g) > 0L) g$gene_name[1] else NA_character_
    if (is.na(sym) || sym == "") NA_character_ else sym
  }, error = function(e) NA_character_)
}

## ----- Per-modality short title -----------------------------------
## Full phenotype IDs are long and overlap data peaks (esp. Isoform).
## Per-modality customization to keep titles compact:
##   Isoform / Alt_TSS  -> transcript ID only          ENST...
##   RNA editing / Alt_polyA / Expression / Stability -> gene (symbol)
##   Splicing           -> intron coords               chr19:start - end Mb
make_title <- function(modality_label, phenotype_id, symbol, ensg) {
  mod <- as.character(modality_label)
  if (mod %in% c("Isoform QTLs", "Alt_TSS QTLs")) {
    m <- regmatches(phenotype_id, regexpr("ENST[0-9]+", phenotype_id))
    short <- if (length(m) > 0L) m else phenotype_id
    return(sprintf("%s: %s", mod, short))
  }
  if (mod == "Splicing QTLs") {
    suffix <- sub(".*__", "", phenotype_id)
    parts  <- strsplit(suffix, "_")[[1]]
    start_v <- suppressWarnings(as.integer(parts[2]))
    end_v   <- suppressWarnings(as.integer(parts[3]))
    if (length(parts) >= 3L && !is.na(start_v) && !is.na(end_v)) {
      return(sprintf("%s: chr%s:%.3f - %.3f Mb", mod,
                     parts[1], start_v / 1e6, end_v / 1e6))
    }
    return(sprintf("%s: %s", mod, phenotype_id))
  }
  # RNA editing, Alt_polyA, Expression, Stability and any fallback:
  # gene-level label
  if (is.na(symbol)) return(sprintf("%s: %s", mod, ensg))
  sprintf("%s: %s (%s)", mod, ensg, symbol)
}

mod_pheno <- modality_data %>%
  distinct(modality, phenotype_id) %>%
  mutate(ensg   = sub("__.*", "", phenotype_id),
         symbol = vapply(ensg, get_symbol, character(1))) %>%
  rowwise() %>%
  mutate(panel_title = make_title(modality, phenotype_id, symbol, ensg)) %>%
  ungroup()

cat("[..] panel titles:\n")
for (i in seq_len(nrow(mod_pheno))) cat("  ", mod_pheno$panel_title[i], "\n")

## ----- Plot helpers -----------------------------------------------
## Generous top/bottom plot margins so y-axis tick labels at panel
## boundaries (the "0" of one and the "8"/"25" of the next) do not
## visually collide when patchwork stacks the panels.
panel_theme <- theme_classic(base_size = 8) +
  theme(
    axis.title    = element_text(size = 8),
    axis.text     = element_text(size = 7, colour = "black"),
    axis.line     = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks    = element_line(linewidth = 0.3, colour = "black"),
    plot.title    = element_blank(),
    plot.margin   = margin(7, 6, 4, 4)
  )

x_indent_mb <- REGION[1] / 1e6 + 0.02 * diff(REGION) / 1e6

build_panel <- function(df, title, color, ylim,
                        label_snps = NULL, show_x = FALSE) {
  if (nrow(df) == 0L) {
    return(ggplot() + ggtitle(paste0(title, " (no data)")) + panel_theme)
  }
  title_y <- ylim[2] * 0.95
  p <- ggplot(df, aes(x = pos / 1e6, y = -log10(p))) +
    geom_point(size = 0.45, alpha = 0.85, colour = color) +
    annotate("text", x = x_indent_mb, y = title_y, label = title,
             hjust = 0, vjust = 1, size = 2.4, fontface = "plain") +
    scale_x_continuous(limits = REGION / 1e6,
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = ylim,
                       expand = expansion(mult = c(0, 0)),
                       oob = scales::oob_squish) +
    ylab(expression(-log[10] ~ pvalue)) +
    panel_theme
  if (!is.null(label_snps) && nrow(label_snps) > 0L) {
    p <- p + ggrepel::geom_text_repel(
      data = label_snps %>% mutate(pos_mb = pos / 1e6, y = -log10(p)),
      aes(x = pos_mb, y = y, label = rsID),
      inherit.aes = FALSE,
      size = 2.4,
      segment.colour = "grey40",
      segment.size = 0.3,
      min.segment.length = 0,
      box.padding = 0.4,
      max.overlaps = Inf,
      seed = 7
    )
  }
  if (!show_x) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())
  } else {
    p <- p + xlab(paste0("Chromosome ", CHR, " (Mb)"))
  }
  p
}

## ----- Build panels -----------------------------------------------
gwas_ymax <- ceiling(max(-log10(gwas_reg$p), na.rm = TRUE) / 5) * 5
gwas_labels <- gwas_reg %>%
  filter(rsID %in% SNPS$rsID) %>%
  distinct(rsID, pos, p)

p_gwas <- build_panel(
  df         = gwas_reg %>% select(pos, p),
  title      = "GWAS meta-analysis",
  color      = "black",
  ylim       = c(0, gwas_ymax),
  label_snps = gwas_labels,
  show_x     = FALSE
)

isoform_df <- modality_data %>% filter(modality == "Isoform QTLs")
isoform_labels <- isoform_df %>%
  filter(pos %in% SNPS$pos[SNPS$rsID == "rs167154"]) %>%
  mutate(rsID = "rs167154") %>%
  distinct(rsID, pos, p)

mod_titles <- setNames(mod_pheno$panel_title, as.character(mod_pheno$modality))

mod_panels <- lapply(seq_len(nrow(mods)), function(i) {
  df  <- modality_data %>% filter(modality == mods$label[i])
  ttl <- mod_titles[[mods$label[i]]]
  if (is.null(ttl) || is.na(ttl)) ttl <- mods$label[i]
  lbl <- if (mods$label[i] == "Isoform QTLs") isoform_labels else NULL
  build_panel(df, ttl, mods$color[i], QTL_YLIM,
              label_snps = lbl, show_x = FALSE)
})
names(mod_panels) <- mods$mod

## ----- Gene track at the bottom -----------------------------------
loc_for_genetrack <- locus(
  data    = gwas_reg %>% transmute(rsid = SNP, chr = chr, pos = pos, p = p),
  gene    = GENE,
  flank   = 1e6,
  ens_db  = ensDb_v113
)
p_genes <- gg_genetracks(loc_for_genetrack,
                         filter_gene_biotype = "protein_coding") +
  xlab(paste0("Chromosome ", CHR, " (Mb)")) +
  theme(plot.margin = margin(7, 6, 4, 4))

## ----- Align all panels so plot regions sit at the same x --------
## Each panel's left-axis label width is different (e.g. "30" wider
## than "8"), so saving panels independently puts the plot region at
## different x positions in each PDF — they don't stack cleanly in
## Illustrator. cowplot::align_plots pads every panel's left margin
## to the global maximum so plot regions line up. We use the aligned
## gtables both for the individual PDFs and the composite, so all
## outputs share the same alignment.
panels    <- c(list(p_gwas), mod_panels, list(p_genes))
n_qtl     <- length(mod_panels)
all_names <- c("gwas", mods$mod, "genetrack")

aligned <- cowplot::align_plots(plotlist = panels, align = "v", axis = "lr")

PANEL_W      <- 6
PANEL_H_QTL  <- 1.5
PANEL_H_GWAS <- 1.9
PANEL_H_GEN  <- 1.9
all_heights  <- c(PANEL_H_GWAS, rep(PANEL_H_QTL, n_qtl), PANEL_H_GEN)

## ----- Save individual aligned panels -----------------------------
## Draw the aligned gtables directly (bypassing ggdraw) so the
## internal layout cowplot computed — same panel-x position across
## panels — is preserved exactly in each saved file. All files share
## PANEL_W so absolute positions match too.
for (i in seq_along(aligned)) {
  pdf_path <- file.path(OUT_DIR, paste0("Supp_Fig3_panel_", all_names[i], ".pdf"))
  png_path <- file.path(OUT_DIR, paste0("Supp_Fig3_panel_", all_names[i], ".png"))

  cairo_pdf(pdf_path, width = PANEL_W, height = all_heights[i])
  grid::grid.newpage(); grid::grid.draw(aligned[[i]])
  dev.off()

  png(png_path, width = PANEL_W, height = all_heights[i],
      units = "in", res = 600)
  grid::grid.newpage(); grid::grid.draw(aligned[[i]])
  dev.off()
}

## ----- Compose vertically (same alignment) ------------------------
composite <- cowplot::plot_grid(plotlist = panels, ncol = 1,
                                align = "v", axis = "lr",
                                rel_heights = c(1.4, rep(1, n_qtl), 1.4))

W_IN <- 6
H_IN <- 14

ggsave(file.path(OUT_DIR, "Supp_Fig3_regional_plots.pdf"), composite,
       width = W_IN, height = H_IN, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Supp_Fig3_regional_plots.png"), composite,
       width = W_IN, height = H_IN, dpi = 600)

cat("[done] wrote:\n",
    "  composite : ", file.path(OUT_DIR, "Supp_Fig3_regional_plots.{pdf,png}"), "\n",
    "  panels    : ", file.path(OUT_DIR, "Supp_Fig3_panel_<mod>.{pdf,png}"), "\n", sep = "")
