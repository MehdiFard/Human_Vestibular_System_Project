## ====================================================================
## Figure 2 — five-panel composite (strict top → bottom)
##   a. MAGMA gene-based Manhattan        (built here)
##   b. enrich-in-GTEx bar plot           (partner PDF, embedded)
##   c. GSEM diagram                      (partner PDF, embedded)
##   d. sfFDR variant-level Manhattan     (built here)
##   e. Regional locuszoom (chr18 52–57Mb)(rendered via locuszoomr)
##
## Run inside conda env:  locuszoomr_env2  (has locuszoomr, ggrepel,
## patchwork, etc.). Existing scripts (man_plots.R, locuszoom_plot.R)
## are NOT modified — code below is self-contained.
## ====================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(data.table)
  library(patchwork)
  library(ggrepel)
  library(cowplot)
  library(grid)
  library(png)
  library(locuszoomr)
  library(ensembldb)
  library(AnnotationHub)
  library(rtracklayer)
  ## load dplyr / tidyr LAST so their verbs shadow Bioconductor masks
  library(dplyr)
  library(tidyr)
})

setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry")

OUT_DIR  <- "figures/paper"
EXT_DIR  <- "figures/paper/external"   # partner PDFs live here
TMP_DIR  <- "figures/paper/external/_render"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

PARTNER_ENRICH <- file.path(EXT_DIR, "enrich_in_gtex.pdf")
PARTNER_GSEM   <- file.path(EXT_DIR, "GSEM.pdf")

## ---- Shared palette / theme (matches Figure 1) -------------------
pal <- list(
  band_dark   = "#1F4E79",
  band_light  = "#9DC3E6",
  sig         = "#E63946",
  thr         = "#E63946",
  label_text  = "grey15"
)

BASE_SIZE <- 7
paper_theme <- theme_classic(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    axis.title         = element_text(size = BASE_SIZE),
    axis.title.x       = element_text(margin = margin(t = 3)),
    axis.text          = element_text(size = BASE_SIZE, colour = "black"),
    axis.line          = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length  = unit(0.08, "cm"),
    plot.title         = element_blank(),
    legend.position    = "none",
    plot.margin        = margin(4, 6, 4, 4)
  )

## ====================================================================
## Helper: convert a single-page PDF to high-DPI PNG (system pdftocairo)
## and return a ggplot panel that draws the PNG.
## ====================================================================
pdf_to_panel <- function(pdf_path, dpi = 600) {
  if (!file.exists(pdf_path))
    stop("Missing partner PDF: ", pdf_path,
         "\nDrop the file there, then re-run.")
  stub <- tools::file_path_sans_ext(basename(pdf_path))
  png_path <- file.path(TMP_DIR, paste0(stub, ".png"))
  needs_render <- !file.exists(png_path) ||
    file.info(pdf_path)$mtime > file.info(png_path)$mtime
  if (needs_render) {
    cmd <- sprintf("pdftocairo -png -r %d -singlefile %s %s",
                   dpi, shQuote(pdf_path),
                   shQuote(file.path(TMP_DIR, stub)))
    sys_status <- system(cmd)
    if (sys_status != 0)
      stop("pdftocairo failed for: ", pdf_path)
  }
  img <- png::readPNG(png_path)
  wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
}

## ====================================================================
## Panel a — MAGMA gene-based Manhattan (in ggplot)
## Source data: data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz
## Style: alternating navy / light-blue chromosomes, red Bonferroni
## line, red highlighted significant genes with ggrepel labels.
## ====================================================================
build_magma <- function() {
  magma <- fread("data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz")
  m <- magma %>%
    transmute(
      SYMBOL = as.character(SYMBOL),
      CHR    = suppressWarnings(as.integer(
                 gsub("^chr", "",
                      ifelse(as.character(CHR) == "X", "23",
                             as.character(CHR)),
                      ignore.case = TRUE))),
      BP     = as.numeric(START),
      P      = as.numeric(P)
    ) %>%
    dplyr::filter(!is.na(CHR), CHR >= 1, CHR <= 22, !is.na(BP),
                  !is.na(P), P > 0, !is.na(SYMBOL), SYMBOL != "") %>%
    arrange(CHR, BP)

  chrlen <- m %>% group_by(CHR) %>%
    summarise(maxbp = max(BP), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(off = lag(cumsum(maxbp), default = 0))
  m <- m %>% inner_join(chrlen, by = "CHR") %>%
    mutate(cum_BP = BP + off,
           log10p = -log10(P),
           band   = factor(CHR %% 2))

  axis_set <- m %>% group_by(CHR) %>%
    summarise(center = mean(cum_BP), .groups = "drop")

  sig <- 0.05 / nrow(m)
  hits <- m %>% dplyr::filter(P <= sig)
  ymax <- max(m$log10p, na.rm = TRUE) * 1.1

  ggplot(m, aes(x = cum_BP, y = log10p)) +
    geom_hline(yintercept = -log10(sig), linetype = 2,
               colour = pal$thr, linewidth = 0.35) +
    geom_point(aes(colour = band), size = 0.45, shape = 16,
               alpha = 0.85, show.legend = FALSE) +
    {if (nrow(hits) > 0)
       geom_point(data = hits, colour = pal$sig,
                  size = 0.9, shape = 16, alpha = 0.95)} +
    {if (nrow(hits) > 0)
       geom_text_repel(data = hits, aes(label = SYMBOL),
                       size = 2.0, fontface = "italic",
                       segment.colour = "grey50",
                       segment.size = 0.25,
                       segment.alpha = 0.7,
                       segment.linetype = "dashed",
                       min.segment.length = 0,
                       box.padding = 0.25,
                       max.overlaps = Inf, seed = 7)} +
    annotate("text",
             x = min(m$cum_BP) + diff(range(m$cum_BP)) * 0.005,
             y = ymax * 0.97, label = "MAGMA",
             hjust = 0, vjust = 1,
             size = (BASE_SIZE + 0.5) / .pt,
             fontface = "bold.italic", colour = pal$label_text) +
    scale_colour_manual(values = c("0" = pal$band_dark,
                                   "1" = pal$band_light)) +
    scale_x_continuous(breaks = axis_set$center,
                       labels = axis_set$CHR,
                       expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
    paper_theme +
    theme(axis.text.x = element_text(size = 3.8))
}

## ====================================================================
## Panel d — sfFDR variant-level Manhattan
## Source: /data/hps/assoc/private/gdml/user/twil15/sffdr_for_dizziness/
##         results/sffdr/results.tsv.gz
## Highlights: rs7098100=MLLT10, rs9536418=OLFM4, rs74564337=C17orf58,
##             rs4456576=DCC, rs4800996=TCF4, rs12980498=ZNF91
## ====================================================================
build_sffdr <- function() {
  sff <- fread("/data/hps/assoc/private/gdml/user/twil15/sffdr_for_dizziness/results/sffdr/results.tsv.gz") %>%
    select(SNP = rsID, CHR, BP, P = fp) %>%
    mutate(SNP = as.character(SNP),
           CHR = suppressWarnings(as.integer(CHR)),
           BP  = suppressWarnings(as.numeric(BP)),
           P   = suppressWarnings(as.numeric(P))) %>%
    dplyr::filter(!is.na(SNP), SNP != "", !is.na(CHR), CHR >= 1, CHR <= 22,
                  !is.na(BP), !is.na(P), P > 0) %>%
    arrange(CHR, BP)

  ## thin baseline for plotting
  sff_log10 <- sff[, log10p := -log10(P)]
  hi <- sff_log10[log10p >= 2]
  lo <- sff_log10[log10p <  2]
  if (nrow(lo) > 80000L) lo <- lo[sample(.N, 80000L)]
  s <- rbindlist(list(hi, lo))

  chrlen <- s[, .(maxbp = max(BP)), by = CHR][order(CHR)]
  chrlen[, off := cumsum(as.numeric(maxbp)) - maxbp]
  s <- merge(s, chrlen[, .(CHR, off)], by = "CHR")
  s[, cum_BP := BP + off]
  s[, band := factor(CHR %% 2)]

  axis_set <- s[, .(center = mean(cum_BP)), by = CHR][order(CHR)]

  snp_gene_map <- data.table(
    SNP = c("rs7098100", "rs9536418", "rs74564337",
            "rs4456576", "rs4800996", "rs12980498"),
    gene = c("MLLT10", "OLFM4", "C17orf58", "DCC", "TCF4", "ZNF91")
  )
  hits <- merge(s, snp_gene_map, by = "SNP")

  thr <- 5e-8
  ymax <- max(s$log10p, na.rm = TRUE) * 1.1

  ggplot(s, aes(x = cum_BP, y = log10p)) +
    geom_hline(yintercept = -log10(thr), linetype = 2,
               colour = pal$thr, linewidth = 0.35) +
    geom_point(aes(colour = band), size = 0.30, shape = 16,
               alpha = 0.85, show.legend = FALSE) +
    geom_point(data = hits, colour = pal$sig,
               size = 0.9, shape = 16, alpha = 0.95) +
    geom_text_repel(data = hits, aes(label = gene),
                    size = 2.0, fontface = "italic",
                    segment.colour = "grey50",
                    segment.size = 0.25,
                    segment.alpha = 0.7,
                    segment.linetype = "dashed",
                    min.segment.length = 0,
                    box.padding = 0.30,
                    max.overlaps = Inf, seed = 7) +
    annotate("text",
             x = min(s$cum_BP) + diff(range(s$cum_BP)) * 0.005,
             y = ymax * 0.97, label = "sfFDR",
             hjust = 0, vjust = 1,
             size = (BASE_SIZE + 0.5) / .pt,
             fontface = "bold.italic", colour = pal$label_text) +
    scale_colour_manual(values = c("0" = pal$band_dark,
                                   "1" = pal$band_light)) +
    scale_x_continuous(breaks = axis_set$center,
                       labels = axis_set$CHR,
                       expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(fp)))) +
    paper_theme +
    theme(axis.text.x = element_text(size = 3.8))
}

## ====================================================================
## Panel e — regional locuszoom (chr18 52–57Mb, sfFDR P-values)
## Renders to a PNG once via locuszoomr; subsequent runs reuse it
## (delete the PNG to force a re-render).
## ====================================================================
LZ_PNG <- file.path(TMP_DIR, "locuszoom_sffdr_chr18.png")

render_locuszoom <- function() {
  if (file.exists(LZ_PNG)) return(invisible(NULL))

  cat("[..] rendering regional locuszoom (chr18 52-57Mb) — first run\n")
  sff <- fread("/data/hps/assoc/private/gdml/user/twil15/sffdr_for_dizziness/results/sffdr/results.tsv.gz") %>%
    transmute(rsid = as.character(rsID),
              chr  = as.integer(CHR),
              pos  = as.numeric(BP),
              p    = as.numeric(fp)) %>%
    dplyr::filter(!is.na(rsid), rsid != "", !is.na(chr), chr >= 1, chr <= 22,
                  !is.na(pos), !is.na(p), p > 0)

  reg <- sff %>% dplyr::filter(chr == 18, pos >= 52e6, pos <= 57e6) %>%
    arrange(p)
  stopifnot(nrow(reg) > 0)

  ah <- AnnotationHub()
  ensDb <- ah[["AH119325"]]

  loc <- locus(data = as.data.frame(reg), seqname = 18,
               xrange = c(52e6, 57e6), ens_db = ensDb,
               chrom = "chr", pos = "pos", p = "p", labs = "rsid",
               index_snp = reg$rsid[which.min(reg$p)])
  loc <- link_LD(loc, token = "7d78d9d5e1ed",
                 pop = "EUR", genome_build = "grch38")
  recomb_gr <- rtracklayer::import.bw("src/ucsc_recomb/hg38/recomb1000GAvg.bw")
  loc <- link_recomb(loc, genome = "hg38", recomb = recomb_gr)

  png(LZ_PNG, width = 180, height = 80, units = "mm", res = 600)
  locus_plot(loc,
             labels = c("index", "rs4800996"),
             ylab = expression(-log[10](SFFDR~P)),
             filter_gene_biotype = "protein_coding",
             highlight = c("DCC", "TCF4"))
  dev.off()
}

build_locuszoom_panel <- function() {
  render_locuszoom()
  img <- png::readPNG(LZ_PNG)
  wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
}

## ====================================================================
## Build all five panels
## ====================================================================
cat("[..] building panel a — MAGMA Manhattan\n")
pA <- build_magma()
cat("[..] building panel d — sfFDR Manhattan\n")
pD <- build_sffdr()
cat("[..] embedding panel b — enrich-in-GTEx (partner PDF)\n")
pB <- pdf_to_panel(PARTNER_ENRICH)
cat("[..] embedding panel c — GSEM (partner PDF)\n")
pC <- pdf_to_panel(PARTNER_GSEM)
cat("[..] building panel e — regional locuszoom\n")
pE <- build_locuszoom_panel()

## ====================================================================
## Save each panel as its own file — for layout in Illustrator
## ====================================================================
W_MAN <- 110; H_MAN <- 60      # Manhattan dimensions (user-specified)

save_fig <- function(p, stub, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stub, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(stub, ".png")), p,
         width = w, height = h, units = "mm", dpi = 600)
}

## a. MAGMA Manhattan
save_fig(pA, "panel_a_magma_manhattan", W_MAN, H_MAN)

## d. sfFDR Manhattan
save_fig(pD, "panel_d_sffdr_manhattan", W_MAN, H_MAN)

## b. enrich-gtex — copy partner PDF unchanged (size is fixed)
file.copy(PARTNER_ENRICH,
          file.path(OUT_DIR, "panel_b_enrich_gtex.pdf"),
          overwrite = TRUE)

## c. GSEM — copy partner PDF unchanged (size is fixed)
file.copy(PARTNER_GSEM,
          file.path(OUT_DIR, "panel_c_gsem.pdf"),
          overwrite = TRUE)

## e. Regional locuszoom — re-render directly to its own PDF + PNG so it
## stays vector-quality where possible. Width matches the Manhattans;
## height is locuszoom-natural (taller because of gene tracks).
LZ_W_MM <- 145
LZ_H_MM <- 90

render_locuszoom_standalone <- function(out_pdf, out_png) {
  ## Reuse the locus object cached during composite build if available
  ## by re-running the same pipeline. Inexpensive after first call:
  ## the AnnotationHub cache and LDlink results are local.
  sff <- fread("/data/hps/assoc/private/gdml/user/twil15/sffdr_for_dizziness/results/sffdr/results.tsv.gz") %>%
    transmute(rsid = as.character(rsID),
              chr  = as.integer(CHR),
              pos  = as.numeric(BP),
              p    = as.numeric(fp)) %>%
    dplyr::filter(!is.na(rsid), rsid != "", !is.na(chr), chr >= 1, chr <= 22,
                  !is.na(pos), !is.na(p), p > 0)
  reg <- sff %>% dplyr::filter(chr == 18, pos >= 52e6, pos <= 57e6) %>%
    arrange(p)
  ah <- AnnotationHub()
  ensDb <- ah[["AH119325"]]
  loc <- locus(data = as.data.frame(reg), seqname = 18,
               xrange = c(52e6, 57e6), ens_db = ensDb,
               chrom = "chr", pos = "pos", p = "p", labs = "rsid",
               index_snp = reg$rsid[which.min(reg$p)])
  loc <- link_LD(loc, token = "7d78d9d5e1ed",
                 pop = "EUR", genome_build = "grch38")
  recomb_gr <- rtracklayer::import.bw("src/ucsc_recomb/hg38/recomb1000GAvg.bw")
  loc <- link_recomb(loc, genome = "hg38", recomb = recomb_gr)

  cairo_pdf(out_pdf, width = LZ_W_MM / 25.4, height = LZ_H_MM / 25.4)
  locus_plot(loc,
             labels = c("index", "rs4800996"),
             ylab = expression(-log[10](SFFDR~P)),
             filter_gene_biotype = "protein_coding",
             highlight = c("DCC", "TCF4"))
  dev.off()
  png(out_png, width = LZ_W_MM, height = LZ_H_MM, units = "mm", res = 600)
  locus_plot(loc,
             labels = c("index", "rs4800996"),
             ylab = expression(-log[10](SFFDR~P)),
             filter_gene_biotype = "protein_coding",
             highlight = c("DCC", "TCF4"))
  dev.off()
}

render_locuszoom_standalone(
  file.path(OUT_DIR, "panel_e_locuszoom.pdf"),
  file.path(OUT_DIR, "panel_e_locuszoom.png")
)

cat("[done] five panels written to ", OUT_DIR, ":\n",
    "    panel_a_magma_manhattan.{pdf,png}   (", W_MAN, " x ", H_MAN, " mm)\n",
    "    panel_d_sffdr_manhattan.{pdf,png}   (", W_MAN, " x ", H_MAN, " mm)\n",
    "    panel_b_enrich_gtex.pdf             (partner PDF, copied)\n",
    "    panel_c_gsem.pdf                    (partner PDF, copied)\n",
    "    panel_e_locuszoom.{pdf,png}         (", LZ_W_MM, " x ", LZ_H_MM, " mm)\n",
    sep = "")
