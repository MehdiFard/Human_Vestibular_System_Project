## ====================================================================
## Figure 2 — panels, each saved to panels/ (PDF + PNG) and assembled in Illustrator
##   a. MAGMA gene-based Manhattan        (built here)
##   b. enrich-in-GTEx bar plot           (partner figure = Supp Fig 2, not built here)
##   c. GSEM diagram                      (partner PDF in inputs/, copied + PNG)
##   d. sfFDR variant-level Manhattan     (built here)
##   e. Regional locuszoom (chr18 52–57Mb)(locuszoomr; LD from NIH LDlink, token in config_lab.yaml)
##
## Environment: locuszoomr (config envs.locuszoom). Panel c is a partner figure (inputs/, not in git);
## d and e use analyses/10_sffdr results.
## ====================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(data.table)
  library(ggrepel)
  library(cowplot)
  library(locuszoomr)
  library(ensembldb)
  library(AnnotationHub)
  library(rtracklayer)
  ## load dplyr / tidyr LAST so their verbs shadow Bioconductor masks
  library(dplyr)
  library(tidyr)
})


cfg <- yaml::read_yaml("config/config.yaml")   # + lab-specific values (config_lab.yaml, not in git)
if (file.exists("config/config_lab.yaml")) cfg <- modifyList(cfg, yaml::read_yaml("config/config_lab.yaml"))
LDLINK_TOKEN <- cfg$apis$ldlink$token
if (is.null(LDLINK_TOKEN) || LDLINK_TOKEN == "") stop("set apis: ldlink: token: in config/config_lab.yaml (personal NIH LDlink token)")

OUT_DIR  <- "manuscript/figures/Fig2/panels"
EXT_DIR  <- "manuscript/figures/Fig2/inputs"   # partner PDF (GSEM path diagram)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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
## Panel a — MAGMA gene-based Manhattan (in ggplot)
## Source data: data/gwas_raw_sumstats/magma_all_sets.csv.gz
## Style: alternating navy / light-blue chromosomes, red Bonferroni
## line, red highlighted significant genes with ggrepel labels.
## ====================================================================
build_magma <- function() {
  magma <- fread("data/gwas_raw_sumstats/magma_all_sets.csv.gz")
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
## Source: analyses/10_sffdr (sfFDR workflow results)
##         results/sffdr/results.tsv.gz
## Highlights: rs7098100=MLLT10, rs9536418=OLFM4, rs74564337=C17orf58,
##             rs4456576=DCC, rs4800996=TCF4, rs12980498=ZNF91
## ====================================================================
build_sffdr <- function() {
  sff <- fread("analyses/10_sffdr/results/sffdr/results.tsv.gz") %>%
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
## Build panels a and d
## ====================================================================
cat("[..] building panel a — MAGMA Manhattan\n")
pA <- build_magma()
cat("[..] building panel d — sfFDR Manhattan\n")
pD <- build_sffdr()

## ====================================================================
## Save each panel as its own file — for layout in Illustrator
## ====================================================================
W_MAN <- 110; H_MAN <- 50      # Manhattan dimensions (user-specified)

save_fig <- function(p, stub, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stub, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(stub, ".png")), p,
         width = w, height = h, units = "mm", dpi = 600)
}

## a. MAGMA Manhattan
save_fig(pA, "panel_a_magma_manhattan_new", W_MAN, H_MAN)

## d. sfFDR Manhattan
save_fig(pD, "panel_d_sffdr_manhattan_new", W_MAN, H_MAN)

## c. GSEM — copy partner PDF unchanged (size is fixed) + 600 dpi PNG
if (!file.exists(PARTNER_GSEM)) stop("Missing partner PDF: ", PARTNER_GSEM)
file.copy(PARTNER_GSEM, file.path(OUT_DIR, "panel_c_gsem.pdf"), overwrite = TRUE)
if (system2("pdftocairo", c("-png", "-r", "600", "-singlefile", shQuote(PARTNER_GSEM),
                            shQuote(file.path(OUT_DIR, "panel_c_gsem")))) != 0)
  stop("pdftocairo failed for: ", PARTNER_GSEM)

## e. Regional locuszoom (chr18 52–57Mb, sfFDR P-values), drawn directly to PDF + PNG.
## Width matches the Manhattans; height is locuszoom-natural (taller because of gene tracks).
LZ_W_MM <- 145
LZ_H_MM <- 90

cat("[..] building panel e — regional locuszoom\n")
sff <- fread("analyses/10_sffdr/results/sffdr/results.tsv.gz") %>%
  transmute(rsid = as.character(rsID),
            chr  = as.integer(CHR),
            pos  = as.numeric(BP),
            p    = as.numeric(fp)) %>%
  dplyr::filter(!is.na(rsid), rsid != "", !is.na(chr), chr >= 1, chr <= 22,
                !is.na(pos), !is.na(p), p > 0)
reg <- sff %>% dplyr::filter(chr == 18, pos >= 52e6, pos <= 57e6) %>%
  arrange(p)
stopifnot(nrow(reg) > 0)
ensDb <- AnnotationHub()[["AH119325"]]
loc <- locus(data = as.data.frame(reg), seqname = 18,
             xrange = c(52e6, 57e6), ens_db = ensDb,
             chrom = "chr", pos = "pos", p = "p", labs = "rsid",
             index_snp = reg$rsid[which.min(reg$p)])
loc <- link_LD(loc, token = LDLINK_TOKEN,
               pop = "EUR", genome_build = "grch38")
recomb_gr <- rtracklayer::import.bw("resources/ucsc_recomb/hg38/recomb1000GAvg.bw")
loc <- link_recomb(loc, genome = "hg38", recomb = recomb_gr)

draw_lz <- function() locus_plot(loc,
                                 labels = c("index", "rs4800996"),
                                 ylab = expression(-log[10](SFFDR~P)),
                                 filter_gene_biotype = "protein_coding",
                                 highlight = c("DCC", "TCF4"))
cairo_pdf(file.path(OUT_DIR, "panel_e_locuszoom.pdf"), width = LZ_W_MM / 25.4, height = LZ_H_MM / 25.4)
draw_lz(); invisible(dev.off())
png(file.path(OUT_DIR, "panel_e_locuszoom.png"), width = LZ_W_MM, height = LZ_H_MM, units = "mm", res = 600)
draw_lz(); invisible(dev.off())

cat("[done] panels written to ", OUT_DIR, ":\n",
    "    panel_a_magma_manhattan_new.{pdf,png}   (", W_MAN, " x ", H_MAN, " mm)\n",
    "    panel_c_gsem.{pdf,png}                  (partner PDF, copied)\n",
    "    panel_d_sffdr_manhattan_new.{pdf,png}   (", W_MAN, " x ", H_MAN, " mm)\n",
    "    panel_e_locuszoom.{pdf,png}             (", LZ_W_MM, " x ", LZ_H_MM, " mm)\n",
    sep = "")


### END
