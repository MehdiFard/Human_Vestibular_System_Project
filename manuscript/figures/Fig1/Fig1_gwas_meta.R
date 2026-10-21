## ====================================================================
## Figure 1 (lower row) — Manhattan + QQ plots, one PDF/PNG per panel
##
## Six standalone files written to manuscript/figures/Fig1/panels/:
##   manhattan_meta_all.{pdf,png}    qq_meta_all.{pdf,png}
##   manhattan_meta_eur.{pdf,png}    qq_meta_eur.{pdf,png}
##   manhattan_gwas_ukb.{pdf,png}    qq_gwas_ukb.{pdf,png}
##
## Source data: data/gwas_raw_sumstats/ (see data/README.md).
## ====================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(data.table)
})


OUT_DIR <- "manuscript/figures/Fig1/panels"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Palette (refined) -------------------------------------------
## Matches the demographic palette family but slightly punchier:
##   - lighter alternating band so the signal dots pop more
##   - vivid red for significant SNPs and threshold line
##   - pale-blue CI ribbon (instead of grey) for visual cohesion
pal <- list(
  band_dark   = "#1F4E79",
  band_light  = "#9DC3E6",
  sig         = "#E63946",
  thr         = "#E63946",
  qq_point    = "#1F4E79",
  qq_ribbon   = "#DCE6F0"
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
## DATA — load each GWAS, keep CHR/BP/P, drop unusable rows
## ====================================================================
gwas_files <- list(
  Meta_All = list(
    path = "data/gwas_raw_sumstats/balance_meta1.tbl.sorted.filtered.gz",
    cols = c(SNP = "MarkerName", CHR = "Chromosome", BP = "Position", P = "P-value"),
    label = "Meta_All",
    file_stub = "meta_all"
  ),
  Meta_Eur = list(
    path = "data/gwas_raw_sumstats/balance_meta_EUR.tbl.sorted.tab.gz",
    cols = c(SNP = "SNP", CHR = "Chromosome", BP = "Position", P = "P-value"),
    label = "Meta_Eur",
    file_stub = "meta_eur"
  ),
  GWAS_UKB = list(
    path = "data/gwas_raw_sumstats/UKB_EUR_balance.txt.gz",
    cols = c(SNP = "SNP", CHR = "CHR", BP = "BP", P = "P"),
    label = "GWAS_UKB",
    file_stub = "gwas_ukb"
  )
)

load_gwas <- function(spec) {
  cols <- spec$cols
  dt <- fread(spec$path, select = unname(cols))
  setnames(dt, unname(cols), names(cols))
  dt <- dt[, .(SNP = as.character(SNP),
               CHR = suppressWarnings(as.integer(CHR)),
               BP  = suppressWarnings(as.numeric(BP)),
               P   = suppressWarnings(as.numeric(P)))]
  dt <- dt[!is.na(P) & P > 0 & P <= 1 & !is.na(CHR) &
             CHR >= 1 & CHR <= 22 & !is.na(BP)]
  dt[, log10p := -log10(P)]
  dt
}

cat("[..] loading GWAS sumstats\n")
gwas_list <- list()
for (nm in names(gwas_files)) {
  cat("    -> ", nm, "\n", sep = "")
  gwas_list[[nm]] <- load_gwas(gwas_files[[nm]])
  cat("       n =", format(nrow(gwas_list[[nm]]), big.mark = ","), "\n")
}

## ---- Common chromosome offsets (union of max BP per chr) ---------
chr_max_global <- rbindlist(lapply(gwas_list,
                                   function(d) d[, .(chr_max = max(BP)), by = CHR]))[
  , .(chr_max = max(chr_max)), by = CHR][order(CHR)]
chr_max_global[, chr_offset := cumsum(as.numeric(chr_max)) - chr_max]
chr_axis <- chr_max_global[, .(CHR,
                                center = chr_offset + chr_max / 2,
                                left   = chr_offset,
                                right  = chr_offset + chr_max)]

GW_THR <- 5e-8
GW_LOG <- -log10(GW_THR)

## ====================================================================
## MANHATTAN
## ====================================================================
set.seed(1)   # the background (p > 0.01) Manhattan points are a random subsample
thin_manhattan <- function(dt, keep_above = 2, base_n = 80000L) {
  hi <- dt[log10p >= keep_above]
  lo <- dt[log10p <  keep_above]
  if (nrow(lo) > base_n) lo <- lo[sample(.N, base_n)]
  rbindlist(list(hi, lo))
}

attach_position <- function(dt) {
  dt <- merge(dt, chr_max_global[, .(CHR, chr_offset)], by = "CHR")
  dt[, cum_BP := BP + chr_offset]
  dt
}

manhattan_plot <- function(dt, label) {
  d <- attach_position(thin_manhattan(dt))
  d[, chr_band := factor(CHR %% 2)]
  d_sig <- d[log10p >= GW_LOG]
  ymax <- max(d$log10p, na.rm = TRUE) * 1.04

  band_cols <- c("0" = pal$band_dark, "1" = pal$band_light)

  p <- ggplot(d, aes(x = cum_BP, y = log10p)) +
    geom_hline(yintercept = GW_LOG, linetype = 2,
               colour = pal$thr, linewidth = 0.35) +
    geom_point(aes(colour = chr_band),
               size = 0.35, shape = 16, alpha = 0.85, show.legend = FALSE)
  if (nrow(d_sig) > 0) {
    p <- p + geom_point(data = d_sig, colour = pal$sig,
                        size = 0.7, shape = 16, alpha = 0.95)
  }
  x_span <- diff(range(chr_axis$left, chr_axis$right))
  p +
    annotate("text",
             x = chr_axis$left[1] + x_span * 0.04,
             y = ymax * 0.95, label = label,
             hjust = 0, vjust = 1,
             size = (BASE_SIZE + 0.5) / .pt, fontface = "bold.italic",
             colour = "grey15") +
    scale_colour_manual(values = band_cols) +
    scale_x_continuous(breaks = chr_axis$center,
                       labels = chr_axis$CHR,
                       expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Chromosome",
         y = expression(-log[10](italic(p)))) +
    paper_theme +
    theme(axis.text.x = element_text(size = 4.5))
}

## ====================================================================
## QQ
## ====================================================================
qq_data <- function(dt, top_keep = 3000L, n_bins = 1500L) {
  obs <- sort(dt$log10p, decreasing = TRUE)
  n   <- length(obs)
  expd <- sort(-log10(ppoints(n)), decreasing = TRUE)
  qq <- data.table(observed = obs, expected = expd, rank = seq_len(n))
  i <- qq$rank
  qq[, ci_lo := -log10(qbeta(0.975, i, n - i + 1))]
  qq[, ci_hi := -log10(qbeta(0.025, i, n - i + 1))]
  hi <- qq[seq_len(min(top_keep, n))]
  if (n > top_keep) {
    rest <- qq[(top_keep + 1):n]
    rest[, bin := cut(expected, breaks = n_bins,
                      labels = FALSE, include.lowest = TRUE)]
    rest <- rest[, .(observed = mean(observed), expected = mean(expected),
                     ci_lo = mean(ci_lo), ci_hi = mean(ci_hi),
                     rank = round(mean(rank))), by = bin][, bin := NULL]
    out <- rbindlist(list(hi[, .(observed, expected, ci_lo, ci_hi, rank)], rest))
  } else {
    out <- hi[, .(observed, expected, ci_lo, ci_hi, rank)]
  }
  out[order(-rank)]
}

genomic_inflation <- function(dt) {
  chisq <- qchisq(1 - dt$P, 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
}

qq_plot <- function(dt) {
  qq <- qq_data(dt)
  lambda <- genomic_inflation(dt)
  xmax <- max(qq$expected, na.rm = TRUE) * 1.02
  ymax <- max(qq$observed, na.rm = TRUE) * 1.02
  ggplot(qq, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = pal$qq_ribbon, alpha = 0.85) +
    geom_abline(intercept = 0, slope = 1, linetype = 2,
                colour = pal$thr, linewidth = 0.35) +
    geom_point(colour = pal$qq_point, size = 0.6, shape = 16, alpha = 0.9) +
    annotate("text", x = xmax * 0.08, y = ymax * 0.95,
             label = sprintf("lambda[GC] == %.2f", lambda),
             parse = TRUE, hjust = 0, vjust = 1,
             size = (BASE_SIZE + 0.5) / .pt, colour = "grey15") +
    scale_x_continuous(limits = c(0, xmax),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = expression(Expected~-log[10](italic(p))),
         y = expression(Observed~-log[10](italic(p)))) +
    paper_theme
}

## ====================================================================
## RENDER & SAVE — six files
## ====================================================================
W_MAN <- 145; H_MAN <- 45     # Manhattan: shorter width, short
W_QQ  <- 45;  H_QQ  <- 45     # QQ: small square

save_fig <- function(p, stub, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stub, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(stub, ".png")), p,
         width = w, height = h, units = "mm", dpi = 600)
}

cat("[..] rendering Manhattan + QQ panels (6 files)\n")
for (nm in names(gwas_files)) {
  spec <- gwas_files[[nm]]
  pM <- manhattan_plot(gwas_list[[nm]], spec$label)
  pQ <- qq_plot(gwas_list[[nm]])
  save_fig(pM, paste0("manhattan_", spec$file_stub), W_MAN, H_MAN)
  save_fig(pQ, paste0("qq_",        spec$file_stub), W_QQ,  H_QQ)
  cat("    saved manhattan_", spec$file_stub, ".{pdf,png} + qq_",
      spec$file_stub, ".{pdf,png}\n", sep = "")
}

cat("[done] all six figures written to ", OUT_DIR, "\n", sep = "")


### END
