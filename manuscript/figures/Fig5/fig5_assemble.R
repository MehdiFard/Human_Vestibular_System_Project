## Fig 5 assembly, fully vector:
##   a  QTL / GWAS tracks (panels/Fig5a.rds, fig5a_suppfig4_regional.R)
##   b  GWAS LocusZoom plot drawn on the page with locuszoomr (scatter_plot + genetracks), LD with rs695017 (1000G EUR);
##      also written on its own as panels/Fig5b.{pdf,png} (3.6 x 3.8 in)
##   c  LocusCompare, GWAS vs RNA editing QTL (panels/Fig5c.rds, fig5c_locuscompare.R)
## Built at 7.2 x 5.6 in, then scaled with ghostscript (vector) to the manuscript size: Fig5_6.45x5.02in.{pdf,png}.
## Inputs: data/cache_znf91_region.rds, data/ld_rs695017_1kgEUR.tsv, the two panel objects. Run from the project root:
##   Rscript manuscript/figures/Fig5/fig5_assemble.R
suppressPackageStartupMessages({ library(locuszoomr); library(ensembldb); library(AnnotationHub); library(rtracklayer)
  library(data.table); library(dplyr); library(png); library(grid); library(ggplot2); library(cowplot) })
F5  <- "manuscript/figures/Fig5/data"
OUT <- "manuscript/figures/Fig5"
FONT <- "Nimbus Sans"
W <- 7.2; H <- 5.6; WB <- W / 2; HB <- 3.8                 # page at 1x, and the panel b slot (top right)
## base-graphics text = 12 pt * BASE_CEX * cex (b is drawn at 1x: BASE_CEX 0.5)
## (device pointsize stays 12: panel a's gene labels inherit it)
BASE_CEX <- 0.5
CEX_TEXT <- 0.75; CEX_AXIS <- 1.15; CEX_LAB <- 1.25          # locus_plot defaults 0.7 / 0.9 / 1 (axes ~25% larger; gene labels as in the paper)
GENE_FRAC <- 0.4                                             # share of the b slot for the gene track (as locus_plot)

## ---- panel b data
gw <- readRDS(file.path(F5, "cache_znf91_region.rds"))$gwas %>%
  transmute(SNP, rsID, chr = CHR, pos = BP, A1, A2, beta = BETA, se = SE, p = P) %>% arrange(p)
ld <- fread(file.path(F5, "ld_rs695017_1kgEUR.tsv"))
gw <- gw %>% left_join(ld %>% select(pos, a1, a2, r2), by = "pos", relationship = "many-to-many") %>%
  filter(is.na(a1) | (A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)) %>%
  group_by(SNP) %>% slice(1) %>% ungroup() %>% select(-a1, -a2) %>% as.data.frame()
ensDb <- AnnotationHub(localHub = TRUE)[["AH119325"]]
loc <- locus(data = gw, gene = "ZNF91", flank = 2.5e5, ens_db = ensDb, labs = "rsID", index_snp = "rs695017")
loc$data$ld <- loc$data$r2
loc <- link_recomb(loc, genome = "hg38", recomb = import.bw("resources/ucsc_recomb/hg38/recomb1000GAvg.bw"))

pa <- readRDS("manuscript/figures/Fig5/panels/Fig5a.rds")
pc <- readRDS("manuscript/figures/Fig5/panels/Fig5c.rds")

## x offset (in) of c's left y-axis from c's left edge, from its 3.6 x 1.8 in render (first long dark vertical line)
axis_in <- function(img, width_in, min_run = 0.3) {
  dark <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3 < 0.35
  which(colMeans(dark) > min_run)[1] / ncol(img) * width_in
}
fc_in <- axis_in(readPNG("manuscript/figures/Fig5/panels/Fig5c.png"), WB)

## b = locus_plot(loc, labels, label_x) split by hand: genetracks (bottom 2/5) + scatter_plot (top 3/5),
## drawn in the device region x0-x1, y0-y1 (fractions of the page)
draw_b <- function(x0, x1, y0, y1) {
  ys <- y0 + (y1 - y0) * GENE_FRAC
  par(fig = c(x0, x1, y0, ys), new = TRUE); par(cex = BASE_CEX)
  genetracks(loc, cex.axis = CEX_AXIS, cex.lab = CEX_LAB, cex.text = CEX_TEXT, gene_col = "blue4", exon_col = "blue4",
             exon_border = "blue4", xticks = TRUE, showRecomb = TRUE, maxrows = 7)
  par(fig = c(x0, x1, ys, y1), new = TRUE); par(cex = BASE_CEX)
  scatter_plot(loc, xticks = FALSE, cex.axis = CEX_AXIS, cex.lab = CEX_LAB, recomb_col = "blue",
               labels = c("rs695017", "rs167154", "rs607168"), label_x = c(-1, -5, 4))
}

compose <- function() {
  par(family = FONT); plot.new()                              # blank page for base graphics; grid draws on top
  draw_b(0.5, 1, (H - HB) / H, 1)
  ## page x (in) of b's y-axis, read from the scatter plot just drawn
  axis_b <- grconvertX(par("usr")[1], "user", "ndc") * W
  ## a (left column) and c (bottom right, its y-axis aligned with b's)
  print(pa, vp = viewport(0, 0, WB, H, just = c(0, 0), default.units = "in"))
  xc <- axis_b - fc_in
  print(pc, vp = viewport(xc, 0, W - xc, H - HB, just = c(0, 0), default.units = "in"))
  grid.text(c("a", "b", "c"), unit(c(0.05, WB + 0.02, WB + 0.02), "in"), unit(c(H - 0.1, H - 0.1, H - HB - 0.08), "in"),
            just = c(0, 1), gp = gpar(fontfamily = FONT, fontface = "bold", fontsize = 10))
  cat(sprintf("b y-axis at %.3f in; c placed at x = %.3f in (axis offset %.3f in)\n", axis_b, xc, fc_in))
}
f1 <- tempfile(fileext = ".pdf")                      # intermediate (vector, 7.2 x 5.6 in)
cairo_pdf(f1, width = W, height = H, family = FONT); compose(); invisible(dev.off())

## ---- panel b on its own, same size and text as on the page
fb <- file.path(OUT, "panels", "Fig5b")
cairo_pdf(paste0(fb, ".pdf"), width = WB, height = HB, family = FONT)
par(family = FONT); plot.new(); draw_b(0, 1, 0, 1); invisible(dev.off())
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=png16m", "-r600", "-dTextAlphaBits=4",
                "-dGraphicsAlphaBits=4", paste0("-sOutputFile=", fb, ".png"), paste0(fb, ".pdf")))

## ---- vector scaling (ghostscript) to the manuscript size, 6.45 x 5.02 in (same aspect ratio), + 600 dpi PNG
## (DPI tag set, so Word places it at exactly w x h in)
scale_to <- function(src, stem, w, h) {
  system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=pdfwrite", "-dFIXEDMEDIA", "-dPDFFitPage",
                  "-dAutoRotatePages=/None", "-dEmbedAllFonts=true", sprintf("-dDEVICEWIDTHPOINTS=%.2f", w * 72),
                  sprintf("-dDEVICEHEIGHTPOINTS=%.2f", h * 72), paste0("-sOutputFile=", stem, ".pdf"), src))
  system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=png16m", "-r600", "-dTextAlphaBits=4",
                  "-dGraphicsAlphaBits=4", paste0("-sOutputFile=", stem, ".png"), paste0(stem, ".pdf")))
}
WORD <- c(6.45, 5.02)
scale_to(f1, file.path(OUT, "Fig5_6.45x5.02in"), WORD[1], WORD[2])
cat("[done]\n")


### END