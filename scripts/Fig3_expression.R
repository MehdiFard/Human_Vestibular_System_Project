## ====================================================================
## Figure 3 — vestibular transcriptomic signatures (2 x 2 grid)
##   a. PCA: GTEx tissues + VESTIB (log2 CPM, same pipeline)
##   b. Highest GTEx median vs VESTIB median expression
##   c. Enrichment dot plot — Disease / Monarch / UniProt
##   d. Enrichment dot plot — GO Biological Process
##
## Heavy intermediates (annotations for panel a, dt_plot for panel b)
## are cached as RDS files under figures/data_tables/_fig3_cache/.
## First run takes a few minutes; subsequent runs are fast.
##
## Source scripts (kept untouched): PCA_VESTIB_vs_GTEx.R,
## highly_expressed_in_vestib.R, gene_enrichment.R.
## Run in env Env1 (default) or any env with the listed packages.
## ====================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(data.table)
  library(patchwork)
  library(jsonlite)
  library(matrixStats)
  library(stringr)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry")

OUT_DIR <- "figures/paper"
CACHE   <- "figures/data_tables/_fig3_cache"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE,   showWarnings = FALSE, recursive = TRUE)

## ---- Paper palette / theme (matches Fig 1 & Fig 2) ----------------
pal <- list(
  navy        = "#1F4E79",
  blue        = "#9DC3E6",
  red         = "#E63946",   # for VESTIB highlight + threshold
  gold        = "#E0A458",
  grey_blue   = "#9CADBF"
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
    legend.title       = element_text(size = BASE_SIZE - 0.5),
    legend.text        = element_text(size = BASE_SIZE - 0.5),
    legend.key.size    = unit(0.28, "cm"),
    legend.margin      = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.05, "cm"),
    plot.margin        = margin(4, 6, 4, 4)
  )

## ====================================================================
## Helper — cache wrapper
## ====================================================================
cache_or_compute <- function(name, builder) {
  fp <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(fp)) {
    cat("[cache] using ", fp, "\n", sep = "")
    return(readRDS(fp))
  }
  cat("[build] computing ", name, " (slow on first run)\n", sep = "")
  obj <- builder()
  saveRDS(obj, fp)
  obj
}

## ====================================================================
## Shared — GTEx color palette (tissue_abbrv -> hex)
## ====================================================================
load_gtex_palette <- function() {
  cache_or_compute("gtex_palette", function() {
    gx <- fromJSON("src/gtex_colors.json")
    annot <- bind_rows(lapply(names(gx), function(t) {
      data.frame(tissue_name  = t,
                 tissue_abbrv = gx[[t]]$tissue_abbrv,
                 tissue_color_hex = paste0("#", gx[[t]]$tissue_color_hex),
                 stringsAsFactors = FALSE)
    }))
    setNames(annot$tissue_color_hex, annot$tissue_abbrv)
  })
}

## ====================================================================
## Panel a — PCA of logCPM (GTEx + VESTIB, same pipeline)
##
## Reuses figures/data_tables/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData
## (created in PCA_VESTIB_vs_GTEx.R). Sample-to-tissue annotation is
## reconstructed from file headers, so we don't need to re-load any
## count matrices.
## ====================================================================
build_annot_gtex_vest <- function() {
  cache_or_compute("annot_gtex_vest", function() {
    gtex_dir <- "src/GTEx/tissues_raw_count_from_Dan_pantry"
    files <- list.files(gtex_dir, pattern = "^gene_counts_.*\\.tsv\\.gz$",
                        full.names = TRUE)
    exclude <- c("TESTIS", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL", "BLDDER")
    files <- files[!grepl(paste(exclude, collapse = "|"), basename(files))]

    gtex_annot <- rbindlist(lapply(files, function(f) {
      tissue  <- sub("^gene_counts_(.*)\\.tsv\\.gz$", "\\1", basename(f))
      samples <- setdiff(names(fread(f, nrows = 0)), "gene_id")
      data.table(human_id = samples, tissue = tissue, batch = "gtex")
    }))

    vest_meta <- fread("src/metadata_rna.tsv",
                       select = c("human_id", "batch")) %>%
      mutate(human_id = as.character(human_id),
             tissue   = "VESTIB") %>%
      distinct(human_id, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(human_id))

    rbindlist(list(gtex_annot, vest_meta[, c("human_id", "tissue", "batch")]),
              use.names = TRUE)
  })
}

build_pca_panel <- function() {
  ## load the saved PCA object (variable name: "pca")
  pca_rdata <- "figures/data_tables/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData"
  load(pca_rdata)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

  pc_df <- data.frame(human_id = rownames(pca$x),
                      PC1 = pca$x[, 1],
                      PC2 = pca$x[, 2],
                      stringsAsFactors = FALSE)

  annot <- build_annot_gtex_vest()

  pc_df_annot <- pc_df %>%
    inner_join(annot, by = "human_id") %>%
    mutate(tissue = dplyr::if_else(grepl("^BRN", tissue), "BRN_TISSUES", tissue))

  ## drop known outlier vestibular samples (original script: 31907773,
  ## 32073088). Also drop any VESTIB sample with PC1 > 0 — the lone
  ## positive-PC1 VESTIB sample was excluded from downstream analysis.
  pc_df_annot <- pc_df_annot %>%
    dplyr::filter(!human_id %in% c("31907773", "32073088")) %>%
    dplyr::filter(!(tissue == "VESTIB" & PC1 > 0))

  ## palette: load GTEx tissue colors, override VESTIB to make it pop
  gtex_palette <- load_gtex_palette()
  gtex_palette["BRN_TISSUES"] <- "#EEEE00"
  gtex_palette["VESTIB"]      <- pal$red
  gtex_palette["GEUVADIS"]    <- "grey30"

  ggplot(pc_df_annot, aes(PC1, PC2, colour = tissue)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.25, colour = "grey60") +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.25, colour = "grey60") +
    geom_point(size = 0.45, alpha = 0.85) +
    ## draw VESTIB on top so it isn't buried
    geom_point(data = dplyr::filter(pc_df_annot, tissue == "VESTIB"),
               size = 0.75, alpha = 0.95, colour = pal$red) +
    scale_colour_manual(values = gtex_palette, guide = "none") +
    labs(x = sprintf("PC1 (%.1f%%)", var_expl[1] * 100),
         y = sprintf("PC2 (%.1f%%)", var_expl[2] * 100)) +
    paper_theme +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
}

## ====================================================================
## Panel b — Highest GTEx median expression vs VESTIB median
## ====================================================================
build_dt_plot <- function() {
  cache_or_compute("dt_plot_panel_b", function() {
    ## protein-coding + lncRNA gene IDs from GTF
    library(rtracklayer)
    gtf <- readGFF("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf")
    pcg_lnc <- gtf %>%
      dplyr::filter(type == "gene",
                    gene_biotype %in% c("protein_coding", "lncRNA")) %>%
      dplyr::pull(gene_id)

    ## VESTIB medians (drop known outliers)
    expr_vest <- fread("Pantry/phenotyping/output/vestibular/unnorm/expression.bed") %>%
      dplyr::select(-`#chr`, -start, -end,
                    -"31947386", -"31907773", -"32073088") %>%
      dplyr::rename(gene_id = phenotype_id) %>%
      dplyr::filter(gene_id %in% pcg_lnc)
    med_v <- data.table(
      gene_id = expr_vest$gene_id,
      VESTIB  = rowMedians(as.matrix(expr_vest[, -1]))
    )

    ## GTEx medians per tissue
    gtex_root <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/phenos/"
    tissues <- setdiff(list.dirs(gtex_root, recursive = FALSE,
                                  full.names = FALSE),
                       c("ref", "TESTIS", "BLDDER", "CVSEND",
                         "CVXECT", "FLLPNT", "KDNMDL"))
    med_tbl <- med_v
    for (tt in tissues) {
      cat("    GTEx median: ", tt, "\n", sep = "")
      expr_t <- fread(file.path(gtex_root, tt, "output/unnorm/expression.bed")) %>%
        dplyr::select(-`#chr`, -start, -end) %>%
        dplyr::rename(gene_id = phenotype_id) %>%
        dplyr::filter(gene_id %in% pcg_lnc)
      med_t <- data.table(
        gene_id = expr_t$gene_id,
        gtex_tissue = rowMedians(as.matrix(expr_t[, -1]))
      )
      med_tbl <- med_tbl %>%
        inner_join(med_t, by = "gene_id") %>%
        dplyr::filter(VESTIB >= 0.1 | gtex_tissue >= 0.1) %>%
        dplyr::rename(!!tt := gtex_tissue)
    }

    ## filter: expressed (>0.1) in >=10 tissues
    med_filt <- med_tbl %>%
      mutate(n_expr = rowSums(across(-gene_id, ~ . > 0.1))) %>%
      dplyr::filter(n_expr >= 10) %>%
      dplyr::select(-n_expr)

    tissue_cols <- setdiff(names(med_filt), c("gene_id", "VESTIB"))
    dt_max <- med_filt[, {
      vals_log <- log2(c(VESTIB, unlist(.SD, use.names = FALSE)) + 1)
      v_log    <- vals_log[1]
      o_log    <- vals_log[-1]
      j        <- which.max(o_log)
      .(
        VESTIB     = v_log,
        max_other  = o_log[j],
        max_tissue = tissue_cols[j],
        log2FC     = v_log - o_log[j]
      )
    }, by = gene_id, .SDcols = tissue_cols]

    ttest <- fread("figures/data_tables/vestib_vs_gtex_t_tests.tsv")[, c(1, 8)]
    dt_max_t <- dt_max %>% left_join(ttest, by = "gene_id")
    dt_max_t <- dt_max_t[max_other > 0]

    ## significant subset (matches original script's filter)
    dt_max_t %>% dplyr::filter(padj_two_sided < 0.01)
  })
}

build_expr_panel <- function() {
  dt_plot <- build_dt_plot()
  ## Three visual classes:
  ##   - VESTIB-enriched: log2FC > 2 (the cut-off used for enrichment)
  ##   - other significant: padj < 0.01 (already filtered in dt_plot)
  ##   - background / non-significant points are not in dt_plot
  dt_plot[, gene_class := dplyr::if_else(log2FC > 2,
                                         "VESTIB-enriched (log2FC > 2)",
                                         "Other significant")]
  dt_plot[, gene_class := factor(gene_class,
                                 levels = c("Other significant",
                                            "VESTIB-enriched (log2FC > 2)"))]

  n_enriched <- sum(dt_plot$gene_class == "VESTIB-enriched (log2FC > 2)")

  ## Range needed for diagonal-cutoff annotation
  rng <- range(c(dt_plot$max_other, dt_plot$VESTIB), na.rm = TRUE)
  xseq <- seq(rng[1], rng[2], length.out = 50)

  ggplot(dt_plot, aes(x = max_other, y = VESTIB)) +
    ## y = x (equal expression) — grey reference line
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                colour = "grey55", linewidth = 0.35) +
    ## y = x + 2 (log2FC > 2 cut-off) — red dashed
    geom_abline(slope = 1, intercept = 2, linetype = 2,
                colour = pal$red, linewidth = 0.4) +
    ## background of significant points (smaller, faded)
    geom_point(data = dplyr::filter(dt_plot,
                                    gene_class == "Other significant"),
               colour = pal$navy, alpha = 0.30, size = 0.30) +
    ## VESTIB-enriched genes on top, in vivid red
    geom_point(data = dplyr::filter(dt_plot,
                                    gene_class == "VESTIB-enriched (log2FC > 2)"),
               colour = pal$red, alpha = 0.85, size = 0.55) +
    annotate("text",
             x = rng[1] + diff(rng) * 0.02,
             y = rng[2] * 0.97,
             label = sprintf("VESTIB-enriched (log₂FC > 2):\n%d genes",
                             n_enriched),
             hjust = 0, vjust = 1,
             size = (BASE_SIZE - 0.5) / .pt, colour = pal$red) +
    labs(
      x = expression("Max GTEx median ("*log[2]*" TPM)"),
      y = expression("VESTIB median ("*log[2]*" TPM)")
    ) +
    paper_theme
}

## ====================================================================
## Panel c — Disease / Monarch / UniProt enrichment dot plot
## Source: src/new_enrichment_all_databases.csv (from cytoscape export).
## Keeps the same database filter used in gene_enrichment.R.
## ====================================================================
build_disease_panel <- function() {
  csv <- fread("src/new_enrichment_all_databases.csv")
  df <- csv %>%
    dplyr::rename(
      bg_genes      = `background genes`,
      gene_count    = genes,
      category_name = category,
      term          = description,
      FDR           = `FDR value`,
      pval          = `p-value`
    ) %>%
    dplyr::filter(!category_name %in% c(
      "GO Biological Process", "COMPARTMENTS",
      "GO Cellular Component", "GO Molecular Function",
      "SMART Domains", "WikiPathways", "STRING Clusters"
    )) %>%
    dplyr::filter(FDR < 0.05) %>%
    arrange(FDR) %>%
    mutate(
      cat_short     = sub(" .*", "", category_name),
      term_category = paste0(term, " (", cat_short, ")"),
      term_category = stringr::str_wrap(term_category, width = 40),
      gene_ratio    = gene_count / bg_genes
    ) %>%
    mutate(term_category = factor(term_category, levels = rev(unique(term_category))))

  ggplot(df, aes(x = gene_ratio, y = term_category)) +
    geom_point(aes(size = gene_count, colour = -log10(FDR))) +
    scale_colour_viridis_c(name = expression(-log[10](adj.~p)),
                           option = "viridis",
                           guide = guide_colorbar(order = 2,
                                                  barheight = unit(10, "mm"),
                                                  barwidth  = unit(2,  "mm"),
                                                  ticks = FALSE,
                                                  title.position = "top")) +
    scale_size_continuous(
      name = "Gene count",
      breaks = c(5, 10, 15, 20, 25, 30),
      range = c(0.8, 3.5),
      guide  = guide_legend(order = 1,
                            keyheight = unit(0.28, "cm"),
                            keywidth  = unit(0.28, "cm"))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    labs(x = "Gene ratio", y = NULL) +
    paper_theme +
    theme(
      panel.grid.major.y = element_line(colour = "grey92",
                                        linewidth = 0.25),
      axis.text.y  = element_text(size = BASE_SIZE - 1),
      axis.title.x = element_text(size = BASE_SIZE),
      ## legend OUTSIDE plot, on the right
      legend.position    = "right",
      legend.justification = "top",
      legend.box         = "vertical",
      legend.box.spacing = unit(2, "pt"),
      legend.spacing.y   = unit(2, "pt"),
      legend.margin      = margin(0, 0, 0, 2)
    )
}

## ====================================================================
## Panel d — GO Biological Process dot plot
## Source: figures/data_tables/enrichment_BP.tsv
## ====================================================================
build_go_panel <- function() {
  bp <- fread("figures/data_tables/enrichment_BP.tsv") %>%
    slice_min(order_by = p.adjust, n = 23) %>%
    mutate(
      Description = stringr::str_wrap(Description, width = 40),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = Count / as.numeric(sub(".*/", "", GeneRatio))
    )

  ggplot(bp, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_viridis_c(name = expression(-log[10](adj.~p)),
                           option = "viridis",
                           guide = guide_colorbar(order = 2,
                                                  barheight = unit(10, "mm"),
                                                  barwidth  = unit(2,  "mm"),
                                                  ticks = FALSE,
                                                  title.position = "top")) +
    scale_size_continuous(
      name = "Gene count",
      breaks = c(5, 10, 15, 20, 25, 30),
      range = c(0.8, 3.5),
      guide  = guide_legend(order = 1,
                            keyheight = unit(0.28, "cm"),
                            keywidth  = unit(0.28, "cm"))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    labs(x = "Gene ratio", y = NULL) +
    paper_theme +
    theme(
      panel.grid.major.y = element_line(colour = "grey92",
                                        linewidth = 0.25),
      axis.text.y  = element_text(size = BASE_SIZE - 1),
      axis.title.x = element_text(size = BASE_SIZE),
      legend.position    = "right",
      legend.justification = "top",
      legend.box         = "vertical",
      legend.box.spacing = unit(2, "pt"),
      legend.spacing.y   = unit(2, "pt"),
      legend.margin      = margin(0, 0, 0, 2)
    )
}

## ====================================================================
## Build all four panels
## ====================================================================
cat("[..] building panel a — PCA\n")
pA <- build_pca_panel()

cat("[..] building panel b — VESTIB vs GTEx median expression\n")
pB <- build_expr_panel()

cat("[..] building panel c — disease / phenotype enrichment\n")
pC <- build_disease_panel()

cat("[..] building panel d — GO BP enrichment\n")
pD <- build_go_panel()

## ====================================================================
## Save each panel as its own PDF + PNG (for AI assembly)
## ====================================================================
save_fig <- function(p, stub, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stub, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(stub, ".png")), p,
         width = w, height = h, units = "mm", dpi = 600)
}

save_fig(pA, "panel_a_pca_vestib_gtex",   90, 90)
save_fig(pB, "panel_b_vestib_vs_gtex",    90, 90)
save_fig(pC, "panel_c_disease_enrich",    100, 100)
save_fig(pD, "panel_d_go_bp_enrich",      100, 100)

## ====================================================================
## Composite 2x2 — for a quick preview
## ====================================================================
fig3 <- ((pA | pB) / (pC | pD)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold",
                                size = BASE_SIZE + 3,
                                family = "sans"),
        plot.tag.position = c(0, 1))

ggsave(file.path(OUT_DIR, "figure3_combined.pdf"), fig3,
       width = 180, height = 180, units = "mm",
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "figure3_combined.png"), fig3,
       width = 180, height = 180, units = "mm",
       dpi = 600, limitsize = FALSE)

cat("[done] Figure 3 written to ", OUT_DIR, "\n",
    "    panel_a_pca_vestib_gtex.{pdf,png}   (90 x 90 mm)\n",
    "    panel_b_vestib_vs_gtex.{pdf,png}    (90 x 90 mm)\n",
    "    panel_c_disease_enrich.{pdf,png}    (90 x 90 mm)\n",
    "    panel_d_go_bp_enrich.{pdf,png}      (90 x 90 mm)\n",
    "    figure3_combined.{pdf,png}          (180 x 180 mm preview)\n",
    sep = "")



##
dim(pc_df_annot)


## ====================================================================
## EXTRA — panels c and d, faithful to gene_enrichment.R styling
## (theme_bw + viridis, original viridis colour scale, no label wrapping).
## Only change vs the original: slightly larger fonts. Saved with a
## "_orig_style" suffix so the redesigned versions above are kept.
##
## Run this section after the main script (or interactively). It does
## not overwrite gene_enrichment.R.
## ====================================================================

## ---- Panel d (GO BP) — original style with larger fonts ----------
bp_orig <- fread("figures/data_tables/enrichment_BP.tsv") %>%
  as.data.frame() %>%
  slice_min(order_by = p.adjust, n = 23) %>%
  mutate(
    Description = factor(Description, levels = rev(Description)),
    GeneRatio_num = Count / as.numeric(sub(".*/", "", GeneRatio))
  )

pD_orig <- ggplot(bp_orig, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_viridis_c(name = "-log10(adj. p)") +
  theme_bw() +
  labs(x = "Gene ratio", y = NULL) +
  theme(
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7)
  ) +
  scale_size_continuous(
    name = "Gene count",
    breaks = c(5, 10, 15, 20, 25, 30),
    guide  = guide_legend(keyheight = unit(0.4, "cm"),
                          keywidth  = unit(0.4, "cm"))
  )

ggsave(file.path(OUT_DIR, "panel_d_go_bp_enrich_orig_style.pdf"),
       pD_orig, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_d_go_bp_enrich_orig_style.png"),
       pD_orig, width = 5, height = 5, dpi = 600)

## ---- Panel c (Disease / Monarch / UniProt) — original style ------
gsea_cyto <- fread("src/new_enrichment_all_databases.csv")
gsea_plot_orig <- as.data.frame(gsea_cyto) %>%
  dplyr::rename(
    bg_genes      = `background genes`,
    gene_count    = genes,
    category_name = category,
    term          = description,
    FDR           = `FDR value`,
    pval          = `p-value`
  ) %>%
  dplyr::filter(!category_name %in% c(
    "GO Biological Process", "COMPARTMENTS",
    "GO Cellular Component", "GO Molecular Function",
    "SMART Domains", "WikiPathways", "STRING Clusters"
  )) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  mutate(
    term          = factor(term, levels = rev(unique(term))),
    category_name = sub(" .*", "", category_name),
    term_category = paste0(term, " (", category_name, ")"),
    gene_ratio    = gene_count / bg_genes
  )

pC_orig <- ggplot(gsea_plot_orig, aes(x = gene_ratio, y = term_category)) +
  geom_point(aes(size = gene_count, color = -log10(FDR))) +
  scale_color_viridis_c(
    name  = "-log10(adj. p)",
    guide = guide_colorbar(order = 2)
  ) +
  scale_size_continuous(
    name   = "Gene count",
    breaks = c(5, 10, 15, 20, 25, 30),
    guide  = guide_legend(order = 1,
                          keyheight = unit(0.3, "cm"),
                          keywidth  = unit(0.3, "cm"))
  ) +
  theme_bw() +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7)
  )

ggsave(file.path(OUT_DIR, "panel_c_disease_enrich_orig_style.pdf"),
       pC_orig, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_c_disease_enrich_orig_style.png"),
       pC_orig, width = 5, height = 5, dpi = 600)

cat("[done] orig-style panels written:\n",
    "    panel_c_disease_enrich_orig_style.{pdf,png}   (5 x 5 in)\n",
    "    panel_d_go_bp_enrich_orig_style.{pdf,png}     (5 x 5 in)\n",
    sep = "")


## ====================================================================
## EXTRA — panel a with tissue labels (faithful to PCA_VESTIB_vs_GTEx.R)
## Labels only the tissues listed in the original script:
##   VESTIB, BRN_TISSUES, PTTARY, WHLBLD, LCL, SPLEEN, GEUVADIS,
##   PANCREAS, ADRENAL, THYROID
## PANCREAS/ADRENAL are not in the data (codes are PNCREAS/ADRNLG) and
## simply won't appear, mirroring the original script's behaviour.
##
## Saved as panel_a_pca_vestib_gtex_labeled.{pdf,png} — does NOT
## overwrite panel_a_pca_vestib_gtex.{pdf,png}.
## ====================================================================

label_tissues <- c("VESTIB", "BRN_TISSUES", "PTTARY", "WHLBLD",
                   "LCL", "SPLEEN", "GEUVADIS",
                   "PANCREAS", "ADRENAL", "THYROID")

## Re-derive pc_df_annot (mirrors build_pca_panel) -------------------
load("figures/data_tables/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData")
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
pc_df <- data.frame(human_id = rownames(pca$x),
                    PC1 = pca$x[, 1],
                    PC2 = pca$x[, 2],
                    stringsAsFactors = FALSE)
annot <- build_annot_gtex_vest()
pc_df_annot <- pc_df %>%
  inner_join(annot, by = "human_id") %>%
  mutate(tissue = dplyr::if_else(grepl("^BRN", tissue), "BRN_TISSUES", tissue)) %>%
  dplyr::filter(!human_id %in% c("31907773", "32073088")) %>%
  dplyr::filter(!(tissue == "VESTIB" & PC1 > 0))

gtex_palette <- load_gtex_palette()
gtex_palette["BRN_TISSUES"] <- "#EEEE00"
gtex_palette["VESTIB"]      <- pal$red
gtex_palette["GEUVADIS"]    <- "grey30"

## Label positions: median of each tissue's cluster + manual nudges ---
label_df <- pc_df_annot %>%
  dplyr::filter(tissue %in% label_tissues) %>%
  group_by(tissue) %>%
  summarise(PC1 = median(PC1, na.rm = TRUE),
            PC2 = median(PC2, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    nudge_x = case_when(
      tissue == "VESTIB"      ~  0.5,
      tissue == "BRN_TISSUES" ~ -0.5,
      TRUE                    ~  0
    ),
    nudge_y = case_when(
      tissue == "WHLBLD" ~  0.5,
      tissue == "LCL"    ~ -0.5,
      TRUE               ~  0
    )
  )

pA_lab <- ggplot(pc_df_annot, aes(PC1, PC2, colour = tissue)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             linewidth = 0.25, colour = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.25, colour = "grey60") +
  geom_point(size = 0.45, alpha = 0.85) +
  geom_point(data = dplyr::filter(pc_df_annot, tissue == "VESTIB"),
             size = 0.75, alpha = 0.95, colour = pal$red) +
  geom_label(
    data = label_df,
    aes(x = PC1 + nudge_x, y = PC2 + nudge_y, label = tissue),
    inherit.aes = FALSE,
    colour = "black", fill = "white",
    label.size = 0.10, label.padding = unit(0.18, "lines"),
    size = 1.8
  ) +
  scale_colour_manual(values = gtex_palette, guide = "none") +
  labs(x = sprintf("PC1 (%.1f%%)", var_expl[1] * 100),
       y = sprintf("PC2 (%.1f%%)", var_expl[2] * 100)) +
  paper_theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave(file.path(OUT_DIR, "panel_a_pca_vestib_gtex_labeled.pdf"),
       pA_lab, width = 90, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_a_pca_vestib_gtex_labeled.png"),
       pA_lab, width = 90, height = 90, units = "mm", dpi = 600)

cat("[done] labeled panel a written:\n",
    "    panel_a_pca_vestib_gtex_labeled.{pdf,png}   (90 x 90 mm)\n",
    sep = "")
