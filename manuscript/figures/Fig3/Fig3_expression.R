## ====================================================================
## Figure 3: otolith-organ transcriptome vs GTEx and vestibular-enriched genes (7.2 in wide)
##   a. PCA of otolith organs + GTEx tissues (Pantry expression, same pipeline), tissue groups labelled
##   b. Otolith median vs highest GTEx median expression per gene; 659 vestibular-enriched genes (delta log2 > 2) in red,
##      inner-ear genes of panel d outlined
##   c. Canonical inner-ear marker genes: otolith organs vs the highest GTEx tissue
##   d. Expression of the genes driving the hearing / inner-ear enrichments across tissues (median log2(TPM+1);
##      GTEx tissues only) + the enriched terms each gene belongs to
##   e. GO Biological Process (clusterProfiler) and disease / phenotype (STRING: Disease, Monarch, UniProt keywords)
##      enrichment of the enriched genes; redundant terms (gene-set Jaccard >= 0.7) dropped; bars = -log10(FDR)
## Inputs: analyses/08_expression_gtex/results (cached tables in _fig3_cache/; the per-tissue median TPM table is
## computed once and cached there, median_tpm_by_tissue.rds).
## Output: manuscript/figures/Fig3/ (Fig3.pdf/png, Fig3_5.8x5.8in.pdf/png, panels/).
## Run from the project root on a compute node (env r_plots):
##   Rscript manuscript/figures/Fig3/Fig3_expression.R
## ====================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(ggrepel); library(cowplot)
  library(matrixStats); library(scales); library(patchwork); library(jsonlite)
})
cfg <- yaml::read_yaml("config/config.yaml")   # + lab-specific values (config_lab.yaml, not in git)
if (file.exists("config/config_lab.yaml")) cfg <- modifyList(cfg, yaml::read_yaml("config/config_lab.yaml"))
OUT   <- "manuscript/figures/Fig3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CACHE <- "analyses/08_expression_gtex/results/_fig3_cache"
FONT  <- "Nimbus Sans"
W <- 7.2; H <- 7.9
RED <- "#C8102E"; GREY <- "grey75"

th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        legend.key.size = unit(7, "pt"), plot.margin = margin(4, 6, 2, 3))

cache_or_compute <- function(name, builder) {
  fp <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(fp)) return(readRDS(fp))
  cat("[build]", name, "\n"); obj <- builder(); saveRDS(obj, fp); obj
}

## ---- gene annotation (Ensembl 110, as the expression quantification)
gene_info <- cache_or_compute("gene_info_ens110", function() {
  g <- rtracklayer::readGFF("resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf",
                            tags = c("gene_id", "gene_name", "gene_biotype"), filter = list(type = "gene"))
  as.data.table(g)[, .(gene_id, gene_name, gene_biotype)]
})
sym <- setNames(gene_info$gene_name, gene_info$gene_id)

## ---- per-tissue median TPM (protein-coding + lncRNA), same samples/tissues as Fig3_expression.R panel b
gtex_root <- cfg$refs$gtex_pantry_phenos   # GTEx tissues run through Pantry
med <- cache_or_compute("median_tpm_by_tissue", function() {
  keep <- gene_info[gene_biotype %in% c("protein_coding", "lncRNA")]$gene_id
  v <- fread("analyses/05_pantry/phenotyping/output/vestibular/unnorm/expression.bed")
  v <- v[phenotype_id %in% keep, !unlist(cfg$samples$outlier_samples), with = FALSE]
  m <- data.table(gene_id = v$phenotype_id, VESTIB = rowMedians(as.matrix(v[, -(1:4)])))
  tissues <- setdiff(list.dirs(gtex_root, recursive = FALSE, full.names = FALSE),
                     c("ref", "TESTIS", "BLDDER", "CVSEND", "CVXECT", "FLLPNT", "KDNMDL"))
  tissues <- tissues[file.exists(file.path(gtex_root, tissues, "output/unnorm/expression.bed"))]
  for (tt in tissues) {
    cat("  median:", tt, "\n")
    e <- fread(file.path(gtex_root, tt, "output/unnorm/expression.bed"))[phenotype_id %in% keep]
    m[, (tt) := rowMedians(as.matrix(e[, -(1:4)]))[match(gene_id, e$phenotype_id)]]
  }
  m
})

## ---- panel b data (cached by Fig3_expression.R): genes with padj < 0.01; enriched = log2FC > 2 (659 genes)
dt <- readRDS(file.path(CACHE, "dt_plot_panel_b.rds"))
dt[, `:=`(gene_name = sym[gene_id], enriched = log2FC > 2)]
n_enr <- sum(dt$enriched)
enr_names <- dt[enriched == TRUE]$gene_name

## ---- enrichment tables
split_genes <- function(x, sep) strsplit(x, sep, fixed = TRUE)
bp <- fread("analyses/08_expression_gtex/results/enrichment_BP.tsv")[order(p.adjust)]
bp[, `:=`(gl = split_genes(geneID, "/"), FDR = p.adjust, term = Description, source = "GO BP")]
cy <- fread("analyses/08_expression_gtex/results/enrichment_stringdb_cytoscape.csv")
dis <- cy[category %in% c("Disease", "Monarch Phenotype", "UniProt Keywords") & `FDR value` < 0.05][order(`FDR value`)]
dis[, `:=`(gl = split_genes(gene_names, "|"), FDR = `FDR value`, Count = lengths(split_genes(gene_names, "|")),
           term = description, source = c(Disease = "DO", "Monarch Phenotype" = "Monarch", "UniProt Keywords" = "UniProt")[category])]

## drop terms whose gene set overlaps a more significant kept term (Jaccard >= 0.7)
nonredundant <- function(d, jac = 0.7) {
  keep <- logical(nrow(d))
  for (i in seq_len(nrow(d))) {
    keep[i] <- !any(vapply(which(keep), function(j) {
      a <- d$gl[[i]]; b <- d$gl[[j]]; length(intersect(a, b)) / length(union(a, b)) >= jac }, logical(1)))
  }
  d[keep]
}
theme_of <- function(term, gl) {
  hist <- vapply(gl, function(g) mean(grepl("^H[1-4]", g)), numeric(1))
  fifelse(hist > 0.5, "Histone / chromatin",
          fifelse(grepl("sound|sens|ear|deaf|hearing|vestib|mechanical|stimulus", term, ignore.case = TRUE),
                  "Sensory / inner ear", "Other"))
}
THEMES <- c("Sensory / inner ear", "Histone / chromatin", "Other")
bp_k  <- head(nonredundant(bp), 10)[, theme := factor(theme_of(term, gl), levels = THEMES)]
dis_k <- head(nonredundant(dis), 10)[, theme := factor(theme_of(term, gl), levels = THEMES)]
## first letter upper case (GO terms are lower case)
cap1 <- function(s) paste0(toupper(substr(s, 1, 1)), substring(s, 2))


cowplot::set_null_device("cairo")   # text measuring device (cairo) for the composite
sq <- theme(aspect.ratio = 1)
tagth <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT), plot.tag.position = c(0, 1))

## ====================================================================
## a — PCA, plain labels
## ====================================================================
load("analyses/08_expression_gtex/results/gtex_vest_logCPM_vestib_gtex_same_pipeline_PCs.RData")   # -> pca
ve <- pca$sdev^2 / sum(pca$sdev^2)
annot <- readRDS(file.path(CACHE, "annot_gtex_vest.rds"))
pcd <- data.table(human_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2]) %>%
  inner_join(annot, by = "human_id") %>%
  mutate(tissue = if_else(grepl("^BRN", tissue), "BRN_TISSUES", tissue)) %>%
  filter(!human_id %in% unlist(cfg$samples$outlier_samples_gtex_pca), !(tissue == "VESTIB" & PC1 > 0)) %>% as.data.table()
gpal <- readRDS(file.path(CACHE, "gtex_palette.rds"))
gpal["BRN_TISSUES"] <- "#D9D900"; gpal["GEUVADIS"] <- "grey40"; gpal["VESTIB"] <- RED
grp_lab <- c(VESTIB = "Otolith\norgans", BRN_TISSUES = "Brain", PTTARY = "Pituitary", WHLBLD = "Whole\nblood",
             LCL = "GTEx LCLs", GEUVADIS = "GEUVADIS LCLs", SPLEEN = "Spleen", MSCLSK = "Muscle", NERVET = "Nerve")
cent <- pcd[tissue %in% names(grp_lab), .(PC1 = median(PC1), PC2 = median(PC2)), by = tissue][, lab := grp_lab[tissue]]
rx <- diff(range(pcd$PC1)); ry <- diff(range(pcd$PC2))
## labels: black text placed by hand in empty space next to each group (u, v = fraction of the PC1 / PC2 data range);
## a short leader line goes to the core sample (closest 40% to the group median) nearest the label
lpos <- data.table(tissue = c("BRN_TISSUES", "PTTARY", "VESTIB", "SPLEEN", "WHLBLD", "LCL", "GEUVADIS", "MSCLSK", "NERVET"),
                   u      = c(0.15,          0.17,     0.42,     0.81,     0.68,     0.70,  0.66,       0.78,     0.30),
                   v      = c(0.74,          0.29,     0.52,     0.42,     0.76,     0.87,  0.99,       0.12,     0.04))
x0 <- min(pcd$PC1); y0 <- min(pcd$PC2)
lpos[, `:=`(lx = x0 + u * rx, ly = y0 + v * ry, lab = grp_lab[tissue])]
anch <- lpos[, {
  g <- pcd[tissue == .BY$tissue]
  dc <- ((g$PC1 - median(g$PC1)) / rx)^2 + ((g$PC2 - median(g$PC2)) / ry)^2
  g <- g[dc <= quantile(dc, 0.4)]                      # dense core only, so lines never end on stray samples
  k <- which.min(((g$PC1 - lx) / rx)^2 + ((g$PC2 - ly) / ry)^2)
  .(PC1 = g$PC1[k], PC2 = g$PC2[k], lx, ly, lab)
}, by = tissue]
pA <- ggplot(pcd[tissue != "VESTIB"], aes(PC1, PC2)) +
  geom_point(aes(colour = tissue), size = 0.2, alpha = 0.5, stroke = 0) +
  geom_point(data = pcd[tissue == "VESTIB"], colour = RED, size = 0.55, alpha = 0.9, stroke = 0) +
  geom_text_repel(data = anch, aes(label = lab), nudge_x = anch$lx - anch$PC1, nudge_y = anch$ly - anch$PC2,
                  force = 0, force_pull = 0, size = 4.8 / .pt, family = FONT, lineheight = 0.85, colour = "grey10",
                  segment.colour = "grey25", segment.size = 0.25, min.segment.length = 0, box.padding = 0.08,
                  point.padding = 0, max.overlaps = Inf, seed = 7) +
  scale_colour_manual(values = gpal, guide = "none") +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2]), tag = "a") + th + tagth +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ====================================================================
## inner-ear enrichment genes (heatmap d; outlined in b)
## ====================================================================
memb <- rbind(
  data.table(set = "Sound perception",      gene = bp[term == "sensory perception of sound"]$gl[[1]]),
  data.table(set = "Inner ear development", gene = bp[term == "inner ear development"]$gl[[1]]),
  data.table(set = "Deafness",    gene = unique(unlist(dis[term %in% c("Nonsyndromic deafness", "Deafness")]$gl))),
  data.table(set = "Vestibular dysfunction", gene = unlist(dis[term == "Vestibular dysfunction"]$gl)))
memb <- memb[gene %in% enr_names]
hg <- unique(memb$gene)

## ====================================================================
## b — otolith vs highest reference median (equal axes, no labels)
## ====================================================================
dt[, cls := factor(fifelse(gene_name %in% hg & enriched, "ear", fifelse(enriched, "enr", "other")), levels = c("other", "enr", "ear"))]
lab_cls <- c(other = "Other genes", enr = sprintf("Vestibular-enriched (n = %d)", n_enr), ear = sprintf("Inner-ear genes in d (n = %d)", length(hg)))
pB <- ggplot(dt[order(cls)], aes(max_other, VESTIB)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 2, linetype = 2, colour = RED, linewidth = 0.35) +
  geom_point(aes(colour = cls, size = cls, fill = cls, shape = cls), stroke = 0.3) +
  scale_colour_manual(values = c(other = alpha("grey55", 0.35), enr = alpha(RED, 0.8), ear = "black"), labels = lab_cls, name = NULL) +
  scale_fill_manual(values = c(other = "grey55", enr = RED, ear = RED), labels = lab_cls, name = NULL) +
  scale_size_manual(values = c(other = 0.15, enr = 0.35, ear = 1.1), labels = lab_cls, name = NULL) +
  scale_shape_manual(values = c(other = 16, enr = 16, ear = 21), labels = lab_cls, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = c(1, 1, 1.4), alpha = 1))) +
  annotate("text", x = 15.6, y = 14.7, label = "1:1", size = 5.5 / .pt, family = FONT, colour = "grey40", hjust = 1) +
  annotate("text", x = 12.2, y = 15.2, label = "+2", size = 5.5 / .pt, family = FONT, colour = RED, hjust = 1) +
  scale_x_continuous(breaks = seq(0, 16, 4), limits = c(0, 16.6), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 16, 4), limits = c(0, 16.6), expand = c(0, 0)) +
  labs(x = expression("Highest reference median, " * log[2] * "(TPM+1)"),
       y = expression("Otolith median, " * log[2] * "(TPM+1)"), tag = "b") + th + tagth +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0), legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5.5), legend.spacing.y = unit(0, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA), legend.margin = margin(1, 1, 1, 1))

## ====================================================================
## c — canonical inner-ear marker genes: otolith vs highest GTEx tissue
## ====================================================================
gx <- fromJSON("resources/gtex_colors.json")
tname <- setNames(names(gx), vapply(gx, function(v) v$tissue_abbrv, ""))
short <- function(a) {
  n <- ifelse(a %in% names(tname), tname[a], a)
  n <- sub(" \\(.*\\)$", "", n)
  n <- sub("Cells - EBV-transformed lymphocytes", "LCL", n); n <- sub("Cells - Cultured fibroblasts", "Fibroblasts", n)
  n <- sub("^Brain - ", "Brain: ", n); n <- sub("^Nerve - Tibial", "Tibial nerve", n)
  n <- sub("^Skin - Not Sun Exposed", "Skin (not sun-exp.)", n); n <- sub("^Skin - Sun Exposed", "Skin (sun-exp.)", n)
  n <- sub("^Esophagus - ", "Esophagus: ", n); n <- sub("^Artery - ", "Artery: ", n); n <- sub("^Adipose - ", "Adipose: ", n)
  n <- sub("^Colon - ", "Colon: ", n); n <- sub("^Heart - ", "Heart: ", n); sub("^Kidney - ", "Kidney: ", n)
}
## c = canonical inner-ear marker genes (chosen a priori by category, not by result): otolith median vs the highest
## GTEx tissue median (GEUVADIS excluded)
GTEX <- setdiff(colnames(med), c("gene_id", "VESTIB", "GEUVADIS"))
mk <- rbind(data.table(cat = "Otolithic membrane", gene = c("OC90", "OTOL1", "OTOG", "OTOGL", "OTOA", "TECTA")),
            data.table(cat = "Hair cell", gene = c("POU4F3", "MYO7A", "USH1C", "OTOF", "TMC1", "LHFPL5", "PCDH15", "CDH23")))
mm <- med[match(gene_info[match(mk$gene, gene_name)]$gene_id, gene_id)]
gx_m <- as.matrix(mm[, ..GTEX])
mk[, `:=`(oto = log2(mm$VESTIB + 1), gtx = log2(rowMaxs(gx_m, na.rm = TRUE) + 1),
          gtx_t = GTEX[max.col(replace(gx_m, is.na(gx_m), -1), ties.method = "first")])]
## single panel (no facets, so row-1 axes stay aligned): bold group headings above each gene group
mk[, cat := factor(cat, levels = unique(cat))]
setorder(mk, cat, -oto)
mk[, oto_gap := oto - gtx]
setorder(mk, cat, -oto_gap)
ypos <- c(); hdr <- data.table(); y <- 0
for (k in levels(mk$cat)) {
  y <- y - 1; hdr <- rbind(hdr, data.table(cat = k, y = y))
  n <- sum(mk$cat == k); ypos <- c(ypos, y - seq_len(n)); y <- y - n - 0.4
}
mk[, y := ypos]
ml <- melt(mk, id.vars = c("cat", "gene", "y"), measure.vars = c("oto", "gtx"), variable.name = "src", value.name = "x")
pCt <- ggplot(mk, aes(y = y)) +
  geom_hline(yintercept = mk$y, colour = "grey94", linewidth = 0.2) +
  geom_segment(aes(x = gtx, xend = oto, yend = y), colour = "grey70", linewidth = 0.6) +
  geom_point(data = ml, aes(x, y, fill = src), shape = 21, size = 1.8, colour = "grey15", stroke = 0.25) +
  geom_text(data = hdr, aes(x = 0, y = y, label = cat), hjust = 0, vjust = 0.5, size = 6 / .pt, family = FONT,
            fontface = "bold") +
  scale_y_continuous(breaks = mk$y, labels = as.character(mk$gene), expand = expansion(add = c(0.6, 0.3))) +
  scale_fill_manual(values = c(oto = RED, gtx = "grey80"), labels = c(oto = "Otolith organs", gtx = "Highest GTEx tissue"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 12, 3), limits = c(0, NA), expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = expression("Median " * log[2] * "(TPM+1)"), y = NULL, tag = "c") + th + tagth +
  theme(axis.text.y = element_text(face = "italic", size = 5.8),
        legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, -4, 0),
        legend.key.size = unit(6, "pt"))

## ====================================================================
## d — heatmap (as v2 c)
## ====================================================================
hm <- med[match(hg, sym[gene_id])]
mat <- log2(as.matrix(hm[, -1]) + 1); rownames(mat) <- hg
gt <- setdiff(colnames(mat), c("VESTIB", "GEUVADIS"))
mat <- mat[, c("VESTIB", gt)]
col_ord <- c("VESTIB", gt[hclust(dist(t(mat[, gt])), "average")$order])
row_ord <- rownames(mat)[order(mat[, "VESTIB"] - rowMaxs(mat[, gt]), decreasing = FALSE)]
hl <- as.data.table(mat, keep.rownames = "gene") %>% melt(id.vars = "gene", variable.name = "tissue", value.name = "expr")
GRP <- c("Otolith", "GTEx tissues", "Enriched term")
hl[, `:=`(gene = factor(gene, levels = row_ord), x = factor(tissue, levels = c(col_ord, levels(factor(memb$set)))),
          grp = factor(fifelse(tissue == "VESTIB", GRP[1], GRP[2]), levels = GRP))]
mb <- memb[, .(gene = factor(gene, levels = row_ord), x = factor(set, levels = levels(hl$x)), grp = factor(GRP[3], levels = GRP))]
mb_bg <- CJ(gene = factor(row_ord, levels = row_ord), x = factor(unique(memb$set), levels = levels(hl$x)))[, grp := factor(GRP[3], levels = GRP)]
pD <- ggplot(hl, aes(x, gene)) + geom_tile(aes(fill = expr), colour = "white", linewidth = 0.1) +
  geom_tile(data = mb_bg, fill = "grey96", colour = "white", linewidth = 0.3) +
  geom_point(data = mb, size = 1.1, colour = "grey15") +
  facet_grid(~grp, scales = "free_x", space = "free_x",
             labeller = as_labeller(c(Otolith = "", "GTEx tissues" = "GTEx tissues", "Enriched term" = "Term"))) +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = expression(atop("Median", log[2] * "(TPM+1)")),
                       guide = guide_colourbar(barheight = unit(36, "pt"), barwidth = unit(4, "pt"))) +
  scale_x_discrete(labels = function(x) ifelse(x == "VESTIB", "Otolith organs", x), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = "d") + th + tagth +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.8),
        axis.text.y = element_text(face = "italic", size = 5.5), axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.2), strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold"), panel.spacing.x = unit(3, "pt"),
        legend.position = "right", plot.margin = margin(2, 2, 2, 9))

## ====================================================================
## e — enrichment: GO BP and disease/phenotype in one lollipop panel
## ====================================================================
thcol <- c("Sensory / inner ear" = RED, "Histone / chromatin" = "#4A6FA5", "Other" = "grey55")
en <- rbind(bp_k[, .(db = "GO biological process", lab = cap1(term), FDR, Count, theme)],
            dis_k[, .(db = "Disease / phenotype", lab = sprintf("%s (%s)", term, source), FDR, Count, theme)])
en[db == "GO biological process", lab := stringr::str_wrap(lab, width = 30)]   # wrap the two long GO terms (left plot gets more width)
en[, `:=`(db = factor(db, levels = unique(db)), y = -log10(FDR))]
en[, lab := factor(paste(db, lab, sep = "||"), levels = rev(paste(db, lab, sep = "||")))]
pE <- ggplot(en, aes(y, lab)) +
  geom_segment(aes(x = 0, xend = y, yend = lab, colour = theme), linewidth = 0.5) +
  geom_point(aes(size = Count, fill = theme), shape = 21, colour = "grey15", stroke = 0.2) +
  facet_wrap(~db, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub(".*\\|\\|", "", x)) +
  scale_colour_manual(values = thcol, guide = "none") +
  scale_fill_manual(values = thcol, name = "Term theme", drop = FALSE,
                    guide = guide_legend(override.aes = list(size = 2), order = 1)) +
  scale_size_area(max_size = 3.2, breaks = c(5, 10, 20, 30), name = "Genes", guide = guide_legend(order = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
  labs(x = expression(-log[10] * FDR), y = NULL, tag = "e") + th + tagth +
  theme(axis.text.y = element_text(size = 5.8, lineheight = 0.85), strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", hjust = 0.5), strip.clip = "off", panel.spacing.x = unit(12, "pt"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.2),
        legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(10, "pt"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 18, 2, 9))

## ====================================================================
## assemble (patchwork within rows: equal square panels, aligned x axes)
## ====================================================================
row1 <- (pA | pB | pCt) + plot_layout(widths = c(1, 1, 0.95)) & theme(plot.margin = margin(4, 8, 2, 5))
fig3 <- plot_grid(row1, pD, pE, ncol = 1, rel_heights = c(2.5, 3.05, 2.35))
ggsave(file.path(OUT, "Fig3.pdf"), fig3, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT, "Fig3.png"), fig3, width = W, height = H, dpi = 600, bg = "white")
## square version for the manuscript (Word: 5.8 x 5.8 in): laid out at 7.2 x 7.2 in (rows a bit shorter), then scaled
## to 5.8 in with ghostscript (vector) + 600 dpi PNG with the DPI tag, so Word places it at exactly 5.8 x 5.8 in
fig3_sq <- plot_grid(row1, pD, pE, ncol = 1, rel_heights = c(2.35, 2.75, 2.1))
sq_src <- file.path(OUT, "panels", "Fig3_square_7.2in.pdf")   # intermediate
ggsave(sq_src, fig3_sq, width = 7.2, height = 7.2, device = cairo_pdf)
sq_out <- file.path(OUT, "Fig3_5.8x5.8in")
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=pdfwrite", "-dFIXEDMEDIA", "-dPDFFitPage",
                "-dAutoRotatePages=/None", "-dEmbedAllFonts=true", "-dDEVICEWIDTHPOINTS=417.6", "-dDEVICEHEIGHTPOINTS=417.6",
                paste0("-sOutputFile=", sq_out, ".pdf"), sq_src))
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=png16m", "-r600", "-dTextAlphaBits=4",
                "-dGraphicsAlphaBits=4", paste0("-sOutputFile=", sq_out, ".png"), paste0(sq_out, ".pdf")))
for (p in list(list(pA, "Fig3a", 2.5, 2.5), list(pB, "Fig3b", 2.5, 2.5), list(pCt, "Fig3c", 2.5, 2.5),
               list(pD, "Fig3d", W, 3.05), list(pE, "Fig3e", W, 2.2)))
  ggsave(file.path(OUT, "panels", paste0(p[[2]], ".pdf")), p[[1]], width = p[[3]], height = p[[4]], device = cairo_pdf)
cat("markers:\n"); print(mk[, .(cat, gene, oto = round(oto, 2), gtx = round(gtx, 2), gtx_t)])
cat(sprintf("enriched genes: %d; heatmap genes: %d; reference columns: %d\n", n_enr, length(hg), length(gt)))
cat("[done]\n")


### END