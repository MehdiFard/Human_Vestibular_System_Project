## Fig 4: multimodal RNA phenotypes, xGenes and xTWAS genes in the otolith organs.
##   a  phenotypes per modality (protein-coding / lncRNA)        (Pantry paper figure1.R, panel B)
##   b  xGenes per modality (cross-modality mapping)              (figure2abc.R, panel A, version 2)
##   c  genes with xTWAS hits per GWAS trait, by modality of the gene's top hit (20 traits with most genes)
##                                                                (figure4.R, panel A right)
## Inputs: analyses/05_pantry/results/xqtl_summary (phenotypes, xQTLs), analyses/06_twas (all TWAS tests, GWAS sizes).
## Writes panels/ (one PDF per panel), source_data/ (data behind each panel) and Fig4.pdf / Fig4.png.
## Run from the project root:  Rscript manuscript/figures/Fig4/fig4.R
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
set.seed(1)   # ties in the top TWAS hit per trait and gene are broken at random (slice_sample)
X <- "analyses/05_pantry/results/xqtl_summary/vestibular"
TW <- "analyses/06_twas"
O <- "manuscript/figures/Fig4"
for (d in c("panels", "source_data")) dir.create(file.path(O, d), recursive = TRUE, showWarnings = FALSE)

modalities <- c(expression = "Expression", isoforms = "Isoform ratio", splicing = "Intron excision",
                alt_TSS = "Alt. TSS", alt_polyA = "Alt. polyA", stability = "RNA stability", RNA_editing = "RNA editing")
modality_colors <- c(Expression = "#e41a1c", `Isoform ratio` = "#377eb8", `Intron excision` = "#4daf4a",
                     `Alt. TSS` = "#984ea3", `Alt. polyA` = "#ff7f00", `RNA stability` = "#fdc11c", `RNA editing` = "#984210")
save_panel <- function(p, stem, w, h) ggsave(file.path(O, "panels", paste0(stem, ".pdf")), p, width = w, height = h)

## ---- a: phenotypes per modality ----------------------------------------------------------------------------------
phenos <- read_tsv(paste0(X, ".phenotypes.tsv.gz"), col_types = "cccc") |>
    mutate(modality = factor(modalities[modality], levels = modalities))
pos1 <- with(phenos, sum(modality == modalities[1] & gene_biotype == "protein_coding") / 2 / 1000)
pos2 <- with(phenos, (sum(modality == modalities[1]) - (sum(modality == modalities[1] & gene_biotype == "lincRNA") / 2)) / 1000)
data_a <- phenos |>
    dplyr::count(modality, gene_biotype) |>
    mutate(thousands = n / 1000, modality = fct_rev(modality))
p_a <- ggplot(data_a, aes(x = thousands, y = modality, fill = modality, alpha = gene_biotype)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    annotate("segment", x = pos1, xend = pos1, y = 7.4, yend = 7.7, linewidth = 0.3) +
    annotate("text", x = pos1 - 3, y = 7.75, hjust = 0, vjust = 0, size = 1.5, label = "protein-coding") +
    annotate("segment", x = pos2, xend = pos2, y = 7.4, yend = 7.5, linewidth = 0.3) +
    annotate("text", x = pos2 - 2, y = 7.55, hjust = 0, vjust = 0, size = 1.5, label = "lncRNA") +
    scale_fill_manual(values = modality_colors) +
    scale_alpha_manual(values = c(0.5, 1)) +
    expand_limits(y = 8) +
    xlab("Phenotypes (×1000)") +
    ylab("Modality") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.title.y = element_text(margin = margin(r = 10)))   # keep "Modality" clear of the tick labels
save_panel(p_a, "Fig4a_phenotypes_per_modality", 4.16, 2.08)
write_tsv(data_a, file.path(O, "source_data", "Fig4a_phenotypes_per_modality.tsv"))

## ---- b: xGenes per modality --------------------------------------------------------------------------------------
qtls <- read_tsv(paste0(X, ".comb.qtls.tsv"), col_types = "cicccd") |>
    mutate(modality = factor(modalities[modality], levels = modalities))
data_b <- qtls |>
    dplyr::count(modality, gene_id, name = "n_QTLs") |>
    mutate(n_QTLs = as_factor(n_QTLs) |> fct_other(keep = as.character(1:4), other_level = "5+"),
           modality = fct_infreq(modality)) |>
    group_by(modality, n_QTLs) |>
    summarise(genes = n(), .groups = "drop")
p_b <- data_b |>
    group_by(modality) |>
    summarise(genes = sum(genes), .groups = "drop") |>
    ggplot(aes(x = genes, y = modality, fill = modality)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = genes), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values = modality_colors) +
    xlab("No. of xGenes") +
    ylab("Modality") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), legend.position = "none")
save_panel(p_b, "Fig4b_xgenes_per_modality", 4.16, 2.08)
write_tsv(data_b, file.path(O, "source_data", "Fig4b_xgenes_per_modality_by_n_xqtls.tsv"))

## ---- c: genes with xTWAS hits per trait, by modality of the top hit ------------------------------------------------
modalities_c <- replace(modalities, "splicing", "Intron excision ratio")
colors_c <- setNames(modality_colors, replace(names(modality_colors), 3, "Intron excision ratio"))
twas <- read_tsv(file.path(TW, "results/vestibular/output/twas_all_tests.tsv.gz"), col_types = "ccccdcddddddd") |>
    mutate(modality = factor(modalities_c[modality], levels = modalities_c)) |>
    filter(TWAS.P < 5e-8 / 7)
ssize <- read_tsv(file.path(TW, "gwas_metadata.tsv"), col_types = cols(trait = "c", sample_size = "i", .default = "-")) |>
    mutate(sample_size = sample_size |> round(-3) |> str_sub(1, -4) |> str_c("K")) |>
    deframe()
trait_labels <- c(
    anx_eur = "Anxiety", balance_meta_eur = "Balance meta eur", balance_meta1 = "Balance meta",
    cognitive_skills_ext = "Cognitive skills", fluid_intelligence_eur = "Fluid intelligence", tinnitus_eur = "Tinnitus",
    frailty_gwas = "Frailty", alzheimer = "Alzheimer", parkinson_gp2_eur = "Parkinson", hearing_hleur = "Hearing",
    icp_fa = "ICP FA", icp_md = "ICP MD", mcp_fa = "MCP FA", mcp_md = "MCP MD",
    non_cognitive_skills_ext = "Non-cognitive skills", depression_mdd = "Depression", scp_fa = "SCP FA", scp_md = "SCP MD",
    cerebellar_volume = "Cerebellar volume", cerebral_volume = "Cerebral volume", subcortical_volume = "Subcortical volume",
    ukbb_falls = "Falls (UKBB)", MVP_hearing_loss = "Hearing loss (MVP)")
trait_labels <- str_glue("{trait_labels} (n={ssize[names(trait_labels)]})") |> purrr::set_names(names(trait_labels))

## genes with expression-only / expression + other / non-expression-only hits (source data)
hit_types <- twas |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(trait, gene_id, modality_type) |>
    summarise(modality_hits = str_glue("{'Expression' %in% modality_type}_{'Other' %in% modality_type}"),
              .by = c(trait, gene_id)) |>
    mutate(modality_hits = c(TRUE_FALSE = "Expression only", TRUE_TRUE = "Expr. & other(s)",
                             FALSE_TRUE = "Non-expr. only")[modality_hits] |>
               fct_relevel("Expression only", "Expr. & other(s)", "Non-expr. only"))
hit_types20 <- hit_types |> filter(trait %in% levels(fct_infreq(trait))[1:20])

## top hit per trait and gene (ties broken at random, as in figure4.R; seeded above)
top_hits <- twas |>
    group_by(trait, gene_id) |>
    slice_min(TWAS.P, n = 1, with_ties = TRUE) |>
    slice_sample(n = 1) |>
    ungroup() |>
    filter(trait %in% levels(fct_infreq(trait))[1:20])
stopifnot(all(top_hits$trait %in% names(trait_labels)))
p_c <- top_hits |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality)) +
    geom_bar(width = 0.8) +
    scale_fill_manual(values = colors_c) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.justification = c(0, 1),
          legend.key.size = unit(10, "pt"),
          legend.position = "inside",
          legend.position.inside = c(0.65, 0.5)) +
    xlab("Genes with xTWAS hit(s)") +
    ylab("GWAS Traits") +
    labs(fill = "Modality of\ngene's top hit")
save_panel(p_c, "Fig4c_xtwas_genes_per_trait_by_top_modality", 5.8, 3.7)
full_join(hit_types20,
          top_hits |> dplyr::select(trait, gene_id, top_hit_modality = modality, top_hit_p_coloc = COLOC.PP4),
          by = c("trait", "gene_id")) |>
    write_tsv(file.path(O, "source_data", "Fig4c_xtwas_genes_per_trait.tsv"))

## ---- Fig 4: a and b stacked on the left, c on the right ------------------------------------------------------------
left <- plot_grid(p_a, p_b, ncol = 1, labels = c("a", "b"), label_size = 12)
fig4 <- plot_grid(left, p_c, nrow = 1, labels = c("", "c"), label_size = 12, rel_widths = c(1, 1.45))
ggsave(file.path(O, "Fig4.pdf"), fig4, width = 11.5, height = 4.8, device = cairo_pdf)
ggsave(file.path(O, "Fig4.png"), fig4, width = 11.5, height = 4.8, dpi = 400, bg = "white",
       device = ragg::agg_png)   # ragg: text measured and drawn with the same fonts (png() caused overlaps)
cat("[done]\n")


### END