# Figure 4: TWAS
rm(list = ls())

library(tidyverse)
library(patchwork)
library(readr)
library(ggplot2)
library(tibble)
library(forcats)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability",
    RNA_editing = "RNA editing"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision ratio` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c",
    `RNA editing` = "#984210"
)

twas <- read_tsv("figures/data_tables/vestibular.twas.tsv.gz", col_types = "ccccdcddddddd") |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(TWAS.P < 5e-8 / 7)
View(twas)
dim(twas)
summary(twas$TWAS.P)
length(unique(twas$phenotype_id))
length(unique(twas$gene_id))

# ssize <- read_tsv("pheast/input/gwas_sumstat/gwas_metadata.tsv",
#                     col_types = cols(trait = "c", sample_size = "i", .default = "-")) |>
#     mutate(
#         # Sample_Size = format(Sample_Size, big.mark=","),
#         # Sample_Size = scales::label_number(scale_cut = scales::cut_short_scale())(Sample_Size), # Currently fails due to bug in scales package
#         sample_size = sample_size |> round(-3) |> str_sub(1, -4) |> str_c("K"),
#     ) |>
#     deframe()

ssize <- read_tsv("Pantry/pheast_no_outlier/input/gwas_sumstat/gwas_metadata.tsv",
                    col_types = cols(trait = "c", sample_size = "i", .default = "-")) |>
    mutate(
        # Sample_Size = format(Sample_Size, big.mark=","),
        # Sample_Size = scales::label_number(scale_cut = scales::cut_short_scale())(Sample_Size), # Currently fails due to bug in scales package
        sample_size = sample_size |> round(-3) |> str_sub(1, -4) |> str_c("K"),
    ) |>
    deframe()

View(ssize)

##################
## Panel A left ## Genes with expression vs other modality TWAS hits
##################

d <- twas |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(trait, gene_id, modality_type) |>
    summarise(
        modality_hits = str_glue("{'Expression' %in% modality_type}_{'Other' %in% modality_type}"),
        .by = c(trait, gene_id)
    ) |>
    mutate(
        modality_hits = c(`TRUE_FALSE` = "Expression only",
                            `TRUE_TRUE` = "Expr. & other(s)",
                            `FALSE_TRUE` = "Non-expr. only")[modality_hits] |>
            fct_relevel("Expression only", "Expr. & other(s)", "Non-expr. only")
    )
View(d)

d20 <- d |>
    filter(trait %in% levels(fct_infreq(trait))[1:20])
View(d20)
table(d20$modality_hits)

trait_labels <- c(
        anx_eur	= "Anxiety",
        balance_meta_eur	= "Balance meta eur",
        balance_meta1 = "Balance meta",
        cognitive_skills_ext = "Cognitive skills",
        fluid_intelligence_eur = "Fluid intelligence",
        tinnitus_eur = "Tinnitus",
        frailty_gwas = "Frailty",
        alzheimer = "Alzheimer",
        parkinson_gp2_eur = "Parkinson ",
        hearing_hleur = "Hearing",
        icp_fa = "ICP FA",
        icp_md = "ICP MD",
        mcp_fa = "MCP FA",
        mcp_md = "MCP MD",
        non_cognitive_skills_ext = "Non-cognitive skills",
        depression_mdd = "Depression",
        scp_fa = "SCP FA",
        scp_md = "SCP MD",
        cerebellar_volume = "Cerebellar volume",
        cerebral_volume = "Cerebral volume",
        subcortical_volume = "Subcortical volume",
        ukbb_falls = "Falls-UKBB"
)

trait_labels <- str_glue("{trait_labels} (n={ssize[names(trait_labels)]})") |> purrr::set_names(names(trait_labels))
stopifnot(all(d20$trait %in% names(trait_labels)))

d20 |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality_hits)) +
    geom_bar(width = 0.8, color = "black") +
    scale_fill_manual(values = c("white", "gray", "#444444")) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.justification = c(0, 1),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.35),
    ) +
    xlab("Genes with xTWAS hit(s)") +
    ylab("Traits (top 20 by # of genes with hits)") +
    labs(fill = "Gene's xTWAS\nhits include")


# ggsave("figures/figure4/figure4a1.png", width = 5.8, height = 3.7, device = png)
ggsave("figures/paper/figure4a1_vestibular.pdf", width = 8.7, height = 5.55)

###################
## Panel A right ## Modality of top hit
###################

top_hits <- twas |>
    group_by(trait, gene_id) |>
    slice_min(TWAS.P, n = 1, with_ties = TRUE) |>
    slice_sample(n = 1) |>
    ungroup() |>
    filter(trait %in% levels(fct_infreq(trait))[1:20])
stopifnot(all(top_hits$trait %in% names(trait_labels)))

top_hits |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality)) +
    geom_bar(width = 0.8) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.justification = c(0, 1),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.5),
    ) +
    xlab("Genes with xTWAS hit(s)") +
    ylab("Traits (top 20 by # of genes with hits)") +
    labs(fill = "Modality of\ngene's top hit")

# ggsave("figures/figure4/figure4a2.png", width = 5.8, height = 3.7, device = png)
ggsave("figures/paper/figure4a2_vestibular.pdf", width = 8.7, height = 5.55)

data_4ab <- full_join(
    d20,
    top_hits |> dplyr::select(trait, gene_id, top_hit_modality = modality, top_hit_p_coloc = COLOC.PP4),
    by = c("trait", "gene_id")
)

write_tsv(data_4ab, "figures/data_tables/Figure_4ab_vestibular.txt")

#############
## Panel B ## Colocalization proportions
#############

top_hits |>
    summarise(prop_coloc = mean(COLOC.PP4 > 0.8),
              .by = c(modality)) |>
    ggplot(aes(x = prop_coloc, y = modality, fill = modality)) +
    geom_col(data = distinct(top_hits, modality) |> mutate(prop_coloc = 1),
             alpha = 0.5, width = 0.8, show.legend = FALSE) +
    geom_col(position = "dodge", width = 0.8, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
    ) +
    xlab("Prop. coloc.") +
    ylab(NULL)

# ggsave("figures/figure4/figure4b1.png", width = 2.5, height = 1.5, device = png)
ggsave("figures/paper/figure4b1_vestibular.pdf", width = 2.5, height = 1.5)

top_hits |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality, alpha = COLOC.PP4 > 0.8)) +
    facet_wrap(~ modality, ncol = 1) +
    geom_bar(width = 1, show.legend = FALSE) +
    # annotate("text", x = 0, y = 10, vjust = -0.5, label = "Traits", angle = 90) +
    geom_text(data = tibble(modality = factor(modalities[1], levels = modalities)),
                mapping = aes(alpha = NULL, y = NULL),
                x = 0, y = 10, vjust = -0.5, label = "Traits", angle = 90, size = 3.5,
                show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5)) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
    ) +
    xlab("Genes with xTWAS hit(s) per trait-modality") +
    ylab(NULL)

# ggsave("figures/figure4/figure4b2.png", width = 4.5, height = 3.5, device = png)
ggsave("figures/paper/figure4b2_vestibular.pdf", width = 4.5, height = 3.5)

#############
## Panel C ## GTEx colocalizations expression vs. other
#############

twas_gtex <- read_tsv("figures/source_data/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    filter(!is.na(COLOC.PP0)) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

expr_other_gtex <- twas_gtex |>
    filter(COLOC.PP4 > 0.8) |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(tissue, trait, gene_id, modality_type) |>
    summarise(
        modality_hits = str_c(sort(modality_type), collapse = "_"),
        .by = c(tissue, trait, gene_id)
    ) |>
    mutate(
        modality_hits = c(`Expression` = "Expression only",
                          `Expression_Other` = "Expr. & other(s)",
                          `Other` = "Non-expr. only")[modality_hits] |>
            fct_relevel("Expression only", "Expr. & other(s)", "Non-expr. only")
    )

data_4c <- expr_other_gtex |>
    mutate(tissue = fct_infreq(tissue)) |>
    summarise(n_pairs = n(),
              thousands = n() / 1000,
              .by = c(tissue, modality_hits))

ggplot(data_4c, aes(x = thousands, y = tissue, fill = modality_hits)) +
    geom_col(width = 1, color = "black") +
    scale_fill_manual(values = c("white", "gray", "#444444")) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.75),
        legend.text = element_text(size = 8),
    ) +
    scale_x_continuous(expand = c(0, 0.1)) +
    labs(fill = "Colocalizations\ninclude") +
    xlab("Colocalized gene-trait pairs (×1000)") +
    ylab("GTEx tissues")

# ggsave("figures/figure4/figure4c.png", width = 3.5, height = 3.5, device = png)
ggsave("figures/figure4/figure4c.pdf", width = 3.5, height = 3.5)

write_tsv(data_4c, "figures/source_data/Figure_4c.txt")


### END


gwas = data.table::fread("data/gwas_raw_sumstats/vertigo_GRCH38_sorted.tsv.gz")
head(gwas)
dim(gwas)
