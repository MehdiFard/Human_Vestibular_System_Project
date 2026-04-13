# Figure 2: QTL overview
rm(list = ls())
# install.packages("tidyverse")

library(tidyverse)
library(dplyr)
library(patchwork)
library(readr)
library(forcats)
library(ggplot2)
library(stringr)


modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability",
    RNA_editing = "RNA editing"
    
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c",
    `RNA editing` = "#984210")


qtls <- read_tsv("figures/data_tables/vestibular.comb.qtls.tsv", col_types = "cicccd") |>
    mutate(modality = factor(modalities[modality], level = modalities),
            rank = fct_other(as.character(rank), keep = as.character(1:9),
                            other_level = "10+"))
dim(qtls)
View(qtls)

#############
## Panel A ## xQTLs per rank per modality using separate mapping
#############

data_2a <- qtls |>
    dplyr::count(modality, gene_id, name = "n_QTLs") |>
    mutate(n_QTLs = as_factor(n_QTLs) |>
            fct_other(keep = as.character(1:4), other_level = "5+"),
            modality = fct_infreq(modality)) |>
    group_by(modality, n_QTLs) |>
    # summarise(thousands = n() / 1000,
    summarise(thousands = n() / 10,
            .groups = "drop")

dim(data_2a)
View(data_2a)

ggplot(data_2a, aes(x = thousands, y = modality, fill = n_QTLs)) +
    geom_col(width = 0.8) +
    # scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_fill_viridis_d() +
    labs(fill = "xQTLs") +
    # xlab("xGenes (×1000)") +
    xlab("xGenes (×10)") +
    ylab("Modality") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.65),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
    )

# Version 2: Just color by modality and ignoring xQTL counts
ggplot(data_2a, aes(x = thousands, y = modality, fill = modality)) +
    geom_col(width = 0.8) +
    # scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_fill_manual(values = modality_colors) +
    labs(fill = "Modality") +
    # xlab("xGenes (×1000)") +
    xlab("xGenes (×10)") +
    ylab("Modality") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.position = "none"
    )


# ggsave("figures/paper/figure2a.png", width = 3.2, height = 1.7, device = png)
# dir.create("figures/paper")
ggsave("figures/paper/figure2a_vestibular.pdf", width = 4.16, height = 2.08)

# ggsave("figures/paper/figure2a_vestibular.pdf",
#         units = "mm",
#         width = 70, height = 35)


# Write out data 
write_tsv(data_2a, "figures/data_tables/Figure_2a_vestibular.txt")
View(qtls)

#############
## Panel B ## QTLs per gene
#############

qtls_per_gene <- qtls |>
    dplyr::count(gene_id)  |>
    mutate(n_qtls = pmin(n, 5) |>
    as.character() |>
    str_replace("5", "5+"))

data_2b <- qtls_per_gene |>
    dplyr::count(n_qtls) |>
    # mutate(thousands = n / 1000)
    mutate(thousands = n / 10)

ggplot(data_2b, aes(x = n_qtls, y = thousands)) +
    geom_col(position = "dodge", width = 0.5, color = "black", fill = NA, linewidth = 0.5) +
    xlab("Total xQTLs per gene") +
    # ylab("Genes (×1000)") +
    ylab("Genes (×10)") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
    )

# ggsave("figures/paper/figure2b.png", width = 2.2, height = 2.1, device = png)
ggsave("figures/paper/figure2b_vestibular.pdf", width = 2.2, height = 2.1)
write_tsv(data_2b, "figures/data_tables/Figure_2b_vestibular.txt")


#############
## Panel C ## Proportion of xQTLs per rank associated with each modality
#############

rank_counts <- qtls |>
    dplyr::count(rank)

rank_counts <- rank_counts[1:5, ]


qtls |> 
    filter(rank %in% 1:5) |>
    ggplot(aes(x = rank, fill = modality)) +
    geom_bar(position = "fill", width = 0.8) +
    geom_text(aes(x = rank, y = 1.27, label = n, fill = NULL), data = rank_counts,
                hjust = 1, angle = 90, size = 3) +
    annotate("text", x = -0.2, y = 1.2, label = "n =", size = 3) +
    scale_fill_manual(values = modality_colors) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    # coord_cartesian(xlim = c(1, 10), ylim = c(0, 1), clip = "off") + # ylim prevents extension of axis line
    coord_cartesian(xlim = c(1, 5), ylim = c(0, 1), clip = "off") + 
    labs(fill = "Modality") +
    xlab("Rank of xQTL within gene") +
    ylab("Proportion of xQTLs") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = NULL,
        plot.margin = margin(30, 5.5, 5.5, 5.5),
    )

ggsave("figures/paper/figure2c_vestibular.png", width = 3.5, height = 2.3, device = png)
ggsave("figures/paper/figure2c_vestibular.pdf", width = 3.5, height = 2.3)

data_2c <- qtls |>
    dplyr::count(rank, modality) |>
    mutate(proportion = n / sum(n),
            .by = rank)

write_tsv(data_2c, "figures/data_tables/Figure_2c_vestibular.txt")

### END