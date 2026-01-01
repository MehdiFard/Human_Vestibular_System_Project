library(tidyverse)
library(readr)
library(dplyr)

get_phenotypes <- function(modality) {
    df <- read_tsv(str_glue("Pantry/pheast_no_outlier/input/vestibular/phenotypes/{modality}.bed.gz"),
                   col_types = cols(phenotype_id = "c", .default = "-"))
    if (modality %in% c("expression", "stability")) {
        df <- mutate(df, gene_id = phenotype_id)
    } else {
        group_file <- str_glue("Pantry/pheast_no_outlier/input/vestibular/phenotypes/{modality}.phenotype_groups.txt")
        groups <- read_tsv(group_file, col_types = "cc",
                           col_names = c("phenotype_id", "gene_id"))
        df <- left_join(df, groups, by = "phenotype_id")
    }
    df
}

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability",
    RNA_editing = "RNA editing"
)


genes <- read_lines("figures/data_tables/genes_pcg_lncrna.txt")
head(genes)

#############
## Samples ##
#############

# "We generated data on these six modalities of transcriptome regulation for the
# 445 samples in Geuvadis7 dataset and all 17,350 samples across 54 tissues in
# the GTEx Project V8 release."

read_lines("data/geuvadis/samples.txt") |> length()

read_tsv("data/gtex/samples_tissues.txt", col_names = c("sample_id", "tissue"),
         col_types = "cc") |>
    nrow()

################
## Phenotypes ##
################

phenos <- tibble(modality = names(modalities)) |>
    reframe(get_phenotypes(modality), .by = modality) |>
    mutate(modality = factor(modalities[modality], level = modalities)) |>
    filter(gene_id %in% genes)
head(phenos)
length(unique(phenos$gene_id)) # No. unique genes with phenotypes
table(phenos$modality)

## Table 1

# No. phenotypes produced:
count(phenos, modality)
nrow(phenos)

# No. genes with phenotypes, mean/SD phenotypes per gene
phenos |>
    count(modality, gene_id) |>
    summarise(
        per_gene_mean = mean(n),
        per_gene_sd = sd(n),
        n_genes = n(),
        .by = modality
    )
phenos |>
    count(gene_id) |>
    summarise(
        per_gene_mean = mean(n),
        per_gene_sd = sd(n),
        n_genes = n()
    )

data_t1 <- phenos |>
    count(modality, gene_id, name = "n_phenos")

write_tsv(data_t1, "figures/data_tables/Table_1.txt")

##########
## QTLs ##
##########

qtls_sep <- read_tsv("figures/data_tables/vestibular.sep.qtls.tsv", col_types = "cicccd") |>
    dplyr::mutate(modality = factor(modalities[modality], level = modalities))

qtls_comb <- read_tsv("figures/data_tables/vestibular.comb.qtls.tsv", col_types = "cicccd") |>
    dplyr::mutate(modality = factor(modalities[modality], level = modalities))

## We start with reporting combined modality QTLs
head(qtls_comb)
View(qtls_comb)
nrow(qtls_comb)
n_distinct(qtls_comb$gene_id)

qtls_comb |>
    filter(modality == "Expression") |>
    distinct(gene_id) |>
    count()
table(qtls_comb$modality)
which(qtls_comb$gene_id == "ENSG00000187323") # DCC
which(qtls_comb$gene_id == "ENSG00000167232") # ZNF91

names(qtls_comb)
a = qtls_comb |>
    count(modality, variant_id) |>
    summarise(
        per_variant_mean = mean(n),
        per_variant_sd = sd(n),
        n_variants = n(),
        .by = modality
    )
a

sum(a$n_variants)
length(unique(qtls_comb$variant_id))

qtls_comb |>
    distinct(modality, variant_id) |>
    count(modality, sort = TRUE)

qtls_comb |>
    group_by(modality, variant_id) |>
    count(modality) #|>

### Venn diagram of variants between modalities
library(dplyr)
library(ggVennDiagram)  # install.packages("ggVennDiagram") if needed

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c",
    `RNA editing` = "#984210")


# Build a named list of unique variant IDs per modality
modality_sets <- qtls_comb %>%
  group_by(modality) %>%
  summarise(variants = list(unique(variant_id)), .groups = "drop") %>%
  tibble::deframe()
head(modality_sets)


# Match colors to the modalities present in the data
set_cols <- modality_colors[names(modality_sets)]

ggVennDiagram(
    modality_sets,
    label = "count",
    set_color = set_cols,
    set_fill = scales::alpha(set_cols, 0.45), # translucent fills to see overlaps
    label_alpha = 0
    ) +
    ggtitle("Variant overlap by modality")



# Or subset to 4 modalities for clarity
subset_sets <- modality_sets[c(
    "Expression",
    "RNA stability",
    "Intron excision",
    "Isoform ratio"
)]

ggVennDiagram(subset_sets, label_alpha = 0,
                set_color = modality_colors[names(subset_sets)],
                set_fill = scales::alpha(modality_colors[names(subset_sets)], 0.45)) +
                # scale_fill_viridis_c(option = "plasma") +
                ggtitle("Variant overlap by modality (subset)")


### It's not what I wanted, but ggvenn would do the job, however it only works for 2-4 sets
library(ggvenn)

ggvenn(
    subset_sets, 
    #   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    fill_color = modality_colors[names(subset_sets)],
    stroke_size = 0.5, set_name_size = 4
    ) +
    # theme_void() +
    ggtitle("Sig. Variant overlap by modality")
    

# Save figure
ggsave("figures/venn_variant_overlap_by_modality.png",
        width = 8,
        height = 8,
        dpi = 300,
        bg = "white")



### UpSet plot of variants between modalities
library(tidyr)
library(ComplexUpset)

variant_df <- qtls_comb %>%
    distinct(variant_id, modality) %>%
    mutate(present = TRUE) %>%
    tidyr::pivot_wider(
        names_from = modality,
        values_from = present,
        values_fill = FALSE
    )

# UpStet plot
ComplexUpset::upset(
    variant_df,
    intersect = names(modality_sets),
    base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size()
    ),
    queries=list(
    upset_query(set='Expression', color='#e41a1c', fill='#e41a1c'),
    upset_query(set='Isoform ratio', color='#377eb8', fill='#377eb8'),
    upset_query(set='Intron excision', color='#4daf4a', fill='#4daf4a'),
    upset_query(set='Alt. TSS', color='#984ea3', fill='#984ea3'),
    upset_query(set='Alt. polyA', color='#ff7f00', fill='#ff7f00'),
    upset_query(set='RNA stability', color='#fdc11c', fill='#fdc11c'),
    upset_query(set='RNA editing', color='#984210', fill='#984210')
    )
) +
    labs(title = "xQTLs overlap across regulatory modalities")

# Save figure
ggsave("figures/upset_xQTLs_overlap_by_modality.png",
        width = 12,
        height = 8,
        dpi = 300,
        # bg = "white"
        )


### Venn diagram of overlapped genes across modalities
# Build a named list of unique xGene IDs per modality
modality_sets_genes <- qtls_comb %>%
    group_by(modality) %>%
    summarise(variants = list(unique(gene_id)), .groups = "drop") %>%
    tibble::deframe()
head(modality_sets_genes)

# Match colors to the modalities present in the data
set_cols <- modality_colors[names(modality_sets_genes)]


subset_sets_genes <- modality_sets_genes[c(
    "Expression",
    "RNA stability",
    "Intron excision",
    "Isoform ratio"
)]

ggvenn(
    subset_sets_genes, 
    #   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    fill_color = modality_colors[names(subset_sets_genes)],
    stroke_size = 0.5, set_name_size = 4
    ) +
    # theme_void() +
    ggtitle("Sig. xGenes overlap by modality")
    

# Save figure
ggsave("figures/venn_xgenes_overlap_by_modality.png",
        width = 8,
        height = 8,
        dpi = 300,
        bg = "white")


### UpSet plot of genes between modalities
library(tidyr)
library(ComplexUpset)

gene_df <- qtls_comb %>%
    distinct(gene_id, modality) %>%
    mutate(present = TRUE) %>%
    tidyr::pivot_wider(
        names_from = modality,
        values_from = present,
        values_fill = FALSE
    )

# UpStet plot
ComplexUpset::upset(
    gene_df,
    intersect = names(modality_sets),
    base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size()
    ),
    queries=list(
    upset_query(set='Expression', color='#e41a1c', fill='#e41a1c'),
    upset_query(set='Isoform ratio', color='#377eb8', fill='#377eb8'),
    upset_query(set='Intron excision', color='#4daf4a', fill='#4daf4a'),
    upset_query(set='Alt. TSS', color='#984ea3', fill='#984ea3'),
    upset_query(set='Alt. polyA', color='#ff7f00', fill='#ff7f00'),
    upset_query(set='RNA stability', color='#fdc11c', fill='#fdc11c'),
    upset_query(set='RNA editing', color='#984210', fill='#984210')
    )
) +
    labs(title = "xGenes overlap across regulatory modalities")
    # ggtitle("Sig. xGenes overlap across regulatory modalities")

# Save figure
ggsave("figures/upset_xGenes_overlap_by_modality.png",
        width = 16,
        height = 8,
        dpi = 300,
        # bg = "white"
        )



### Seperated mapping QTLs
nrow(qtls_sep)
n_distinct(qtls_sep$gene_id)
qtls_sep |>
    filter(modality == "Expression") |>
    distinct(gene_id) |>
    count()
table(qtls_sep$modality)

qtls_sep |>
    distinct(modality, gene_id) |>
    count(modality, sort = TRUE)

which(qtls_sep$gene_id == "ENSG00000187323")
which(qtls_sep$gene_id == "ENSG00000167232")

sig_iso <- qtls_sep[which(qtls_sep$gene_id == "ENSG00000167232"), ]
View(sig_iso)

# # A tibble: 1 × 6
#   gene_id            rank     modality        phenotype_id                            variant_id          pval_beta
# 1 ENSG00000167232     1       Isoform_ratio   ENSG00000167232__ENST00000300619     19_23314211_G_C         0.000346
# Its present in fine mapping 

# "cis-eQTLs found in 7215 genes, more than 3.2 times the second-most abundant xQTL group isoform ratio."
qtls_comb |>
    distinct(modality, gene_id) |>
    count(modality, sort = TRUE)
qtls_comb |>
    mutate(modality = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(gene_id, modality) |>
    summarise(modalities = str_c(sort(modality), collapse = "_"),
                .by = gene_id) |>
    count(modalities)

qtls_comb |>
    count(gene_id) |>
    count(n, name = "nn") |>
    with(sum(nn[n > 1]))

## Abstract
# "We applied Pantry to Geuvadis and GTEx data, and found that 4,768 of the genes
# with no identified expression QTL in Geuvadis had QTLs in at least one other
# transcriptional modality, resulting in a 66% increase in genes over expression
# QTL mapping."
# Also later:
# "The 66% increase in the number of xGenes highlights..."
qtls_comb |>
    summarise(has_expr = "Expression" %in% modality,
              .by = gene_id) |>
    count(has_expr)
# Mention the denominator?
egenes <- unique(qtls_comb$gene_id[qtls_comb$modality == "Expression"])
length(genes[!(genes %in% egenes)])

## "We discovered comparable numbers of xQTLs as for Geuvadis, which varied across
# tissues due to factors such as sample size, but generally found non-expression
# xQTLs in thousands of genes per tissue for which no eQTLs were found in our data,
# resulting in a 61% increase of xGenes over eGenes alone on average

# See figureS2.R

## "Cross-modality mapping resulted in fewer total xQTLs per gene, from 2.94 on average to 1.76"

# See figure3ab.R

## Cover letter?
# "Applying pantry to data from Geuvadis and GTEx data, we identify 736,810
# conditionally independent QTL associations...
# Analyzing expression alone would result in 441,280 eQTLs, 91,841 of which
# are better described by another transcriptional phenotype."
qtls_sep_gtex <- read_tsv("data/processed/gtex.sep.qtls.tsv.gz", col_types = "ccicccid") |>
    mutate(modality = factor(modalities[modality], levels = modalities))
qtls_comb_gtex <- read_tsv("data/processed/gtex.comb.qtls.tsv.gz", col_types = "ccicccid") |>
    mutate(modality = factor(modalities[modality], levels = modalities))
nrow(qtls_comb) + nrow(qtls_comb_gtex)
sum(qtls_sep_gtex$modality == "Expression") + sum(qtls_sep$modality == "Expression")

qtls_all_comb <- bind_rows(
    qtls_comb |> mutate(tissue = "Geuvadis"),
    qtls_comb_gtex
)
qtls_all_sep <- bind_rows(
    qtls_sep |> mutate(tissue = "Geuvadis"),
    qtls_sep_gtex
)
qtls_venn_sep_expr <- bind_rows(
    qtls_all_comb |> mutate(method = "combined"),
    qtls_all_sep |> mutate(method = "separate")
) |>
    group_by(tissue, gene_id) |>
    filter(length(unique(method)) == 2) |> # Include only genes found with both methods to estimate "consolidation" of redundant QTLs
    ungroup() |>
    group_by(tissue, modality, gene_id) |>
    mutate(methods = str_c(sort(unique(method)), collapse = "_")) |>
    ungroup() |>
    filter(methods == "separate",
           modality == "Expression")
nrow(qtls_venn_sep_expr)

##########
## TWAS ##
##########

twas <- read_tsv("data/processed/geuvadis.twas.tsv.gz", col_types = "ccccdcddddddd") |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(TWAS.P < 5e-8 / 6)

# "We found 10,065 significant hits..."
nrow(twas)

# "...across 80 traits..."
n_distinct(twas$trait)

# "...involving 4,304 unique RNA phenotypes..."
twas |>
    distinct(gene_id, phenotype_id) |>
    nrow()

# "...for 1,934 genes."
n_distinct(twas$gene_id)

# "Of the 4,487 unique trait-gene pairs among these hits..."
twas |>
    distinct(trait, gene_id) |>
    nrow()

# "...51.3% involved only non-expression RNA phenotypes"
twas |>
    summarise(has_expr = "Expression" %in% modality,
              .by = c(trait, gene_id)) |>
    with(mean(!has_expr))

# Abstract: "...and enhances identification of regulatory mechanisms underlying GWAS signal in a large fraction of previously associated gene-trait pairs"
twas |>
    filter("Expression" %in% modality,
           .by = c(trait, gene_id)) |>
    summarise(has_other = any(modality != "Expression"),
              .by = c(trait, gene_id)) |>
    with(mean(has_other))

# "our in-depth analysis of the transcriptome doubles the number of
# gene-to-phenotype discoveries, ???,"
twas_gtex <- read_tsv("data/processed/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    mutate(modality = factor(modalities[modality], levels = modalities))

twas_all <- bind_rows(select(mutate(twas, tissue = "GEUVADIS"), trait, tissue, gene_id, modality, COLOC.PP4),
                      select(twas_gtex, trait, tissue, gene_id, modality, COLOC.PP4))
twas_all |>
    filter(COLOC.PP4 >= 0.8) |>
    filter(modality == "Expression") |>
    distinct(trait, tissue, gene_id) |>
    nrow()
twas_all |>
    filter(COLOC.PP4 >= 0.8) |>
    distinct(trait, tissue, gene_id) |>
    nrow()

## GTEx
# "identifying colocalizing hits for thousands more trait-tissue-gene triplets than would be found using expression alone"
twas_gtex |>
    filter(COLOC.PP4 >= 0.8) |>
    summarise(has_expr = "Expression" %in% modality,
              .by = c(trait, tissue, gene_id)) |>
    summarise(n_non_expr = sum(!has_expr),
              n_expr = sum(has_expr)) |>
    mutate(factor_increase = (n_non_expr + n_expr) / n_expr)

## Discussion: "Notably, for more than two-fifths of the gene-trait pairs with
# previous TWAS hits from gene expression analysis, we now identify at least one
# other regulation modality connecting them."
# counting trait-gene pairs:
twas |>
    filter("Expression" %in% modality,
           .by = c(trait, gene_id)) |>
    summarise(has_non_expr = any(modality != "Expression"),
              .by = c(trait, gene_id)) |>
    with(mean(has_non_expr))
