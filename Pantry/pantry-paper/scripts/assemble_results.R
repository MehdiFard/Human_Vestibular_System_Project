rm(list = ls())

# Set dir
setwd("Pantry/")

library(tidyverse)
library(data.table)
library(readr)
library(dplyr)
library(tibble)
library(stringr)

modalities <- c("alt_polyA", "alt_TSS", "expression", "isoforms", "RNA_editing", "splicing", "stability")

# tissues <- read_lines("data/gtex/tissues.txt")
tissues <- c("saccule", "utricle")

traits <- read_tsv("pheast/input/gwas_sumstat/gwas_metadata.tsv",
            col_types = cols(trait = "c", .default = "-")) |>
            pull()

###########
## Genes ##
###########

genes <- rtracklayer::import("phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene",
            gene_biotype %in% c("protein_coding", "lncRNA")) |>
    pull(gene_id)

# dir.create("figures/source_data/", recursive = TRUE)
write_lines(genes, "figures/source_data/genes_pcg_lncrna.txt")

###################
## Geuvadis QTLs ##
###################

## Geuvadis top associations per gene per modality

# geuv_assoc <- tibble(modality = modalities) |>
#     reframe({
#         fname <- str_glue("pheast/output/saccule_w8/qtl/{modality}.cis_qtl.txt.gz")
#         if (modality %in% c("expression", "stability")) {
#             read_tsv(fname, col_types = "c-----cc-------cc-") |>
#                 mutate(gene_id = phenotype_id)
#         } else {
#             read_tsv(fname, col_types = "c-----cc-------cc-c-") |>
#                 mutate(gene_id = group_id) |>
#                 select(-group_id)
#         }
#     }, .by = modality) |>
#     filter(gene_id %in% genes) |>
#     relocate(gene_id, .before = 1) |>
#     arrange(gene_id, modality)

# write_tsv(geuv_assoc, "data/processed/geuvadis.sep.assoc.tsv")

geuv_assoc <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/output/saccule/qtl/{modality}.cis_qtl.txt.gz")
        df <- read_tsv(fname, col_types = cols(.default = col_character()))
        if ("group_id" %in% colnames(df)) {
            df <- df |>
                mutate(gene_id = group_id) |>
                select(-group_id)
        } else if ("phenotype_id" %in% colnames(df)) {
            df <- df |>
                mutate(gene_id = phenotype_id)
        } else {
            stop("Neither group_id nor phenotype_id found in ", fname)
        }
        df
    }, .by = modality) |>
    filter(gene_id %in% genes) |>
    relocate(gene_id, .before = 1) |>
    arrange(gene_id, modality) |>
    mutate(qval = as.numeric(qval)) # Ensure qval is numeric

# # To get the counts of significant associations per modality
# modality_counts <- geuv_assoc %>%
#     filter(modality %in% modalities, qval < 0.05) %>%
#     group_by(modality) %>%
#     summarise(n = n(), .groups = "drop") %>%
#     complete(modality = modalities, fill = list(n = 0)) 
# print(modality_counts)

# Write out the data
write_tsv(geuv_assoc, "figures/source_data/geuvadis.sep.assoc.tsv")



## Geuvadis QTLs (mapped separately per modality)

# geuv_qtls_sep <- tibble(modality = modalities) |>
#     reframe({
#         fname <- str_glue("pheast/output/saccule_w8/qtl/{modality}.cis_independent_qtl.txt.gz")
#         if (modality %in% c("expression", "stability")) {
#             read_tsv(fname, col_types = "c-----cc-------ci") |>
#                 mutate(gene_id = phenotype_id)
#         } else {
#             read_tsv(fname, col_types = "c-----cc-------cc-i") |>
#                 mutate(gene_id = group_id) |>
#                 select(-group_id)
#         }
#     }, .by = modality) |>
#     filter(gene_id %in% genes) |>
#     select(gene_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
#     arrange(gene_id, modality, rank)

# write_tsv(geuv_qtls_sep, "data/processed/geuvadis.sep.qtls.tsv")


geuv_qtls_sep <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/output/saccule/qtl/{modality}.cis_independent_qtl.txt.gz")
        if (!file.exists(fname)) {
            warning("Missing file: ", fname)
            return(tibble())}
        df <- read_tsv(fname, col_types = cols(.default = col_character(), rank = col_integer()))
        if ("group_id" %in% colnames(df)) {
            df <- df |> mutate(gene_id = group_id) |> select(-group_id)
        } else if ("phenotype_id" %in% colnames(df)) {
            df <- df |> mutate(gene_id = phenotype_id)
        } else {
            warning("No ID column found in: ", fname)
            return(tibble())}
        df
    }, .by = modality) |>
    filter(gene_id %in% genes) |>
    select(any_of(c("gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    arrange(gene_id, modality, rank)
View(geuv_qtls_sep)
dim(geuv_qtls_sep)

# Write out the data
write_tsv(geuv_qtls_sep, "figures/source_data/geuvadis.sep.qtls.tsv")


## Geuvadis QTLs (combined-modality mapping)

# geuv_qtls_comb <- read_tsv("data/geuvadis/qtl/all.cis_independent_qtl.txt.gz",
#                       col_types = "c-----cc-------cc-i") |>
#     separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
#     select(gene_id = group_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
#     filter(gene_id %in% genes) |>
#     arrange(gene_id, rank)
# write_tsv(geuv_qtls_comb, "data/processed/geuvadis.comb.qtls.tsv")


# Modified
geuv_qtls_comb <- read_tsv("pheast/output/saccule/qtl/cross_modality.cis_independent_qtl.txt.gz",
        col_types = cols(.default = col_character(), rank = col_integer())) |>
    separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    mutate(gene_id = if ("group_id" %in% names(pick(everything()))) group_id else phenotype_id) |>
    select(any_of(c("gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    filter(gene_id %in% genes) |>
    arrange(gene_id, rank)
dim(geuv_qtls_comb)
View(geuv_qtls_comb)

# Write out the data
write_tsv(geuv_qtls_comb, "figures/source_data/geuvadis.comb.qtls.tsv")


###############
## GTEx QTLs ##
###############

## GTEx top associations per gene per modality

# gtex_assoc <- crossing(tissue = tissues,
#                         modality = modalities) |>
#     reframe({
#         fname <- str_glue("data/gtex/qtl/{tissue}/{modality}.cis_qtl.txt.gz")
#         if (file.exists(fname)) {
#             if (modality %in% c("expression", "stability")) {
#                 read_tsv(fname, col_types = "c-----cc-------cc-") |>
#                     mutate(gene_id = phenotype_id)
#             } else {
#                 read_tsv(fname, col_types = "c-----cc-------cc-c-") |>
#                     mutate(gene_id = group_id) |>
#                     select(-group_id)
#             }
#         } else {
#             tibble()
#         }
#     }, .by = c(tissue, modality)) |>
#     filter(gene_id %in% genes) |>
#     relocate(gene_id, .before = 2) |>
#     arrange(tissue, gene_id, modality)

# write_tsv(gtex_assoc, "data/processed/gtex.sep.assoc.tsv.gz")


gtex_assoc <- crossing(tissue = tissues, modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/output/{tissue}/qtl/{modality}.cis_qtl.txt.gz")
        
        if (!file.exists(fname)) return(tibble())
        
        df <- read_tsv(fname, col_types = cols(.default = col_character()))
        
        if ("group_id" %in% names(df)) {
            df <- df |> mutate(gene_id = group_id) |> select(-group_id)
        } else if ("phenotype_id" %in% names(df)) {
            df <- df |> mutate(gene_id = phenotype_id)
        } else {
            warning("No gene ID column found in: ", fname)
            return(tibble())
        }
        
        df
    }, .by = c(tissue, modality)) |>
    filter(gene_id %in% genes) |>
    relocate(gene_id, .before = 2) |>
    arrange(tissue, gene_id, modality)
# Write out the data
write_tsv(gtex_assoc, "figures/source_data/gtex.sep.assoc.tsv.gz")


## GTEx QTLs (mapped separately per modality)

# gtex_qtls_sep <- crossing(tissue = tissues,
#                           modality = modalities) |>
#     reframe({
#         fname <- str_glue("data/gtex/qtl/{tissue}/{modality}.cis_independent_qtl.txt.gz")
#         if (file.exists(fname)) {
#             if (modality %in% c("expression", "stability")) {
#                 read_tsv(fname, col_types = "c-----cc-------ci") |>
#                     mutate(gene_id = phenotype_id)
#             } else {
#                 read_tsv(fname, col_types = "c-----cc-------cc-i") |>
#                     mutate(gene_id = group_id) |>
#                     select(-group_id)
#             }
#         } else {
#             tibble()
#         }
#     }, .by = c(tissue, modality)) |>
#     filter(gene_id %in% genes) |>
#     select(tissue, gene_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
#     arrange(tissue, gene_id, modality, rank)

# write_tsv(gtex_qtls_sep, "data/processed/gtex.sep.qtls.tsv.gz")


gtex_qtls_sep <- crossing(tissue = tissues, modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/output/{tissue}/qtl/{modality}.cis_independent_qtl.txt.gz")
        
        if (!file.exists(fname)) return(tibble())

        df <- read_tsv(fname, col_types = cols(.default = col_character(), rank = col_integer()))

        if ("group_id" %in% names(df)) {
            df <- df |> mutate(gene_id = group_id) |> select(-group_id)
        } else if ("phenotype_id" %in% names(df)) {
            df <- df |> mutate(gene_id = phenotype_id)
        } else {
            warning("Missing 'group_id' and 'phenotype_id' in: ", fname)
            return(tibble())
        }

        df
    }, .by = c(tissue, modality)) |>
    filter(gene_id %in% genes) |>
    select(any_of(c("tissue", "gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    arrange(tissue, gene_id, modality, rank)
# Write out the data
write_tsv(gtex_qtls_sep, "figures/source_data/gtex.sep.qtls.tsv.gz")


## GTEx QTLs (combined-modality mapping)

# gtex_qtls_comb <- tibble(tissue = tissues) |>
#     reframe({
#         fname <- str_glue("data/gtex/qtl/{tissue}/all.cis_independent_qtl.txt.gz")
#         if (file.exists(fname)) {
#             read_tsv(fname, col_types = "c-----cc-------cc-i")
#         } else {
#             tibble()
#         }
#     }, .by = tissue) |>
#     separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
#     select(tissue, gene_id = group_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
#     filter(gene_id %in% genes) |>
#     arrange(tissue, gene_id, rank)
# write_tsv(gtex_qtls_comb, "data/processed/gtex.comb.qtls.tsv.gz")


gtex_qtls_comb <- tibble(tissue = tissues) |>
    reframe({
        fname <- str_glue("pheast/output/{tissue}/qtl/cross_modality.cis_independent_qtl.txt.gz")
        if (!file.exists(fname)) return(tibble())
        df <- read_tsv(fname, col_types = cols(.default = col_character(), rank = col_integer()))
        # Assign gene_id based on available column
        if ("group_id" %in% names(df)) {
            df <- df |> mutate(gene_id = group_id)
        } else if ("phenotype_id" %in% names(df)) {
            df <- df |> mutate(gene_id = phenotype_id)
        } else {
            warning("Missing gene ID columns in: ", fname)
            return(tibble())
        }
        df
    }, .by = tissue) |>
    separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    select(any_of(c("tissue", "gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    filter(gene_id %in% genes) |>
    arrange(tissue, gene_id, rank)
# Write out the data
write_tsv(gtex_qtls_comb, "figures/source_data/gtex.comb.qtls.tsv.gz")



###################
## Geuvadis TWAS ##
###################

## Geuvadis heritability/r2/pvals

# geuv_hsq <- tibble(modality = modalities) |>
#     reframe(
#         read_tsv(str_glue("data/geuvadis/twas_models/{modality}.profile"),
#                     col_types = "cccccccc-cccc-c"),
#         .by = modality
#     ) |>
#     rename(phenotype_id = id) |>
#     mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
#             .before = 2) |>
#     filter(gene_id %in% genes)
# write_tsv(geuv_hsq, "data/processed/geuvadis.hsq.tsv.gz")


geuv_hsq <- tibble(modality = modalities) |>
    reframe(
        read_tsv(str_glue("pheast/intermediate/saccule/twas/{modality}.profile"),
                    col_types = "cccccccc-cccc-c"),
        .by = modality
    ) |>
    rename(phenotype_id = id) |>
    mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
            .before = 2) |>
    filter(gene_id %in% genes)
# Write out the data
write_tsv(geuv_hsq, "figures/source_data/geuvadis.hsq.tsv")



## Geuvadis TWAS

# geuv_twas <- crossing(modality = modalities,
#                         trait = traits) |>
#     reframe(
#         read_tsv(str_glue("data/geuvadis/twas/{modality}/fusion.Geuvadis.{modality}.{trait}.tsv"),
#                     col_types = "c---c--------c--ccccccc"),
#         .by = c(modality, trait)
#     ) |>
#     filter(!is.na(TWAS.P)) |>
#     mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
#     filter(gene_id %in% genes) |>
#     select(trait, gene_id, modality, phenotype_id = ID, HSQ:COLOC.PP4) |>
#     arrange(trait, gene_id, modality, phenotype_id)
# write_tsv(geuv_twas, "data/processed/geuvadis.twas.tsv.gz")


geuv_twas <- crossing(modality = modalities,
                        trait = traits) |>
    reframe(
        read_tsv(str_glue("pheast/twas_assoc/saccule/output/{modality}/{modality}.{trait}.tsv"),
                    col_types = "c---c--------c--ccccccc"),
        .by = c(modality, trait)
    ) |>
    filter(!is.na(TWAS.P)) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
    filter(gene_id %in% genes) |>
    select(trait, gene_id, modality, phenotype_id = ID, HSQ:COLOC.PP4) |>
    arrange(trait, gene_id, modality, phenotype_id)
# Write out the data
write_tsv(geuv_twas, "figures/source_data/geuvadis.twas.tsv.gz")


###############
## GTEx TWAS ##
###############

## GTEx heritability/r2/pvals

# gtex_hsq <- crossing(tissue = tissues,
#                         modality = modalities) |>
#     reframe({
#         fname <- str_glue("data/gtex/twas_models/{tissue}/{modality}.profile")
#         if (file.exists(fname)) {
#             read_tsv(fname, col_types = "cccccccc-cccc-c")
#         } else {
#             tibble()
#         }
#     }, .by = c(tissue, modality)) |>
#     rename(phenotype_id = id) |>
#     mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
#             .before = 3) |>
#     filter(gene_id %in% genes)
# write_tsv(gtex_hsq, "data/processed/gtex.hsq.tsv.gz")


gtex_hsq <- crossing(tissue = tissues,
                        modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/intermediate/{tissue}/twas/{modality}.profile")
        if (file.exists(fname)) {
            read_tsv(fname, col_types = "cccccccc-cccc-c")
        } else {
            tibble()
        }
    }, .by = c(tissue, modality)) |>
    rename(phenotype_id = id) |>
    mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
            .before = 3) |>
    filter(gene_id %in% genes)
# Write out the data
write_tsv(gtex_hsq, "figures/source_data/gtex.hsq.tsv.gz")



## GTEx TWAS (significant hits only)

# gtex_twas <- tibble(tissue = tissues) |>
#     reframe({
#         fname <- str_glue("data/gtex/twas/{tissue}.twas_hits.tsv")
#         if (file.exists(fname)) {
#             read_tsv(fname, col_types = "ccccccccccccccccccccccccc")
#         } else {
#             tibble()
#         }
#     }, .by = tissue) |>
#     mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
#     filter(gene_id %in% genes) |>
#     select(tissue, trait = TRAIT, gene_id, modality = MODALITY, phenotype_id = ID, HSQ, MODEL, TWAS.Z, TWAS.P, COLOC.PP0, COLOC.PP1, COLOC.PP2, COLOC.PP3, COLOC.PP4) |>
#     arrange(tissue, trait, gene_id, modality, phenotype_id)
# write_tsv(gtex_twas, "data/processed/gtex.twas_hits.tsv.gz")


gtex_twas <- tibble(tissue = tissues) |>
    reframe({
        fname <- str_glue("pheast/twas_assoc/{tissue}/output/alt_polyA/twas_hits.tsv") # For a given tissue, all modalites have the same "twas_hits.tsv" file, so I just pick "alt_polyA".
        if (file.exists(fname)) {
            read_tsv(fname, col_types = "ccccccccccccccccccccccccc")
        } else {
            tibble()
        }
    }, .by = tissue) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
    filter(gene_id %in% genes) |>
    select(tissue, trait = TRAIT, gene_id, modality = MODALITY, phenotype_id = ID, HSQ, MODEL, TWAS.Z, TWAS.P, COLOC.PP0, COLOC.PP1, COLOC.PP2, COLOC.PP3, COLOC.PP4) |>
    arrange(tissue, trait, gene_id, modality, phenotype_id)
# Write out the data
write_tsv(gtex_twas, "figures/source_data/gtex.twas_hits.tsv.gz")

### END