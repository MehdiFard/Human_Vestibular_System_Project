rm(list = ls())

# Set dir
setwd("Pantry/")

library(tidyverse)
library(data.table)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)

modalities <- c("alt_polyA", "alt_TSS", "expression", "isoforms", "RNA_editing", "splicing", "stability")

# tissues <- read_lines("data/gtex/tissues.txt")
# tissues <- c("saccule", "utricle", "vestibular")
tissues <- c("vestibular")

traits <- read_tsv("Pantry/pheast_no_outlier/input/gwas_sumstat/gwas_metadata.tsv",
            col_types = cols(trait = "c", .default = "-")) |>
            pull()


###########
## Genes ##
###########

genes <- rtracklayer::import("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene",
            gene_biotype %in% c("protein_coding", "lncRNA")) |>
    pull(gene_id)

# dir.create("figures/data_tables/", recursive = T, showWarnings = F)
write_lines(genes, "figures/data_tables/genes_pcg_lncrna.txt")

###################
## Geuvadis QTLs ##
###################

## Geuvadis top associations per gene per modality
vestib_assoc <- tibble(modality = modalities) |>
    reframe({
        # fname <- str_glue("pheast/output/saccule/qtl/{modality}.cis_qtl.txt.gz")
        fname <- str_glue("Pantry/pheast_no_outlier/output/vestibular/qtl/{modality}.cis_qtl.txt.gz")
        df <- read_tsv(fname, col_types = cols(.default = col_character()))
        if ("group_id" %in% colnames(df)) {
            df <- df |>
                mutate(gene_id = group_id) |>
                dplyr::select(-group_id)
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
dim(vestib_assoc)
View(vestib_assoc)

# To get the counts of significant variant–phenotype pairs per modality
egene_counts <- vestib_assoc %>%
    dplyr::filter(qval < 0.05) %>%
    dplyr::distinct(modality, gene_id) %>%
    dplyr::count(modality) %>%
    tidyr::complete(modality = modalities, fill = list(n = 0))

print(egene_counts)
sum(egene_counts$n)

# eGenes unique
total_egenes <- vestib_assoc %>%
    filter(qval < 0.05) %>%
    distinct(gene_id) %>%
    nrow()
print(total_egenes)

# Write out the data
write_tsv(vestib_assoc, "figures/data_tables/vestibular.sep.assoc.tsv")


## Geuvadis QTLs (mapped separately per modality)
vestib_qtls_sep <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("Pantry/pheast_no_outlier/output/vestibular/qtl/{modality}.cis_independent_qtl.txt.gz")
        if (!file.exists(fname)) {
            warning("Missing file: ", fname)
            return(tibble())}
        df <- read_tsv(fname, col_types = cols(.default = col_character(), rank = col_integer()))
        if ("group_id" %in% colnames(df)) {
            df <- df |> mutate(gene_id = group_id) |> dplyr::select(-group_id)
        } else if ("phenotype_id" %in% colnames(df)) {
            df <- df |> mutate(gene_id = phenotype_id)
        } else {
            warning("No ID column found in: ", fname)
            return(tibble())}
        df
    }, .by = modality) |>
    dplyr::filter(gene_id %in% genes) |>
    dplyr::select(any_of(c("gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    dplyr::arrange(gene_id, modality, rank)
View(vestib_qtls_sep)
dim(vestib_qtls_sep)

# Write out the data
write_tsv(vestib_qtls_sep, "figures/data_tables/vestibular.sep.qtls.tsv")


# Modified
vestib_qtls_comb <- read_tsv("Pantry/pheast_no_outlier/output/vestibular/qtl/cross_modality.cis_independent_qtl.txt.gz",
        col_types = cols(.default = col_character(), rank = col_integer())) |>
    tidyr::separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    mutate(gene_id = if ("group_id" %in% names(pick(everything()))) group_id else phenotype_id) |>
    dplyr::select(any_of(c("gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    filter(gene_id %in% genes) |>
    arrange(gene_id, rank)
dim(vestib_qtls_comb)
View(vestib_qtls_comb)

# Write out the data
write_tsv(vestib_qtls_comb, "figures/data_tables/vestibular.comb.qtls.tsv")


###############
## GTEx QTLs ##
###############

## GTEx top associations per gene per modality
gtex_assoc <- crossing(tissue = tissues, modality = modalities) |>
    reframe({
        fname <- str_glue("Pantry/pheast_no_outlier/output/{tissue}/qtl/{modality}.cis_qtl.txt.gz")
        
        if (!file.exists(fname)) return(tibble())
        
        df <- read_tsv(fname, col_types = cols(.default = col_character()))
        
        if ("group_id" %in% names(df)) {
            df <- df |> mutate(gene_id = group_id) |> dplyr::select(-group_id)
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
write_tsv(gtex_assoc, "figures/data_tables/gtex.sep.assoc.tsv.gz")


## GTEx QTLs (mapped separately per modality)
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
write_tsv(gtex_qtls_sep, "figures/data_tables/gtex.sep.qtls.tsv.gz")


## GTEx QTLs (combined-modality mapping)
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
write_tsv(gtex_qtls_comb, "figures/data_tables/gtex.comb.qtls.tsv.gz")


###################
## Geuvadis TWAS ##
###################

## Vestibualr heritability/r2/pvals
vestib_hsq <- tibble(modality = modalities) |>
    reframe(
        read_tsv(str_glue("Pantry/pheast_no_outlier/intermediate/vestibular/twas/{modality}.profile"),
                    col_types = "cccccccc-cccc-c"),
        .by = modality
    ) |>
    dplyr::rename(phenotype_id = id) |>
    mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
            .before = 2) |>
    filter(gene_id %in% genes)
dim(vestib_hsq)
View(vestib_hsq)

# Write out the data
write_tsv(vestib_hsq, "figures/data_tables/vestibular.hsq.tsv")


## Vestibualr TWAS
vestib_twas <- crossing(modality = modalities,
                        trait = traits) |>
    reframe(
        read_tsv(str_glue("Pantry/pheast_no_outlier/twas/vestibular/output/{modality}/{modality}.{trait}.tsv"),
                    col_types = "c---c--------c--ccccccc"),
        .by = c(modality, trait)
    ) |>
    filter(!is.na(TWAS.P)) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
    filter(gene_id %in% genes) |>
    dplyr::select(trait, gene_id, modality, phenotype_id = ID, HSQ:COLOC.PP4) |>
    arrange(trait, gene_id, modality, phenotype_id)


# Write out the data
write_tsv(vestib_twas, "figures/data_tables/vestibular.twas.tsv.gz")


###############
## GTEx TWAS ##
###############

## GTEx heritability/r2/pvals
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
write_tsv(gtex_hsq, "figures/data_tables/gtex.hsq.tsv.gz")



## GTEx TWAS (significant hits only)
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
write_tsv(gtex_twas, "figures/data_tables/gtex.twas_hits.tsv.gz")


### END