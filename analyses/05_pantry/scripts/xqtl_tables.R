### xQTL summary tables for the paper figures (Fig 4, Supp Fig 3) from the Pantry phenotypes and the pheast xQTLs.
### Protein-coding and lncRNA genes only (Ensembl 110). Code from the Pantry paper's figure1.R / assemble_results.R.
### Run from analyses/05_pantry:  Rscript scripts/xqtl_tables.R [dataset=vestibular]
### Writes results/xqtl_summary/:
###   <ds>.phenotypes.tsv.gz  phenotypes per modality with their gene and biotype (Fig 4a, Supp Fig 3a)
###   <ds>.sep.assoc.tsv      top association per gene and modality, separate mapping (genes tested: enrichment background)
###   <ds>.comb.qtls.tsv      conditionally independent xQTLs of the cross-modality mapping (Fig 4b, Supp Fig 3b-e)
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "vestibular"
out <- "results/xqtl_summary"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
modalities <- c("alt_polyA", "alt_TSS", "expression", "isoforms", "RNA_editing", "splicing", "stability")

genes <- rtracklayer::import("../../resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene", gene_biotype %in% c("protein_coding", "lncRNA")) |>
    dplyr::select(gene_id, gene_biotype)

## ---- phenotypes per modality (expression and stability are gene-level; the others map to genes via their groups file)
get_phenotypes <- function(modality) {
    df <- read_tsv(str_glue("phenotyping/output/{ds}/{modality}.bed.gz"),
                   col_types = cols(phenotype_id = "c", .default = "-"))
    if (modality %in% c("expression", "stability")) return(mutate(df, gene_id = phenotype_id))
    groups <- read_tsv(str_glue("phenotyping/output/{ds}/{modality}.phenotype_groups.txt"),
                       col_types = "cc", col_names = c("phenotype_id", "gene_id"))
    left_join(df, groups, by = "phenotype_id")
}
phenos <- tibble(modality = modalities) |>
    reframe(get_phenotypes(modality), .by = modality) |>
    inner_join(genes, by = "gene_id", relationship = "many-to-one")
write_tsv(phenos, file.path(out, str_glue("{ds}.phenotypes.tsv.gz")))
cat("phenotypes per modality:\n"); print(table(phenos$modality))

## ---- top association per gene and modality (separate mapping per modality)
read_ids <- function(fname, df) {
    if ("group_id" %in% colnames(df)) return(df |> mutate(gene_id = group_id) |> dplyr::select(-group_id))
    if ("phenotype_id" %in% colnames(df)) return(df |> mutate(gene_id = phenotype_id))
    stop("Neither group_id nor phenotype_id found in ", fname)
}
sep_assoc <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("pheast/output/{ds}/qtl/{modality}.cis_qtl.txt.gz")
        read_ids(fname, read_tsv(fname, col_types = cols(.default = col_character())))
    }, .by = modality) |>
    filter(gene_id %in% genes$gene_id) |>
    relocate(gene_id, .before = 1) |>
    arrange(gene_id, modality) |>
    mutate(qval = as.numeric(qval))
write_tsv(sep_assoc, file.path(out, str_glue("{ds}.sep.assoc.tsv")))
cat("xGenes per modality (separate mapping, q < 0.05):\n")
print(sep_assoc |> filter(qval < 0.05) |> distinct(modality, gene_id) |> dplyr::count(modality))

## ---- independent xQTLs, cross-modality mapping (phenotype_id = "<modality>:<phenotype>")
comb_qtls <- read_tsv(str_glue("pheast/output/{ds}/qtl/cross_modality.cis_independent_qtl.txt.gz"),
                      col_types = cols(.default = col_character(), rank = col_integer())) |>
    separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    mutate(gene_id = if ("group_id" %in% names(pick(everything()))) group_id else phenotype_id) |>
    dplyr::select(any_of(c("gene_id", "rank", "modality", "phenotype_id", "variant_id", "tss_distance", "pval_beta"))) |>
    filter(gene_id %in% genes$gene_id) |>
    arrange(gene_id, rank)
write_tsv(comb_qtls, file.path(out, str_glue("{ds}.comb.qtls.tsv")))
cat("xQTLs per modality (cross-modality mapping):\n"); print(table(comb_qtls$modality))
cat("written:", out, "\n")
