### All TWAS tests (every trait x modality x phenotype with a result) in one table, protein-coding and lncRNA genes
### only (Fig 4c). TWAS.Z of enet/top1 models is sign-corrected for the old-weight modalities as in annotate_twas_hits.R.
### Code from the Pantry paper's assemble_results.R.
### Run from analyses/06_twas:
###   Rscript scripts/combine_twas_all_tests.R [dataset=vestibular] [old_weight_modalities=alt_polyA,...,stability]
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "vestibular"
old_mods <- if (length(args) >= 2) args[2] else "alt_polyA,alt_TSS,expression,isoforms,splicing,stability"
old_mods <- if (old_mods == "none") character(0) else strsplit(old_mods, ",")[[1]]
modalities <- c("alt_polyA", "alt_TSS", "expression", "isoforms", "RNA_editing", "splicing", "stability")
traits <- read_tsv("gwas_metadata.tsv", col_types = cols(trait = "c", .default = "-")) |> pull()

genes <- rtracklayer::import("../../resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene", gene_biotype %in% c("protein_coding", "lncRNA")) |>
    pull(gene_id)

twas <- crossing(modality = modalities, trait = traits) |>
    reframe(read_tsv(str_glue("results/{ds}/output/{modality}/{modality}.{trait}.tsv"),
                     col_types = "c---c--------c--ccccccc"),
            .by = c(modality, trait)) |>
    filter(!is.na(TWAS.P)) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+"),
           TWAS.Z = as.numeric(TWAS.Z) * if_else(modality %in% old_mods & MODEL %in% c("enet", "top1"), -1, 1)) |>
    filter(gene_id %in% genes) |>
    dplyr::select(trait, gene_id, modality, phenotype_id = ID, HSQ:COLOC.PP4) |>
    arrange(trait, gene_id, modality, phenotype_id)

out <- str_glue("results/{ds}/output/twas_all_tests.tsv.gz")
write_tsv(twas, out)
cat("written:", out, "(", nrow(twas), "tests; TWAS.Z flipped for enet/top1 of", paste(old_mods, collapse = ", "), ")\n")
