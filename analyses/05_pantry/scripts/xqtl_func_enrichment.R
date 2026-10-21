### Functional annotation enrichment of the cross-modality xQTLs (Supp Fig 3e). Background: all variants in the cis
### windows (TSS +- 1 Mb) of the genes tested per modality (separate mapping). Annotation: GTEx v8 WGS VEP / regulatory
### features (src/gtex_wgs_vep_annotation, from https://storage.googleapis.com/adult-gtex/references/v8/reference-tables/
### WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz). Code from the Pantry paper's func_anno_enrich.R.
### Run from analyses/05_pantry (after xqtl_tables.R):  Rscript scripts/xqtl_func_enrichment.R [dataset=vestibular]
suppressPackageStartupMessages({ library(GenomicRanges); library(tidyverse) })
args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "vestibular"
out <- "results/xqtl_summary"

snps_in_windows <- function(genes, snps, snps_rng) {
    cis_rng <- with(genes, GRanges(chrom, IRanges(TSS - 1e6, TSS + 1e6)))
    snps[countOverlaps(snps_rng, cis_rng) > 0]
}

background_snps <- function(gene_ids, genes, snps, snps_rng) {
    genes |>
        filter(gene_id %in% gene_ids) |>
        snps_in_windows(snps, snps_rng)
}

count_in_snps <- function(snps, anno) {
    map_int(anno, \(anno_snps) sum(snps %in% anno_snps)) |>
        enframe(name = "category", value = "count")
}

background_counts <- function(tested, genes, anno, anno_snps) {
    snp_info <- tibble(SNP = anno_snps) |>
        mutate(chrom = str_extract(SNP, "^chr([^_]+)_", group = 1),
               pos = str_extract(SNP, "^[^_]+_([^_]+)_", group = 1) |> as.integer())
    snps_rng <- with(snp_info, GRanges(chrom, IRanges(pos, pos)))
    tested |>
        reframe({
            bg_snps <- background_snps(gene_id, genes, snp_info$SNP, snps_rng)
            count_in_snps(bg_snps, anno) |>
                dplyr::rename(count_bg = count) |>
                mutate(total_bg = length(bg_snps))
        }, .by = c(tissue, modality))
}

enrichment_table <- function(qtls, anno, count_bg) {
    count_qtls <- qtls |>
        distinct(tissue, modality, variant_id) |>
        reframe(
            count_in_snps(variant_id, anno) |>
                dplyr::rename(count_qtls = count) |>
                mutate(total_qtls = length(variant_id)),
            .by = c(tissue, modality)
        )
    count_qtls |>
        left_join(count_bg, by = c("tissue", "modality", "category"), relationship = "one-to-one") |>
        mutate(frac_qtls = count_qtls / total_qtls,
               frac_bg = count_bg / total_bg,
               log2_enrich = log2(frac_qtls / frac_bg))
}

## variant IDs as in the annotation file: chr<chrom>_<pos>_<ref>_<alt>_b38
to_gtex_id <- function(df) df |> mutate(tissue = "VESTIB", variant_id = paste0("chr", variant_id, "_b38"))

tested <- read_tsv(file.path(out, str_glue("{ds}.sep.assoc.tsv")), col_types = "ccc-----") |> to_gtex_id()

genes <- rtracklayer::import("../../resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start)) |>
    dplyr::select(gene_id, chrom = seqnames, TSS)

anno <- read_tsv("src/gtex_wgs_vep_annotation/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz",
                 col_types = cols(SNP = "c", .default = "l"))
anno_snps <- anno$SNP
anno <- colnames(anno)[-1] |>
    set_names() |>
    map(\(category) anno$SNP[anno[[category]]])

count_bg <- background_counts(tested, genes, anno, anno_snps)

qtls_comb <- read_tsv(file.path(out, str_glue("{ds}.comb.qtls.tsv")), col_types = "ccccccc-") |> to_gtex_id()
enrich_comb <- enrichment_table(qtls_comb, anno, count_bg)
write_tsv(enrich_comb, file.path(out, str_glue("enrich.{ds}.comb.tsv")))
cat("written:", file.path(out, str_glue("enrich.{ds}.comb.tsv")), "\n")
