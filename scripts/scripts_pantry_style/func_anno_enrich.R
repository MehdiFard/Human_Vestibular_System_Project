rm(list = ls())
library(GenomicRanges)
library(tidyverse)
library(tibble)

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
        left_join(count_bg,
                    by = c("tissue", "modality", "category"),
                    relationship = "one-to-one") |>
        mutate(frac_qtls = count_qtls / total_qtls,
                frac_bg = count_bg / total_bg,
                log2_enrich = log2(frac_qtls / frac_bg))
}



# All genes tested for QTLs per tissue-modality, to use set of cis-window variants as background
# Use separately-mapped data since those actually give the genes tested per modality

# tested <- read_tsv("figures/source_data/gtex.sep.assoc.tsv.gz", col_types = "ccc-----")
tested <- read_tsv("figures/data_tables/vestibular.sep.assoc.tsv", col_types = "ccc-----")

# Adding tissue name
tested$tissue <- as.character(rep("VESTIB", nrow(tested)))

# Adding "_b38" suffix to variant IDs to match those in annotation file
tested <- tested |>
    dplyr::mutate(
        variant_id = paste0("chr" ,variant_id, "_b38"))
dim(tested)
View(tested[1:20, ])

# Getting gene annotation from GTF file
genes <- rtracklayer::import("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start)) |>
    select(gene_id, chrom = seqnames, TSS)


# Download annotation file from: https://storage.googleapis.com/adult-gtex/references/v8/reference-tables/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz
anno <- read_tsv("src/GTEx/wgs_vep_annotaion/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz",
                    col_types = cols(SNP = "c", .default = "l"))
class(anno)
dim(anno)
View(anno[1:20, ])

anno_snps <- anno$SNP
library(purrr)
anno <- colnames(anno)[-1] |>
    set_names() |>
    map(\(category) anno$SNP[anno[[category]]])
class(anno)
names(anno)

library(tibble)
count_bg <- background_counts(tested, genes, anno, anno_snps)
dim(count_bg)
View(count_bg)

qtls_comb <- read_tsv("figures/data_tables/vestibular.comb.qtls.tsv", col_types = "ccccccc-")
dim(qtls_comb)
View(qtls_comb[1:20, ])

# Adding tissue name
qtls_comb$tissue <- as.character(rep("VESTIB", nrow(qtls_comb)))

# Adding "chr" and "_b38" suffix to variant IDs to match those in annotation file
qtls_comb <- qtls_comb |>
    dplyr::mutate(
        variant_id = paste0("chr" ,variant_id, "_b38"))

# Run enrichment analysis for combined QTLs
enrich_comb <- enrichment_table(qtls_comb, anno, count_bg)
dim(enrich_comb)
View(enrich_comb[1:20, ])

# Save enrichment results
write_tsv(enrich_comb, "figures/data_tables/enrich.vestibular.comb.tsv")


qtls_sep <- read_tsv("figures/data_tables/vestibular.sep.qtls.tsv", col_types = "cccccccc")
dim(qtls_sep)
View(qtls_sep[1:20, ])

# Adding tissue name
qtls_sep$tissue <- as.character(rep("VESTIB", nrow(qtls_sep)))

# Adding "chr" and "_b38" suffix to variant IDs to match those in annotation file
qtls_sep <- qtls_sep |>
    dplyr::mutate(
        variant_id = paste0("chr" ,variant_id, "_b38"))


enrich_sep <- enrichment_table(qtls_sep, anno, count_bg)
dim(enrich_sep)
View(enrich_sep[1:20, ])

write_tsv(enrich_sep, "figures/data_tables/enrich.vestibular.sep.tsv")
    
# Compare enrichment in two sets of QTLs per tissue-modality. Looking at genes
# found at least once in any separate and any combined mapping in any modality
# in a tissue:
# 1. Those in genes found only using separate mapping for the modality
# 2. Those in genes found using both methods for the modality
qtls_venn <- bind_rows(
    qtls_comb |> mutate(method = "combined"),
    qtls_sep |> mutate(method = "separate")
) |>
    group_by(tissue, gene_id) |>
    filter(length(unique(method)) == 2) |> # Include only genes found with both methods to estimate "consolidation" of redundant QTLs
    ungroup() |>
    group_by(tissue, modality, gene_id) |>
    mutate(methods = str_c(sort(unique(method)), collapse = "_")) |>
    ungroup()

dim(qtls_venn)
View(qtls_venn[1:20, ])
table(qtls_venn$methods)

enrich_sep_only <- qtls_venn |>
    filter(methods == "separate") |>
    enrichment_table(anno, count_bg)
dim(enrich_sep_only)
View(enrich_sep_only[1:20, ])

# Save enrichment results for genes found only in separate mapping
write_tsv(enrich_sep_only, "figures/data_tables/enrich.vestibular.sep_only.tsv")


enrich_in_both <- qtls_venn |>
    filter(methods == "combined_separate") |>
    enrichment_table(anno, count_bg)

dim(enrich_in_both)
View(enrich_in_both[1:20, ])

# Save enrichment results for genes found in both separate and combined mapping
write_tsv(enrich_in_both, "figures/data_tables/enrich.vestibular.in_sep_and_comb.tsv")
