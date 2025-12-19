rm(list = ls())
library(GenomicRanges)
library(tidyverse)
library(data.table)


nearest_pos <- function(genes, posn, posns) {
    # For each posn, return closest position in posns for the same gene. `posns`
    # is tibble with columns `gene_id` and `pos`. Returned vector corresponds to
    # `genes` and `posn`, and includes NA for positions with no overlapping
    # feature.
    all_genes <- unique(c(genes, posns$gene_id)) # To avoid GRanges warning about non-shared seqnames
    pos_rng <- GRanges(genes, IRanges(posn, posn), seqinfo = all_genes)
    feat_rng <- with(posns, GRanges(gene_id, IRanges(pos, pos), seqinfo = all_genes))
    hits <- nearest(pos_rng, feat_rng, select = "all") |>
        as_tibble() |>
        group_by(queryHits) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(hitpos = posns$pos[subjectHits])
    tibble(qpos = posn) |>
        mutate(index = 1:n()) |>
        left_join(hits, by = c("index" = "queryHits"), relationship = "one-to-one") |>
        pull(hitpos)
}

feature_rel_pos <- function(genes, pos, features) {
    # For each pos, find overlapping feature for same gene (choose at random if
    # multiple) and return relative position. Returned vector corresponds to
    # `genes` and `pos`, and includes NA for positions with no overlapping
    # feature.
    all_genes <- unique(c(genes, features$gene_id)) # To avoid GRanges warning about non-shared seqnames
    pos_rng <- GRanges(genes, IRanges(pos, pos), seqinfo = all_genes)
    feat_rng <- with(features, GRanges(
        gene_id,
        IRanges(pmin(start, end), pmax(start, end)),
        seqinfo = all_genes
    ))
    hits <- findOverlaps(pos_rng, feat_rng, select = "all") |>
        as_tibble() |>
        group_by(queryHits) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(start = features$start[subjectHits],
                end = features$end[subjectHits])
    tibble(pos = pos) |>
        mutate(index = 1:n()) |>
        left_join(hits, by = c("index" = "queryHits"), relationship = "one-to-one") |>
        mutate(rel_pos = (pos - start) / (end - start)) |>
        pull(rel_pos)
}

anno <- rtracklayer::import("Pantry/phenotyping/input/ref/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    dplyr::rename(chrom = seqnames)

genes <- anno |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start),
            TES = if_else(strand == "-", start, end)) |>
    dplyr::select(gene_id, chrom, TSS, TES)

# Record all exons per isoform in chromosome coordinate order
exons_chrom <- anno |>
    filter(type == "exon") |>
    dplyr::select(gene_id, transcript_id, exon_start = start, exon_end = end, strand, exon_number) |>
    arrange(gene_id, transcript_id, exon_start)

# Oriented and ordered in direction of gene, i.e. for intron on - strand, start > end, just like TSS/TES
introns <- exons_chrom |>
    group_by(transcript_id) |>
    mutate(intron_start = exon_end + 1,
            intron_end = lead(exon_start)) |>
    ungroup() |>
    filter(!is.na(intron_end)) |>
    mutate(start = if_else(strand == '+', intron_start, intron_end),
            end = if_else(strand == '+', intron_end, intron_start)) |>
    arrange(gene_id, transcript_id, exon_number) |>
    dplyr::select(gene_id, transcript_id, start, end)

# Oriented and ordered in direction of gene, i.e. for intron on - strand, start > end, just like TSS/TES
exons <- exons_chrom |>
    mutate(start = if_else(strand == '+', exon_start, exon_end),
            end = if_else(strand == '+', exon_end, exon_start)) |>
    dplyr::select(gene_id, transcript_id, exon_number, start, end) |>
    arrange(gene_id, transcript_id, exon_number)

exon_intron <- exons |>
    group_by(transcript_id) |>
    filter(exon_number != max(exon_number)) |> # Last exon has no exon-intron boundary
    ungroup() |>
    # Gene-oriented, so end pos is always exon-intron boundary
    dplyr::select(gene_id, transcript_id, exon_number, pos = end) |>
    distinct(gene_id, pos)

intron_exon <- exons |>
    group_by(transcript_id) |>
    filter(exon_number > 1) |> # First exon has no intron-exon boundary
    ungroup() |>
    # Gene-oriented, so start pos is always intron-exon boundary
    dplyr::select(gene_id, transcript_id, exon_number, pos = start) |>
    distinct(gene_id, pos)

# Reading SNPs map file to get chr and pos
vestib_bim <- fread("Pantry/pheast_no_outlier/input/vestibular/vestibular_genotype.bim")[,c(1,2,4)]
names(vestib_bim) <- c("chrom", "variant_id", "pos")
head(vestib_bim)

# For combined modalities (cross-modalities results)
qtls <- read_tsv("figures/data_tables/vestibular.comb.qtls.tsv", col_types = "cicccd") |>
    inner_join(vestib_bim, by = "variant_id") |>
    left_join(genes, by = "gene_id") |>
    mutate(
        rel_pos_gene = (pos - TSS) / (TES - TSS),
        rel_pos_TSS = if_else(TSS < TES, pos - TSS, TSS - pos),
        rel_pos_TES = if_else(TSS < TES, pos - TES, TES - pos),
        rel_pos_exon = feature_rel_pos(gene_id, pos, exons),
        rel_pos_intron = feature_rel_pos(gene_id, pos, introns),
        ex_in_bnd = nearest_pos(gene_id, pos, exon_intron),
        in_ex_bnd = nearest_pos(gene_id, pos, intron_exon),
        rel_pos_ex_in_bnd = if_else(TSS < TES, pos - ex_in_bnd, ex_in_bnd - pos),
        rel_pos_in_ex_bnd = if_else(TSS < TES, pos - in_ex_bnd, in_ex_bnd - pos),
    )
stopifnot(all(qtls$chrom.x == qtls$chrom.y))

posns <- qtls |>
    dplyr::select(gene_id, modality, rel_pos_gene, rel_pos_TSS, rel_pos_TES,
            rel_pos_exon, rel_pos_intron, rel_pos_ex_in_bnd, rel_pos_in_ex_bnd)

# Wite it out
write_tsv(posns, "figures/data_tables/vestibular.comb.qtls.rel_pos.tsv")


# For modalities analyzed separately
qtls <- read_tsv("figures/data_tables/vestibular.sep.qtls.tsv", col_types = "cicccd") |>
    inner_join(vestib_bim, by = "variant_id") |>
    left_join(genes, by = "gene_id") |>
    mutate(
        rel_pos_gene = (pos - TSS) / (TES - TSS),
        rel_pos_TSS = if_else(TSS < TES, pos - TSS, TSS - pos),
        rel_pos_TES = if_else(TSS < TES, pos - TES, TES - pos),
        rel_pos_exon = feature_rel_pos(gene_id, pos, exons),
        rel_pos_intron = feature_rel_pos(gene_id, pos, introns),
        ex_in_bnd = nearest_pos(gene_id, pos, exon_intron),
        in_ex_bnd = nearest_pos(gene_id, pos, intron_exon),
        rel_pos_ex_in_bnd = if_else(TSS < TES, pos - ex_in_bnd, ex_in_bnd - pos),
        rel_pos_in_ex_bnd = if_else(TSS < TES, pos - in_ex_bnd, in_ex_bnd - pos),
    )
stopifnot(all(qtls$chrom.x == qtls$chrom.y))

posns <- qtls |>
    dplyr::select(gene_id, modality, rel_pos_gene, rel_pos_TSS, rel_pos_TES,
            rel_pos_exon, rel_pos_intron, rel_pos_ex_in_bnd, rel_pos_in_ex_bnd)

write_tsv(posns, "figures/data_tables/vestibular.sep.qtls.rel_pos.tsv")

### END