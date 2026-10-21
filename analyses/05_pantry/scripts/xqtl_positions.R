### Position of each cross-modality xQTL relative to its gene (Supp Fig 3d): normalised to gene length (0 = gene start,
### 1 = gene end) and in bp from the TSS and TES, in the gene's orientation. Code from the Pantry paper's
### QTL_pos_norm_gene_data.R (the exon / intron positions it also computed are not used and were dropped).
### Run from analyses/05_pantry (after xqtl_tables.R):  Rscript scripts/xqtl_positions.R [dataset=vestibular]
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "vestibular"
out <- "results/xqtl_summary"

genes <- rtracklayer::import("../../resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf") |>
    as_tibble() |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start),
           TES = if_else(strand == "-", start, end)) |>
    dplyr::select(gene_id, chrom = seqnames, TSS, TES)

bim <- fread(str_glue("pheast/input/{ds}/{ds}_genotype.bim"))[, c(1, 2, 4)]
names(bim) <- c("chrom", "variant_id", "pos")

qtls <- read_tsv(file.path(out, str_glue("{ds}.comb.qtls.tsv")), col_types = "cicccd") |>
    inner_join(bim, by = "variant_id") |>
    left_join(genes, by = "gene_id") |>
    mutate(rel_pos_gene = (pos - TSS) / (TES - TSS),
           rel_pos_TSS = if_else(TSS < TES, pos - TSS, TSS - pos),
           rel_pos_TES = if_else(TSS < TES, pos - TES, TES - pos))
stopifnot(all(qtls$chrom.x == qtls$chrom.y))

qtls |>
    dplyr::select(gene_id, modality, rel_pos_gene, rel_pos_TSS, rel_pos_TES) |>
    write_tsv(file.path(out, str_glue("{ds}.comb.qtls.rel_pos.tsv")))
cat("written:", file.path(out, str_glue("{ds}.comb.qtls.rel_pos.tsv")), "(", nrow(qtls), "xQTLs )\n")
