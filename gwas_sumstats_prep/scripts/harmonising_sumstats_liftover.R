#!/usr/bin/env Rscript

# Do NOT setwd() when called by Snakemake.
suppressPackageStartupMessages({
    library(data.table)
    library(MungeSumstats)
    library(dplyr)
    library(readr)
    library(optparse)
    library(stringr)
})

option_list <- list(
    make_option("--sumstat",   type="character", help="Input GWAS summary stats file"),
    make_option("--ref_build", type="character", help="Reference genome build (GRCh37 or GRCh38)"),
    make_option("--out",       type="character", help="Output path for formatted sumstats (.tsv.gz)")
)
opt <- parse_args(OptionParser(option_list=option_list))

cat("Running harmonising_sumstats_liftover.R\n")
cat("Sumstat:",   opt$sumstat,   "\n")
cat("Ref build:", opt$ref_build, "\n")
cat("Output:",    opt$out,       "\n")

sumstats_path <- opt$sumstat
out_file      <- opt$out
build         <- opt$ref_build

# Allow the chain path to be relative to project root; adjust if needed.
chain_path <- "src/hg19ToHg38.over.chain"

# Run MungeSumstats
if (build == "GRCh37") {
    gwas_sumstat <- format_sumstats(
        path                  = sumstats_path,
        ref_genome            = "GRCh37",
        convert_ref_genome    = "GRCh38",
        local_chain           = chain_path,
        nThread               = 16,
        force_new             = TRUE,
        bi_allelic_filter     = FALSE,
        flip_frq_as_biallelic = TRUE,
        rmv_chr               = c("X", "Y", "MT"),
        save_path             = out_file
    )
} else {
    gwas_sumstat <- format_sumstats(
        path                  = sumstats_path,
        ref_genome            = "GRCh38",
        nThread               = 16,
        force_new             = TRUE,
        bi_allelic_filter     = FALSE,
        flip_frq_as_biallelic = TRUE,
        rmv_chr               = c("X", "Y", "MT"),
        save_path             = out_file
    )
}

# Getting trait name
trait <- sub("_formatted_hg38.tsv.gz", "", basename(out_file))

# Add chr_pos_ref_alt ID and write companion file
if(trait %in% c("cerebellar_volume", "cerebral_volume", "subcortical_volume")){
    sum_cln <- fread(out_file) %>%
        dplyr::rename(rsID = SNP) %>%
        mutate(SNP = str_c(CHR, BP, A1, A2, sep = "_")) %>%
        mutate(Z = if_else(is.na(SE) | SE == 0, NA_real_, BETA / SE)) %>%
        select(SNP, everything())  # keep rsID if you want; remove with `, -rsID` if not

} else {

    sum_cln <- fread(out_file)

    # If Z_STAT exists and Z doesn't, rename it to Z
    if ("Z_STAT" %in% names(sum_cln) && !"Z" %in% names(sum_cln)) {
        sum_cln <- sum_cln %>% dplyr::rename(Z = Z_STAT)
    }

    sum_cln <- sum_cln %>%
        dplyr::rename(rsID = SNP) %>%
        mutate(SNP = str_c(CHR, BP, A1, A2, sep = "_")) %>%
        select(SNP, everything())  
}

# Write out and gzip
out_path_cln <- sub("\\.tsv\\.gz$", "_ChrPosId.tsv.gz", out_file)
fwrite(sum_cln, file = out_path_cln, sep = "\t", quote = FALSE, compress = "gzip")


### END