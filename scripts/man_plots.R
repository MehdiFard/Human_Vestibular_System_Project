### Making manhattan plots

rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(qqman)


### Manhattan plot - hg38 result
library(qqman)
gwas_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")
gwas_cln_hg38 <- gwas_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))

# Threshold
fdr_threshold <- 0.05/dim(gwas_hg38)[1]

# Manhattan plot
pdf("figures/man_plot_meta_balance_hg38.pdf", height = 6, width = 10)

manhattan(
    gwas_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()

gwas_cln_hg38_sig <- gwas_cln_hg38 %>%
    dplyr::filter(CHR %in% c(17, 18, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_cln_hg38_sig)
head(gwas_cln_hg38_sig)
fwrite(gwas_cln_hg38_sig, "figures/data_tables/balance_meta_top_hits_hg38.tsv", sep = "\t")



###----------- Manhattan plot - hg38 result
gwas_hg19 <- fread("data/gwas_raw_sumstats/balance_meta1.tbl.sorted.filtered.gz")

# Clean and filter the GWAS data
gwas_cln_hg19 <- gwas_hg19 %>%
    select(CHR = Chromosome, BP = Position, P = "P-value", SNP = MarkerName) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))

# Threshold
fdr_threshold <- 0.05/dim(gwas_hg19)[1]

# Manhattan plot
pdf("figures/man_plot_meta_balance_hg19.pdf", height = 6, width = 10)

manhattan(
    gwas_cln_hg19, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()


###----------- Manhattan plot - meta_Eur result - hg38
gwas_meta_eur_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta_eur_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_meta_eur_hg38)
dim(gwas_meta_eur_hg38)

# Clean and filter the GWAS data
gwas_meta_eur_cln_hg38 <- gwas_meta_eur_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_meta_eur_cln_hg38)

# Threshold
fdr_threshold <- 0.05/dim(gwas_meta_eur_hg38)[1]

# Manhattan plot
pdf("figures/man_plot_balance_meta_eur_hg38.pdf", height = 6, width = 10)

manhattan(
    gwas_meta_eur_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()

gwas_meta_eur_sig_hg38 <- gwas_meta_eur_cln_hg38 %>%
    dplyr::filter(CHR %in% c(17, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_meta_eur_sig_hg38)
head(gwas_meta_eur_sig_hg38)
fwrite(gwas_meta_eur_sig_hg38, "figures/data_tables/balance_meta_eur_top_hits_hg38.tsv", sep = "\t")



###----------- Manhattan plot - UKBB result - hg19
gwas_ukb <- fread("data/gwas_raw_sumstats/Royce_files/UKB_EUR_balance.txt.gz")
head(gwas_ukb)

# Clean and filter the GWAS data
gwas_ukb_cln_hg19 <- gwas_ukb %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_ukb_cln_hg19)

# Threshold
fdr_threshold <- 0.05/dim(gwas_ukb)[1]

# Manhattan plot
pdf("figures/man_plot_balance_ukb_hg19.pdf", height = 6, width = 10)

manhattan(
    gwas_ukb_cln_hg19, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()

gwas_ukb_cln_sig <- gwas_ukb_cln_hg19 %>%
    dplyr::filter(CHR %in% c(3, 13, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_ukb_cln_sig)
head(gwas_ukb_cln_sig)
fwrite(gwas_ukb_cln_sig, "figures/data_tables/balance_ukb_top_hits.tsv", sep = "\t")


###----------- Manhattan plot - UKBB result - hg38
gwas_ukb_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/UKB_EUR_balance_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_ukb_hg38)

# Clean and filter the GWAS data
gwas_ukb_cln_hg38 <- gwas_ukb_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_ukb_cln_hg38)

# Threshold
fdr_threshold <- 0.05/dim(gwas_ukb_hg38)[1]

# Manhattan plot
pdf("figures/man_plot_balance_ukb_hg38.pdf", height = 6, width = 10)

manhattan(
    gwas_ukb_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()

gwas_ukb_cln_hg38_sig <- gwas_ukb_cln_hg38 %>%
    dplyr::filter(CHR %in% c(13, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_ukb_cln_hg38_sig)
head(gwas_ukb_cln_hg38_sig)
fwrite(gwas_ukb_cln_hg38_sig, "figures/data_tables/balance_ukb_top_hits_hg38.tsv", sep = "\t")

