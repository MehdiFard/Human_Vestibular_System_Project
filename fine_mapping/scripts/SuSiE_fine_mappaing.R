### Runing fine mapping using susieR
rm(list = ls())

library(susieR)
library(dplyr)
library(data.table)
library(BEDMatrix)
library(ggplot2)


# Preparing LD matrix (R) using BEDMatrix
bed <- "Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bed"  # needs geno.bim and geno.fam in same folder

G <- BEDMatrix(bed)         # samples x SNPs (lazy / memory-mapped)
bim <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bim", sep = "\t", stringsAsFactors = FALSE)
fam <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.fam", sep = "\t", stringsAsFactors = FALSE)
colnames(bim) <- c("chr","id","cm","pos","a1","a2")

# choose locus
chr <- 19
start_bp <- 23e6
end_bp <- 24e6

# get indices of SNPs in the locus
idx <- which(bim$chr == chr & bim$pos >= start_bp & bim$pos <= end_bp)

# extract locus genotypes (can be large; still manageable for typical cis windows)
X <- as.matrix(G[, idx])     # numeric 0/1/2 with possible NAs
dim(X)
X[1:20, 1:3]

# handle missing values (choose ONE)
# Option 1: mean-impute per SNP
for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j])
    if (any(nas)) X[nas, j] <- mean(X[, j], na.rm = TRUE)
}

# compute LD correlation matrix (SNP x SNP)
R <- cor(X)
dim(R)
R[1:5, 1:5]


# basic sanity
R <- (R + t(R)) / 2
diag(R) <- 1
dim(R)
R[1:5, 1:5]


###--------------------- SusieR
## Read GWAS sumstats
ss <- read.table("gwas_sumstats_prep/output/balance_meta1_hg38_twas_input.sumstats.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(ss)

## Get the SNP IDs in the LD matrix order
snps_ld <- bim$id[idx]   # same order as columns in X and in R
length(snps_ld)

## Match sumstats to LD SNP order
m <- match(snps_ld, ss$SNP)
keep <- which(!is.na(m))


## subset to SNPs present in both
snps_ld2 <- snps_ld[keep]
R2 <- R[keep, keep, drop = FALSE]
ss2 <- ss[m[keep], ]
z <- ss2$Z

## Basic sanity checks
stopifnot(length(z) == ncol(R2), ncol(R2) == nrow(R2))

# (optional) check no NA
stopifnot(!any(is.na(z)), !any(is.na(R2)))

## Choose effective N
# If N varies by SNP, use median (common practice for susie_rss)
n_eff <- as.integer(median(ss2$N, na.rm = TRUE))

## Stabilize LD if needed (ridge)
R2 <- (R2 + t(R2)) / 2
diag(R2) <- 1
R2 <- R2 + diag(1e-5, nrow(R2))  # helps if near-singular


## Run SuSiE
fit <- susie_rss(
    z = z,
    R = R2,
    n = n_eff,
    L = 10,            # max signals; try 5-20 depending on locus
    max_iter = 200
)


## Results
fit$converged
fit$sets$cs          # credible sets (indices)
fit$sets$purity      # purity metrics per CS (if computed)
pip <- fit$pip


## list of integer indices (one vector per CS)
cs_idx <- fit$sets$cs
cs_snps <- lapply(cs_idx, function(ii) snps_ld2[ii]) # SNP IDs per credible set

## Make a result table for credible sets
cs_df <- do.call(rbind, lapply(seq_along(cs_idx), function(k) {
    ii <- cs_idx[[k]]
    data.frame(
        cs_id = k,
        SNP   = snps_ld2[ii],
        Z     = z[ii],
        PIP   = pip[ii],
        stringsAsFactors = FALSE
    )
}))

# sort by cs_id and PIP within each CS
cs_df <- cs_df[order(cs_df$cs_id, -cs_df$PIP), ]

# sanity check: Z in cs_df matches z by index
stopifnot(
    all(cs_df$Z == z[ match(cs_df$SNP, snps_ld2) ]))

# Spot-check
cs_df[1:10, c("cs_id", "SNP", "Z")]
z[ match(cs_df$SNP[1:10], snps_ld2) ]


### Adding rsIDs
gwas_res <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")

# Merge
cs_df_rsid <- merge(cs_df, gwas_res[, c("SNP", "rsID", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "P")], by = "SNP", sort = F)
head(cs_df_rsid)

# Save
fwrite(cs_df_rsid, file = "fine_mapping/results/susie_results_balance_meta1.tsv", sep = "\t")


### Plotting the results
# Locus plot: position vs -log10(P) (shape = CS, size = PIP
cs_df_rsid$logp <- -log10(cs_df_rsid$P)
ggplot(cs_df_rsid, aes(x = BP, y = logp)) +
    geom_point(aes(color = factor(cs_id), size = PIP), alpha = 0.8) +
    theme_classic() +
    labs(
        x = "Position (bp)",
        y = "-log10(P)",
        color = "Credible set",
        size = "PIP"
    )
ggsave("fine_mapping/plots/plot1_susie_results_balance_meta1_locus_plot.png", width = 8, height = 8)


# PIP vs |Z|, label the top variant per CS (by PIP)
top_per_cs <- cs_df_rsid %>%
    group_by(cs_id) %>%
    slice_max(PIP, n = 1, with_ties = FALSE) %>%
    ungroup()

ggplot(cs_df_rsid, aes(x = abs(Z), y = PIP, color = factor(cs_id))) +
    geom_point(alpha = 0.7) +
    geom_text(data = top_per_cs, aes(label = rsID), vjust = -0.6, show.legend = FALSE) +
    theme_classic() +
    labs(x = "|Z|", y = "PIP", color = "CS")

ggsave("fine_mapping/plots/plot2_susie_results_balance_meta1_pip_vs_z.png", width = 8, height = 8)


### END


### Run susie on QTL mapping
qtls <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/cross_modality.cis_qtl.txt.gz")
qtls <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/isoforms.cis_qtl.txt.gz")
dim(qtls)

## Get the SNP IDs in the LD matrix order
snps_ld <- bim$id[idx]   # same order as columns in X and in R
length(snps_ld)


## Match sumstats to LD SNP order
m <- match(snps_ld, qtls$variant_id)
length(m)
keep <- which(!is.na(m))
length(keep)

## subset to SNPs present in both
snps_ld2 <- snps_ld[keep]
R2 <- R[keep, keep, drop = FALSE]
qtls2 <- qtls[m[keep], ]
qtls2$z <- qtls2$slope / qtls2$slope_se
z <- qtls2$z
View(qtls2)

## Basic sanity checks
stopifnot(length(z) == ncol(R2), ncol(R2) == nrow(R2))

# (optional) check no NA
stopifnot(!any(is.na(z)), !any(is.na(R2)))

## Choose effective N
# If N varies by SNP, use median (common practice for susie_rss)
n_eff <- as.integer(median(qtls2$ma_samples, na.rm=TRUE))

## Stabilize LD if needed (ridge)
R2 <- (R2 + t(R2)) / 2
diag(R2) <- 1
R2 <- R2 + diag(1e-5, nrow(R2))  # helps if near-singular
dim(R2)

## Run SuSiE
fit <- susie_rss(
    z = z,
    R = R2,
    n = n_eff,
    L = 10,            # max signals; try 5-20 depending on locus
    max_iter = 200
)


## Results
fit$converged
fit$sets$cs          # credible sets (indices)
fit$sets$purity      # purity metrics per CS (if computed)
pip <- fit$pip


## list of integer indices (one vector per CS)
cs_idx <- fit$sets$cs
cs_snps <- lapply(cs_idx, function(ii) snps_ld2[ii]) # SNP IDs per credible set

## Make a result table for credible sets
cs_df <- do.call(rbind, lapply(seq_along(cs_idx), function(k) {
    ii <- cs_idx[[k]]
    data.frame(
        cs_id = k,
        SNP   = snps_ld2[ii],
        Z     = z[ii],
        PIP   = pip[ii],
        stringsAsFactors = FALSE
    )
}))

# sort by cs_id and PIP within each CS
cs_df <- cs_df[order(cs_df$cs_id, -cs_df$PIP), ]

# sanity check: Z in cs_df matches z by index
stopifnot(
    all(cs_df$Z == z[ match(cs_df$SNP, snps_ld2) ]))

# Spot-check
cs_df[1:10, c("cs_id", "SNP", "Z")]
z[ match(cs_df$SNP[1:10], snps_ld2) ]
View(cs_df)

### Adding rsIDs
gwas_res <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")

# Merge
cs_df_rsid <- merge(cs_df, gwas_res[, c("SNP", "rsID", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "P")], by = "SNP", sort = F)
dim(cs_df_rsid)
View(cs_df_rsid)

# Save
fwrite(cs_df_rsid, file = "fine_mapping/results/susie_results_QTLs_cross_modalities.tsv", sep = "\t")
fwrite(cs_df_rsid, file = "fine_mapping/results/susie_results_QTLs_isoforms.tsv", sep = "\t")


### Plotting the results
# Locus plot: position vs -log10(P) (shape = CS, size = PIP
cs_df_rsid$logp <- -log10(cs_df_rsid$P)
ggplot(cs_df_rsid, aes(x = BP, y = logp)) +
    geom_point(aes(color = factor(cs_id), size = PIP), alpha = 0.8) +
    theme_classic() +
    labs(
        x = "Position (bp)",
        y = "-log10(P)",
        color = "Credible set",
        size = "PIP"
    )
ggsave("fine_mapping/plots/plot3_susie_results_QTLs_cross_modalities_locus_plot.png", width = 8, height = 8)


# PIP vs |Z|, label the top variant per CS (by PIP)
top_per_cs <- cs_df_rsid %>%
    group_by(cs_id) %>%
    slice_max(PIP, n = 1, with_ties = FALSE) %>%
    ungroup()

ggplot(cs_df_rsid, aes(x = abs(Z), y = PIP, color = factor(cs_id))) +
    geom_point(alpha = 0.7) +
    geom_text(data = top_per_cs, aes(label = rsID), vjust = -0.6, show.legend = FALSE, size = 3) +
    theme_classic() +
    labs(x = "|Z|", y = "PIP", color = "CS")

ggsave("fine_mapping/plots/plot4_susie_results_QTLs_cross_modalities_pip_vs_z.png", width = 8, height = 8)


### END
