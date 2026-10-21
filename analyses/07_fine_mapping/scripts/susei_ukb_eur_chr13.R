### Runing fine mapping using susieR
rm(list = ls())

library(susieR)
library(dplyr)
library(data.table)
library(BEDMatrix)
library(ggplot2)


bed_13 <- "resources/ldref_1kg_hg38/1kg_ldref_hg38_chr13.bed"  # needs geno.bim and geno.fam in same folder

G <- BEDMatrix(bed_13)         # samples x SNPs (lazy / memory-mapped)
bim <- read.table("resources/ldref_1kg_hg38/1kg_ldref_hg38_chr13.bim", sep = "\t", stringsAsFactors = FALSE)
fam <- read.table("resources/ldref_1kg_hg38/1kg_ldref_hg38_chr13.fam", sep = "\t", stringsAsFactors = FALSE)
colnames(bim) <- c("chr","id","cm","pos","a1","a2")
head(bim)

# choose locus
chr <- 13
start_bp <- 53045982 - 100000
end_bp <- 53472473 + 100000

# get indices of SNPs in the locus
idx <- which(bim$chr == chr & bim$pos >= start_bp & bim$pos <= end_bp)
length(idx)

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
ss <- read.table("analyses/04_gwas_sumstats/results/twas_input/UKB_EUR_balance_hg38_twas_input.sumstats.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
    max_iter = 1000
)


## Results
fit$converged
fit$sets$cs          # credible sets (indices)
fit$sets$purity      # purity metrics per CS (if computed)
pip <- fit$pip


## list of integer indices (one vector per CS)
cs_idx <- fit$sets$cs
if (length(cs_idx) == 0) {   # no 95% credible set at this locus: write an empty result table and stop
    cat("No credible set found at chr13 52.94-53.57 Mb; writing an empty result table\n")
    fwrite(data.frame(cs_id = integer(), SNP = character(), Z = numeric(), PIP = numeric(), rsID = character(),
                      CHR = integer(), BP = integer(), A1 = character(), A2 = character(), P_gwas = numeric()),
           file = "analyses/07_fine_mapping/results/susie_results_balance_UKB_EUR_chr13_52.94M-53.57M.tsv", sep = "\t")
    quit(save = "no")
}
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
dim(cs_df)

# sanity check: Z in cs_df matches z by index
stopifnot(
    all(cs_df$Z == z[ match(cs_df$SNP, snps_ld2) ]))

# Spot-check
cs_df[1:10, c("cs_id", "SNP", "Z")]
z[ match(cs_df$SNP[1:10], snps_ld2) ]


### Adding rsIDs
gwas_res <- fread("analyses/04_gwas_sumstats/results/formatted_sumstats/UKB_EUR_balance_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_res)

# Merge
cs_df_rsid <- cs_df %>%
    left_join(
        gwas_res %>%
        select(SNP, rsID, CHR, BP, A1, A2, P) %>%
        dplyr::rename(
            # FRQ_gwas  = FRQ,
            # BETA_gwas = BETA,
            P_gwas    = P
        ),
        by = "SNP"
    )


# Save
fwrite(cs_df_rsid, file = "analyses/07_fine_mapping/results/susie_results_balance_UKB_EUR_chr13_52.94M-53.57M.tsv", sep = "\t")

