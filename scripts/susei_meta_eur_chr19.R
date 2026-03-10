
### Runing fine mapping using susieR
rm(list = ls())

library(susieR)
library(dplyr)
library(data.table)
library(BEDMatrix)
library(ggplot2)


bed_19 <- "Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bed"  # needs geno.bim and geno.fam in same folder

G <- BEDMatrix(bed_19)         # samples x SNPs (lazy / memory-mapped)
bim <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bim", sep = "\t", stringsAsFactors = FALSE)
fam <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.fam", sep = "\t", stringsAsFactors = FALSE)
colnames(bim) <- c("chr","id","cm","pos","a1","a2")
head(bim)

# choose locus
chr <- 19
start_bp <- 21270296 - 100000
end_bp <- 23912809 + 100000

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
ss <- read.table("gwas_sumstats_prep/output/balance_meta_eur_hg38_twas_input.sumstats.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
View(cs_df)

# sanity check: Z in cs_df matches z by index
stopifnot(
    all(cs_df$Z == z[ match(cs_df$SNP, snps_ld2) ]))

# Spot-check
cs_df[1:10, c("cs_id", "SNP", "Z")]
z[ match(cs_df$SNP[1:10], snps_ld2) ]


### Adding rsIDs
gwas_res <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta_eur_formatted_hg38_ChrPosId.tsv.gz")

# Merge
cs_df_rsid <- cs_df %>%
    left_join(
        gwas_res %>%
        select(SNP, rsID, CHR, BP, A1, A2, FRQ, BETA, P) %>%
        dplyr::rename(
            FRQ_gwas  = FRQ,
            BETA_gwas = BETA,
            P_gwas    = P
        ),
        by = "SNP"
    )
View(cs_df_rsid)


# Save
fwrite(cs_df_rsid, file = "fine_mapping/results/susie_results_balance_meta_EUR_chr19_21.27M-24.01M.tsv", sep = "\t")

