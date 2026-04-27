
### Runing fine mapping using susieR
rm(list = ls())

library(susieR)
library(dplyr)
library(data.table)
library(BEDMatrix)
library(ggplot2)


## Read nominal isoform qtl mapping results
library(arrow)
sqtl <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/isoforms/vestibular_isoforms.cis_qtl_pairs.19.parquet")
dim(sqtl)
View(sqtl[1:20, ])

# I need to subset to one isoform,
target_isoform <- "ENSG00000167232__ENST00000300619"
length(which(sqtl$phenotype_id == target_isoform))

# Extract relevant columns and compute Z-scores
sqtl_cln <- sqtl %>%
    filter(phenotype_id == target_isoform) %>%
    select(phenotype_id, SNP = variant_id, pval_nominal, slope, slope_se) %>%
    mutate(
        Z = slope / slope_se,
        pos = as.numeric(sub("^[^_]+_([^_]+)_.*$", "\\1", SNP)),
    )

dim(sqtl_cln)
length(unique(sqtl_cln$SNP))
head(sqtl_cln)
range(sqtl_cln$Z)
range(sqtl_cln$pos)
summary(sqtl_cln$pos)



bed_19 <- "Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bed"  # needs geno.bim and geno.fam in same folder

G <- BEDMatrix(bed_19)         # samples x SNPs (lazy / memory-mapped)
bim <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.bim", sep = "\t", stringsAsFactors = FALSE)
fam <- read.table("Pantry/pheast_no_outlier/input/LDREF/1kg_ldref_hg38_chr19.fam", sep = "\t", stringsAsFactors = FALSE)
colnames(bim) <- c("chr","id","cm","pos","a1","a2")
head(bim)

# choose locus
chr <- 19
start_bp <- 22430078 - 100000
end_bp <- 24390886 + 100000

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
library(matrixStats)
sd <- matrixStats::colSds(X)
X <- X[, sd > 0, drop = FALSE]
R <- crossprod(scale(X)) / (nrow(X) - 1)
diag(R) <- 1
R[1:5, 1:5]

# Fix SNP IDs in R to match bim
new_names <- sub("_[^_]+$", "", colnames(R))
colnames(R) <- new_names
rownames(R) <- sub("_[^_]+$", "", rownames(R))
R[1:5, 1:5]

###--------------------- SusieR
## Get the SNP IDs in the LD matrix order
snps_ld <- bim$id[idx]   # same order as columns in X and in R
length(snps_ld)

## Match sumstats to LD SNP order
m <- match(snps_ld, sqtl_cln$SNP)
keep <- which(!is.na(m))
length(keep)


## subset to SNPs present in both
snps_ld2 <- snps_ld[keep]
R2 <- R[keep, keep, drop = FALSE]
sqtl_cln2 <- sqtl_cln[m[keep], ]
z <- sqtl_cln2$Z
# Check alignment
table(colnames(R2) == sqtl_cln2$SNP)
table(rownames(R2) == sqtl_cln2$SNP)

## Basic sanity checks
stopifnot(length(z) == ncol(R2), ncol(R2) == nrow(R2))

# (optional) check no NA
stopifnot(!any(is.na(z)), !any(is.na(R2)))

## Choose effective N
# Use the actual isoform-QTL sample size.
n_eff <- 107L

## Stabilize LD if needed (ridge)
R2 <- (R2 + t(R2)) / 2
diag(R2) <- 1
# R2 <- R2 + diag(1e-5, nrow(R2))  # helps if near-singular


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


## Diagnostic output
cat("\n========== SuSiE Fine-Mapping Results ==========\n")
cat("Matched SNPs:", length(snps_ld2), "\n")
cat("SuSiE converged:", fit$converged, "\n")
cat("Max PIP:", round(max(pip, na.rm = TRUE), 4), "\n")
cat("Number of credible sets:", length(fit$sets$cs), "\n\n")

## Create all-SNPs PIP table (for all matched SNPs)
pip_df <- data.frame(
    SNP   = snps_ld2,
    Z     = z,
    PIP   = pip,
    stringsAsFactors = FALSE
)
pip_df <- pip_df[order(-pip_df$PIP), ]  # sort by PIP descending
head(pip_df, 30)

## list of integer indices (one vector per CS)
cs_idx <- fit$sets$cs
cs_snps <- lapply(cs_idx, function(ii) snps_ld2[ii]) # SNP IDs per credible set

## Make a result table for credible sets
if (length(cs_idx) > 0) {
    cat("Creating credible set results table...\n")
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
    cat("Credible sets table:\n")
    print(head(cs_df, 20))
    
    # sanity check: Z in cs_df matches z by index
    stopifnot(
        all(cs_df$Z == z[ match(cs_df$SNP, snps_ld2) ]))

} else {
    cat("WARNING: No credible sets found. PIPs may be too small or signal spread across many SNPs.\n")
    cat("Proceeding to save all-SNP PIP table only.\n\n")
    cs_df <- NULL
}


### Adding rsIDs
gwas_res <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")

## Save all SNPs with PIPs
pip_df_rsid <- pip_df %>%
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

fwrite(pip_df_rsid, file = "fine_mapping/results/susie_results_isoformQTL_nominal_chr19_21.42M-23.7M_all_snps.tsv", sep = "\t")
cat("Saved all-SNP PIP table to: fine_mapping/results/susie_results_isoformQTL_nominal_chr19_21.42M-23.7M_all_snps.tsv\n\n")

## Save credible sets (only if found)
if (!is.null(cs_df)) {
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
    
    fwrite(cs_df_rsid, file = "fine_mapping/results/susie_results_isoformQTL_nominal_chr19_21.42M-23.7M.tsv", sep = "\t")
    cat("Saved credible sets to: fine_mapping/results/susie_results_isoformQTL_nominal_chr19_21.42M-23.7M.tsv\n")
} else {
    cat("(Credible sets file not generated—no CS found)\n")
}

cat("\n========== Analysis Complete ==========\n")


