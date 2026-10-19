import pandas as pd

p = snakemake.params

# Read tau coefficients — one row per annotation in regression order
results = pd.read_csv(snakemake.input.results, sep="\t")
tau = results["Coefficient"].values  # shape: (n_annot,)
n_annot = len(tau)

# Read median N from the munged sumstats file.
# The LDSC model is E[chi2_j] = N * Σ_c(tau_c * l2_cj) + intercept,
# so without N the Coefficient values (~1e-8) produce near-zero signal.
sumstats = pd.read_csv(snakemake.input.sumstats, sep="\t", usecols=["SNP", "A1", "A2"])

# Read and concatenate per-chromosome annotation LD scores
# File format: CHR  SNP  BP  <annot1_L2>  <annot2_L2>  ...
ldscore_dfs = [
    pd.read_csv(f"{p.annot_prefix}{chrom}.l2.ldscore.gz", sep="\t", compression="gzip")
    for chrom in p.chroms
]
ld = pd.concat(ldscore_dfs, ignore_index=True)

# Annotation columns are the n_annot columns immediately following CHR, SNP, BP
annot_cols = ld.columns[3 : 3 + n_annot].tolist()

assert len(annot_cols) == n_annot, (
    f"Annotation count mismatch: {n_annot} tau coefficients but "
    f"{len(annot_cols)} annotation columns in LD score files"
)

# E[chi2_j] = N * Σ_c(tau_c * l2_cj) + 1
# Intercept assumed ≈ 1 (well-controlled GWAS). For studies with substantial
# population stratification, parse the intercept from the .log file instead.
K = snakemake.params.case_proportion
N = 4 * K * (1 - K) * snakemake.params.samples
ld["chi2_expected"] = 1.0 + N * (ld[annot_cols].values @ tau)
# Duplicate for the sake of merging with sffdr

ld = ld.merge(sumstats, on="SNP")

ld[["SNP", "CHR", "BP", "A1", "A2", "chi2_expected"]].to_csv(
    snakemake.output[0], sep="\t", index=False, compression="gzip"
)
