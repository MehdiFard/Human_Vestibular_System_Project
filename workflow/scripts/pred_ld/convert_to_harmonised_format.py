import numpy as np
import pandas as pd
from scipy.stats import norm

p = snakemake.params

dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input]
df = pd.concat(dfs, ignore_index=True)

df["p_value"] = 2 * norm.sf(np.abs(df["z"]))

out = df.rename(columns={
    "chr":  p.chr_col,
    "pos":  p.bp_col,
    "A2":   p.ref_col,
    "A1":   p.alt_col,
    "snp":  p.rsid_col,
    "SE":   p.se_col,
})[[ p.rsid_col, p.chr_col, p.bp_col, p.alt_col, p.ref_col,
     p.beta_col, p.se_col, p.p_col,
     "imputed", "R2" ]]

out = out.sort_values([p.chr_col, p.bp_col])

out.to_csv(snakemake.output[0], sep="\t", index=False, compression="gzip")
