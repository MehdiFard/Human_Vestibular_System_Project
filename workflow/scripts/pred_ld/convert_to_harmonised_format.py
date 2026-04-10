import numpy as np
import pandas as pd
import dask.dataframe as dd
from scipy.stats import norm

p = snakemake.params

df = dd.read_csv(list(snakemake.input.imputed), sep="\t").compute(scheduler="threads", num_workers=snakemake.threads)

df["p_value"] = 2 * norm.sf(np.abs(df["z"]))

original = pd.read_csv(str(snakemake.input.original), sep="\t", usecols=[p.n_col], nrows=1)
df[p.n_col] = original[p.n_col].iloc[0]

out = df.rename(columns={
    "chr":     p.chr_col,
    "pos":     p.bp_col,
    "A2":      p.ref_col,
    "A1":      p.alt_col,
    "snp":     p.rsid_col,
    "beta":    p.beta_col,
    "SE":      p.se_col,
    "p_value": p.p_col,
})[[ p.rsid_col, p.chr_col, p.bp_col, p.alt_col, p.ref_col,
     p.beta_col, p.se_col, p.p_col, p.n_col,
     "imputed", "R2" ]]

out = out.sort_values([p.chr_col, p.bp_col])

out.to_csv(snakemake.output[0], sep="\t", index=False, compression="gzip")
