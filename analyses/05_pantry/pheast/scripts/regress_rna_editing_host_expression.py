#!/usr/bin/env python3
# Regress host-gene expression out of RNA editing phenotypes before QTL mapping (approach of Li et al. 2022):
#   for each editing site, regress its normalized editing level on the expression (log2 TPM + 1) of the gene(s) it
#   lies in, keep the residuals and re-normalize them (rank-based inverse normal, = 'quantile_int').
#   Sites in several genes use all of them; sites without an expressed host gene are left unchanged.
# Pantry specifics:
#   - input: fixed RNA editing phenotypes for the paper's 107 individuals (already quantile + INT normalized)
#   - host genes of a site: all genes it is assigned to in Pantry (phenotype IDs gene__site; one row per gene,
#     every copy gets the same residual)
#   - expression: analyses/05_pantry/phenotyping/output/vestibular/unnorm/expression.bed (kallisto gene TPM, sample IDs)
# Input:  input/vestibular/phenotypes_unregressed/RNA_editing.bed.gz (0/0-fixed, individual IDs, 107 individuals)
# Output: input/vestibular/phenotypes/RNA_editing.{bed.gz,phenotype_groups.txt}
# Run from analyses/05_pantry/pheast (part of prepare_inputs.sh):  python3 scripts/regress_rna_editing_host_expression.py
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm

SRC = "input/vestibular/phenotypes_unregressed"
TPM = "../phenotyping/output/vestibular/unnorm/expression.bed"
MAP = "input/vestibular/vestibular.sample_map.txt"
OUT = "input/vestibular/phenotypes"
os.makedirs(OUT, exist_ok=True)

bed = pd.read_csv(f"{SRC}/RNA_editing.bed.gz", sep="\t")
subjects = list(bed.columns[4:])
smap = pd.read_csv(MAP, sep=r"\s+", header=None, dtype=str)
s2i = dict(zip(smap[0], smap[1]))

tpm = pd.read_csv(TPM, sep="\t", index_col=3).iloc[:, 3:]
tpm.columns = [s2i.get(str(c), None) for c in tpm.columns]
tpm = tpm.loc[:, [c for c in tpm.columns if c in subjects]]
missing = [s for s in subjects if s not in tpm.columns]
assert not missing, f"no expression for {missing[:5]}"
expr = np.log2(tpm[subjects] + 1)
expr = expr[expr.std(axis=1) > 0]
print(f"{len(subjects)} individuals, expression for {len(expr):,} genes")

ids = bed["phenotype_id"]
site = ids.str.split("__", n=1).str[1]
gene = ids.str.split("__", n=1).str[0]
site_genes = pd.DataFrame({"site": site, "gene": gene}).query("gene in @expr.index").groupby("site").gene.apply(lambda g: sorted(set(g)))

values = bed[subjects].to_numpy(float)
E = expr.to_numpy()
row = {g: i for i, g in enumerate(expr.index)}
n_genes, r2, done = [], [], {}
for i, s in enumerate(site):
    g = site_genes.get(s, [])
    n_genes.append(len(g))
    if not g:
        continue
    if s not in done:                                   # same site in several genes: compute once
        y = values[i]
        X = np.column_stack([np.ones(len(y))] + [E[row[x]] for x in g])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        done[s] = (norm.ppf(rankdata(resid) / (len(resid) + 1)), 1 - resid.var() / y.var())
    values[i] = done[s][0]
    r2.append(done[s][1])

bed[subjects] = values
tmp = f"{OUT}/RNA_editing.bed"
bed.to_csv(tmp, sep="\t", index=False, float_format="%g")
subprocess.run(["bgzip", "-f", tmp], check=True)
subprocess.run(["tabix", "-f", "-p", "bed", tmp + ".gz"], check=True)
shutil.copy2(f"{SRC}/RNA_editing.phenotype_groups.txt", f"{OUT}/RNA_editing.phenotype_groups.txt")

n_genes = np.array(n_genes); r2 = np.array(r2)
print(f"phenotypes: {len(n_genes):,} | 1 host gene: {(n_genes == 1).sum():,} | 2+: {(n_genes >= 2).sum():,} | "
      f"no expressed host gene (unchanged): {(n_genes == 0).sum():,}")
print(f"variance of editing explained by host-gene expression: median {np.median(r2):.3f}, "
      f">10% for {np.mean(r2 > 0.1):.1%}, >25% for {np.mean(r2 > 0.25):.1%} of phenotypes")
