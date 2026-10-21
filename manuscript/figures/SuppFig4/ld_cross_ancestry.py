#!/usr/bin/env python3
# Pairwise LD (r2) between the three ZNF91-locus lead variants in each 1000 Genomes superpopulation:
#   rs695017 (GWAS lead, 19_23421111), rs167154 (isoform QTL, 19_23307401), rs607168 (RNA editing QTL, 19_23459336).
# Reference: resources/ldref_1kg_hg38 (1000 Genomes, GRCh38); superpopulations from the IGSR sample table.
# Public reference data only. Output: data/ld_cross_ancestry_3snps.tsv (read by suppfig_coloc_ld.R)
# Run from project root (compute node): python3 manuscript/figures/SuppFig4/ld_cross_ancestry.py
import os
import yaml
import itertools
import numpy as np
import pandas as pd

G = "resources/ldref_1kg_hg38/1kg_ldref_hg38_chr19"
cfg = yaml.safe_load(open("config/config.yaml"))   # + lab-specific values (config_lab.yaml, not in git)
if os.path.exists("config/config_lab.yaml"):
    lab = yaml.safe_load(open("config/config_lab.yaml"))
    for k, v in lab.items():
        cfg[k] = {**cfg.get(k, {}), **v} if isinstance(v, dict) else v
IGSR = cfg["refs"]["igsr_samples"]   # 1000 Genomes sample -> population
OUT = "manuscript/figures/SuppFig4/data/ld_cross_ancestry_3snps.tsv"
SNPS = {"rs695017": 23421111, "rs167154": 23307401, "rs607168": 23459336}

fam = pd.read_csv(G + ".fam", sep=r"\s+", header=None)[1]
pop = fam.map(pd.read_csv(IGSR, sep="\t").set_index("Sample name")["Superpopulation code"]).to_numpy()
bim = pd.read_csv(G + ".bim", sep="\t", header=None, names=["chr", "snp", "cm", "pos", "a1", "a2"], dtype={"snp": str})
n = len(fam)
bed = np.memmap(G + ".bed", dtype=np.uint8, mode="r", offset=3, shape=(len(bim), (n + 3) // 4))
CODE = np.array([2.0, np.nan, 1.0, 0.0])
geno = {}
for rs, p in SNPS.items():
    i = np.flatnonzero(bim.pos.to_numpy() == p)
    assert len(i) == 1, f"{rs}: {len(i)} rows at {p}"
    geno[rs] = CODE[((bed[i[0]][:, None] >> np.array([0, 2, 4, 6])) & 3).ravel()[:n]]

rows = []
for sp in ["EUR", "AMR", "EAS", "SAS", "AFR"]:
    m = pop == sp
    for a, b in itertools.combinations(SNPS, 2):
        x, y = geno[a][m], geno[b][m]
        ok = ~np.isnan(x) & ~np.isnan(y)
        rows.append(dict(pop=sp, n=int(ok.sum()), snp1=a, snp2=b, r2=round(float(np.corrcoef(x[ok], y[ok])[0, 1] ** 2), 3)))
d = pd.DataFrame(rows)
d.to_csv(OUT, sep="\t", index=False)
print(d.pivot(index="pop", columns=["snp1", "snp2"], values="r2").loc[["EUR", "AMR", "EAS", "SAS", "AFR"]])
print(d.groupby("pop").n.first())


### END