#!/usr/bin/env python3
# r2 of every 1000 Genomes SNP in chr19:22.95-23.85 Mb with rs695017 (19_23421111_A_G), EUR samples only.
# Replaces LDlink colouring for Fig 5b (LDlink returned no r2 for rs167154 / rs607168 because of rsID mismatches).
# Reference: resources/ldref_1kg_hg38 (1000 Genomes, GRCh38, the TWAS LD reference);
# EUR = superpopulation in the IGSR sample table. Public reference data only.
# Run from project root (compute node): python3 manuscript/figures/Fig5/ld_rs695017_1kg_eur.py
import os
import yaml
import numpy as np
import pandas as pd

G = "resources/ldref_1kg_hg38/1kg_ldref_hg38_chr19"
cfg = yaml.safe_load(open("config/config.yaml"))   # + lab-specific values (config_lab.yaml, not in git)
if os.path.exists("config/config_lab.yaml"):
    lab = yaml.safe_load(open("config/config_lab.yaml"))
    for k, v in lab.items():
        cfg[k] = {**cfg.get(k, {}), **v} if isinstance(v, dict) else v
IGSR = cfg["refs"]["igsr_samples"]   # 1000 Genomes sample -> population
OUT = "manuscript/figures/Fig5/data/ld_rs695017_1kgEUR.tsv"
INDEX_POS, REGION = 23421111, (22_950_000, 23_850_000)

fam = pd.read_csv(G + ".fam", sep=r"\s+", header=None)[1]
pop = pd.read_csv(IGSR, sep="\t").set_index("Sample name")["Superpopulation code"]
eur = (fam.map(pop) == "EUR").to_numpy()
pd.Series(fam[eur]).to_csv("manuscript/figures/Fig5/data/1kg_EUR_samples.txt", index=False, header=False)

bim = pd.read_csv(G + ".bim", sep="\t", header=None, names=["chr", "snp", "cm", "pos", "a1", "a2"], dtype={"snp": str})
n = len(fam)
bed = np.memmap(G + ".bed", dtype=np.uint8, mode="r", offset=3, shape=(len(bim), (n + 3) // 4))
CODE = np.array([2.0, np.nan, 1.0, 0.0])


def geno(i):
    return CODE[((bed[i][:, None] >> np.array([0, 2, 4, 6])) & 3).ravel()[:n]][eur]


idx = np.flatnonzero(bim.pos.between(*REGION).to_numpy())
ref_i = int(np.flatnonzero(bim.pos.to_numpy() == INDEX_POS)[0])
x = geno(ref_i)
r2 = []
for i in idx:
    y = geno(i)
    m = ~np.isnan(x) & ~np.isnan(y)
    r2.append(np.corrcoef(x[m], y[m])[0, 1] ** 2 if y[m].std() > 0 else np.nan)
out = bim.iloc[idx][["chr", "pos", "snp", "a1", "a2"]].assign(r2=np.round(r2, 4))
out.to_csv(OUT, sep="\t", index=False)
print(f"EUR samples: {eur.sum()} | SNPs: {len(out)} | index {bim.snp[ref_i]}")
print(out[out.pos.isin([23421111, 23307401, 23459336])].to_string(index=False))
