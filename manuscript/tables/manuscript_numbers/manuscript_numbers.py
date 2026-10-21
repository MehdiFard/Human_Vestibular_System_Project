#!/usr/bin/env python3
"""Every number quoted in the manuscript text, recomputed from the results (107 individuals; RNA editing with the 0/0 fix
and host-gene expression regressed out). One value per row:
    section  quantity  value  previous_draft
previous_draft = the value in the draft written before the RNA editing fix (Manuscript_Balance_HGG_Final.pdf), kept to
show which numbers in the text change; empty if the draft did not quote it.
Run from the project root:  python3 manuscript/tables/manuscript_numbers/manuscript_numbers.py
"""
import re

import pandas as pd

OUT = "manuscript/tables/manuscript_numbers/manuscript_numbers.tsv"
GTF = "resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf"
P = "analyses/05_pantry/pheast"
X = "analyses/05_pantry/results/xqtl_summary/vestibular"
TW = "analyses/06_twas/results/vestibular/output/twas_hits_annot.tsv"
LP = ["lncRNA", "protein_coding"]
MODS = ["expression", "isoforms", "splicing", "alt_TSS", "alt_polyA", "stability", "RNA_editing"]
LABEL = {"expression": "Expression", "isoforms": "Isoform ratio", "splicing": "Intron excision", "alt_TSS": "Alt. TSS",
         "alt_polyA": "Alt. polyA", "stability": "RNA stability", "RNA_editing": "RNA editing"}

PREVIOUS = {   # draft values (before the RNA editing fix), by quantity
    "phenotypes, all": 486415, "phenotypes, lncRNA + protein-coding": 467739, "phenotypes, protein-coding": 389941,
    "unique genes, all": 46129, "unique genes, lncRNA + protein-coding": 34495, "phenotypes, all except RNA editing": 346794,
    "phenotypes, lncRNA + protein-coding, except RNA editing": 330101,
    "xQTL pairs, all": 3164, "unique xQTL variants, all": 2920, "unique xGenes, all": 2900, "% xGenes with >1 xQTL, all": 7.0,
    "xQTL pairs, lncRNA + protein-coding": 2627, "unique xQTL variants, lncRNA + protein-coding": 2493,
    "unique xGenes, lncRNA + protein-coding": 2403, "% xGenes with >1 xQTL, lncRNA + protein-coding": 7.4,
    "expression: xQTL pairs": 940, "expression: xQTL variants": 890, "expression: xGenes": 927,
    "splicing: xQTL pairs": 674, "splicing: xQTL variants": 626, "splicing: xGenes": 625,
    "isoforms: xQTL pairs": 362, "isoforms: xQTL variants": 360, "isoforms: xGenes": 351,
    "RNA_editing: xQTL pairs": 137, "RNA_editing: xQTL variants": 130, "RNA_editing: xGenes": 132,
    "expression / splicing ratio, xQTL pairs": 1.39, "expression / splicing ratio, xGenes": 1.48,
    "TWAS hits": 1195, "traits with TWAS hits": 19, "unique TWAS phenotypes": 802, "unique TWAS phenotypes (modality:ID)": 806,
    "unique TWAS genes, all": 282, "unique TWAS genes, lncRNA + protein-coding": 245,
    "trait-gene pairs, lncRNA + protein-coding": 394, "trait-gene pairs with non-expression hits only": 302,
    "% trait-gene pairs with non-expression hits only": 76.6, "TWAS hits with coloc PP4 > 0.8": 142,
    "% TWAS hits with coloc PP4 > 0.8": 12,
}
for tr in ["balance_meta1", "balance_meta_eur"]:
    PREVIOUS.update({f"{tr}: TWAS hits": 10, f"{tr}: RNA_editing hits": 5, f"{tr}: isoforms hits": 4,
                     f"{tr}: expression hits": 1})

# gene biotypes (Ensembl 110, as the quantification)
bt = {}
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if f[2] == "gene":
            bt.setdefault(re.search(r'gene_id "([^"]+)"', f[8]).group(1), re.search(r'gene_biotype "([^"]+)"', f[8]).group(1))
bt = pd.Series(bt)

rows = []
add = lambda section, quantity, value: rows.append((section, quantity, value))

# ---- phenotypes --------------------------------------------------------------------------------------------------------
ph = pd.concat([pd.read_csv(f"{P}/input/vestibular/phenotypes/{m}.bed.gz", sep="\t", usecols=[3]).iloc[:, 0].to_frame("id")
                .assign(m=m) for m in MODS])
ph["gene"] = ph.id.str.split("__").str[0]
ph["bt"] = ph.gene.map(bt)
S = "Phenotypes"
add(S, "phenotypes, all", len(ph))
add(S, "phenotypes, lncRNA + protein-coding", int(ph.bt.isin(LP).sum()))
add(S, "phenotypes, protein-coding", int((ph.bt == "protein_coding").sum()))
add(S, "unique genes, all", ph.gene.nunique())
add(S, "unique genes, lncRNA + protein-coding", ph[ph.bt.isin(LP)].gene.nunique())
add(S, "phenotypes, all except RNA editing", int((ph.m != "RNA_editing").sum()))
add(S, "phenotypes, lncRNA + protein-coding, except RNA editing", int((ph.bt.isin(LP) & (ph.m != "RNA_editing")).sum()))

# ---- cross-modality xQTLs ------------------------------------------------------------------------------------------------
x = pd.read_csv(f"{P}/output/vestibular/qtl/cross_modality.cis_independent_qtl.txt.gz", sep="\t")
x["m"] = x.phenotype_id.str.split(":").str[0]
S = "xQTLs (cross-modality mapping)"
for lab, d in [("all", x), ("lncRNA + protein-coding", x[x.group_id.map(bt).isin(LP)])]:
    add(S, f"xQTL pairs, {lab}", len(d))
    add(S, f"unique xQTL variants, {lab}", d.variant_id.nunique())
    add(S, f"unique xGenes, {lab}", d.group_id.nunique())
    add(S, f"% xGenes with >1 xQTL, {lab}", round(100 * (d.groupby("group_id").size() > 1).mean(), 1))
d = x[x.group_id.map(bt).isin(LP)]
for m in ["expression", "splicing", "isoforms", "RNA_editing"]:
    dm = d[d.m == m]
    add(S, f"{m}: xQTL pairs", len(dm)); add(S, f"{m}: xQTL variants", dm.variant_id.nunique())
    add(S, f"{m}: xGenes", dm.group_id.nunique())
e, s = d[d.m == "expression"], d[d.m == "splicing"]
add(S, "expression / splicing ratio, xQTL pairs", round(len(e) / len(s), 2))
add(S, "expression / splicing ratio, xGenes", round(e.group_id.nunique() / s.group_id.nunique(), 2))

# ---- TWAS ------------------------------------------------------------------------------------------------------------
t = pd.read_csv(TW, sep="\t")
S = "TWAS"
add(S, "TWAS hits", len(t)); add(S, "traits with TWAS hits", t.TRAIT.nunique())
add(S, "unique TWAS phenotypes", t.ID.nunique()); add(S, "unique TWAS phenotypes (modality:ID)", (t.MODALITY + ":" + t.ID).nunique())
add(S, "unique TWAS genes, all", t.ENSG.nunique()); add(S, "unique TWAS genes, lncRNA + protein-coding", t[t.gene_type.isin(LP)].ENSG.nunique())
tl = t[t.gene_type.isin(LP)]
g = tl.groupby(["TRAIT", "ENSG"]).MODALITY.agg(lambda s: ("expression" in set(s), bool(set(s) - {"expression"})))
ne = sum(1 for a, b in g if b and not a)
add(S, "trait-gene pairs, lncRNA + protein-coding", len(g))
add(S, "trait-gene pairs with non-expression hits only", ne)
add(S, "% trait-gene pairs with non-expression hits only", round(100 * ne / len(g), 1))
pp = int((t["COLOC.PP4"] > 0.8).sum())
add(S, "TWAS hits with coloc PP4 > 0.8", pp); add(S, "% TWAS hits with coloc PP4 > 0.8", round(100 * pp / len(t)))
for tr in ["balance_meta1", "balance_meta_eur"]:
    S = f"Dizziness TWAS ({tr})"
    b = t[t.TRAIT == tr]
    add(S, f"{tr}: TWAS hits", len(b))
    for m, n in b.MODALITY.value_counts().items():
        add(S, f"{tr}: {m} hits", int(n))
    for _, r in b[b["COLOC.PP4"] > 0.8].iterrows():
        add(S, f"{tr}: coloc PP4 > 0.8, {r.MODALITY} {r.gene_name} {r.ID.split('__')[-1]}", r["COLOC.PP4"])

# ---- figure numbers (xQTL summary tables of 05_pantry; lncRNA + protein-coding) ---------------------------------------
q = pd.read_csv(f"{X}.comb.qtls.tsv", sep="\t")
pos = pd.read_csv(f"{X}.comb.qtls.rel_pos.tsv", sep="\t")
S = "Figures"
inside = int(pos.rel_pos_gene.between(-1, 2).sum())
add(S, "Supp Fig 3d: xQTLs within one gene length of the gene start/end", inside)
add(S, "Supp Fig 3d: % of xQTLs within one gene length of the gene start/end", round(100 * inside / len(pos), 1))
n_per_gene = q.groupby("gene_id").size()
add(S, "Supp Fig 3b: xGenes with >1 xQTL", int((n_per_gene > 1).sum()))
add(S, "Supp Fig 3b: % xGenes with >1 xQTL", round(100 * (n_per_gene > 1).mean(), 1))
rank = q["rank"].where(q["rank"] < 10, 10)
for k, n in rank.value_counts().sort_index().items():
    add(S, f"Supp Fig 3c: xQTLs of rank {'10+' if k == 10 else k}", int(n))
for m, n in q.groupby("modality").gene_id.nunique().items():
    add(S, f"Fig 4b: xGenes, {LABEL[m]}", int(n))

fmt = lambda v: "" if v is None or v != v else (str(int(v)) if float(v).is_integer() and not isinstance(v, float) or
                                                 (isinstance(v, float) and v.is_integer() and abs(v) >= 2) else str(v))
out = pd.DataFrame(rows, columns=["section", "quantity", "value"])
prev = out.quantity.map(PREVIOUS)
changed = out[prev.notna() & (out.value.astype(float) != prev.astype(float))]
out["value"] = [fmt(v) for v in out.value]
out["previous_draft"] = [fmt(v) for v in prev]
out.to_csv(OUT, sep="\t", index=False)
print(f"written {OUT}: {len(out)} numbers, {len(changed)} differ from the draft")


### END