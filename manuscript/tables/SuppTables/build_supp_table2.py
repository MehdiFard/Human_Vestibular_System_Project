#!/usr/bin/env python3
# Build the updated Supplementary Table 2 (QTL + TWAS) with the fixed RNA editing results.
#   Input : Supplementary_Table2_QTL_TWAS.xlsx (paper version)
#   Output: manuscript/tables/SuppTables/Supplementary_Table2_QTL_TWAS.xlsx
# Updated sheets (recomputed from files, same definitions as the paper version, which they reproduce exactly
# when run on the paper's files):
#   S1  phenotypes per modality (lncRNA + protein-coding genes)       RNA editing from the fixed run
#   S2  cross-modality independent QTLs, sorted by pval_beta          fixed cross-modality run
#   S3  significant pairs / genes / variants (lncRNA + protein-coding) separate: fixed RNA editing; cross: fixed run
#   S6  TWAS hits (enet/top1 sign-corrected), sorted by COLOC.PP4      RNA editing hits from the fixed run
#   S7  TWAS trait x gene: expression vs other modalities             from the new S6 (lncRNA + protein-coding)
# Unchanged data: S4, S5, S8, S9. All sheets restyled (Arial, header style, frozen header, filters).
# Trait name balance_meta1 -> balance_meta in S6/S7 to match S5.
# Run from project root with a Python that has openpyxl + pandas.
import os
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = os.getcwd()
SRC = "Supplementary_Table2_QTL_TWAS.xlsx"
OUT = "manuscript/tables/SuppTables/Supplementary_Table2_QTL_TWAS.xlsx"
PAPER = "analyses/05_pantry/pheast"
## RNA editing (and cross-modality, TWAS) with host-gene expression regressed out (2026-10-14); same layout as before
FIX = "analyses/05_pantry/pheast"
GTF = "resources/annotation/ensembl110/Homo_sapiens.GRCh38.110.chr.gtf"
MODS = ["expression", "splicing", "isoforms", "alt_polyA", "stability", "alt_TSS", "RNA_editing"]
S1_LABEL = {"expression": "Expression", "isoforms": "Isoform ratio", "splicing": "Intron excision",
            "alt_TSS": "Alt. TSS", "alt_polyA": "Alt. polyA", "stability": "RNA stability",
            "RNA_editing": "RNA editing"}
S1_ORDER = ["expression", "isoforms", "splicing", "alt_TSS", "alt_polyA", "stability", "RNA_editing"]
TRAIT_RENAME = {"balance_meta1": "balance_meta"}
KEEP = ["lncRNA", "protein_coding"]


def run_of(mod):
    """Fixed run for RNA editing (and cross_modality), paper run for everything else."""
    return FIX if mod in ("RNA_editing", "cross_modality") else PAPER


# ---------------------------------------------------------------- gene biotypes (Ensembl 110)
def load_biotypes():
    ids, types = [], []
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] != "gene":
                continue
            a = p[8]
            ids.append(a.split('gene_id "', 1)[1].split('"', 1)[0])
            types.append(a.split('gene_biotype "', 1)[1].split('"', 1)[0] if 'gene_biotype "' in a else None)
    return pd.Series(types, index=ids)


bt = load_biotypes()
keep = lambda genes: genes.map(bt).isin(KEEP)

# ---------------------------------------------------------------- S1: phenotypes per modality and biotype
s1 = []
for m in S1_ORDER:
    ids = pd.read_csv(f"{run_of(m)}/input/vestibular/phenotypes/{m}.bed.gz", sep="\t", usecols=[3]).iloc[:, 0]
    g = ids.str.split("__").str[0].map(bt)
    for b in KEEP:
        s1.append((S1_LABEL[m], b, int((g == b).sum())))

# ---------------------------------------------------------------- S3: significant pairs / genes / variants
def indep(mod):
    d = pd.read_csv(f"{run_of(mod)}/output/vestibular/qtl/{mod}.cis_independent_qtl.txt.gz", sep="\t")
    if "group_id" not in d:
        d["group_id"] = d.phenotype_id
    return d[keep(d.group_id)]


cross = indep("cross_modality")
cross["mod"] = cross.phenotype_id.str.split(":").str[0]
s3 = {}
for m in MODS:
    sep, cr = indep(m), cross[cross["mod"] == m]
    s3[m] = dict(pairs=(len(sep), len(cr)), genes=(sep.group_id.nunique(), cr.group_id.nunique()),
                 variants=(sep.variant_id.nunique(), cr.variant_id.nunique()))

# ---------------------------------------------------------------- S2: cross-modality independent QTLs
s2 = pd.read_csv(f"{FIX}/output/vestibular/qtl/cross_modality.cis_independent_qtl.txt.gz", sep="\t")
s2 = s2.sort_values("pval_beta", kind="stable")

# ---------------------------------------------------------------- S6 / S7: TWAS
tw = pd.read_csv(f"{FIX}/twas/vestibular/output/twas_hits_annot.tsv", sep="\t")  # sign-corrected, see README
tw["TRAIT"] = tw.TRAIT.replace(TRAIT_RENAME)
s6 = tw.rename(columns={"TRAIT": "GWAS.TRAIT", "MODALITY": "RNA.MODALITY", "gene_name": "GENE.NAME",
                        "gene_type": "GENE.TYPE", "EQTL.ID": "xQTL.ID"})
s6 = s6[["GWAS.TRAIT", "RNA.MODALITY", "ID", "ENSG", "GENE.NAME", "GENE.TYPE", "CHR", "BEST.GWAS.ID",
         "BEST.GWAS.Z", "xQTL.ID", "TWAS.Z", "TWAS.P", "COLOC.PP4"]]
s6 = s6.sort_values(["COLOC.PP4", "TWAS.P"], ascending=[False, True], kind="stable")

t = tw[tw.gene_type.isin(KEEP)]
g = t.groupby(["TRAIT", "ENSG"]).MODALITY.agg(lambda s: ("expression" in set(s), bool(set(s) - {"expression"})))
s7 = g.map(lambda x: "Expr. & other(s)" if x[0] and x[1] else ("Expression only" if x[0] else "Non-expr. only"))
s7 = s7.reset_index().rename(columns={"TRAIT": "trait", "ENSG": "gene_id", "MODALITY": "modality_hits"})
s7 = s7.sort_values(["trait", "gene_id"]).reset_index(drop=True)

# ================================================================ workbook
wb = openpyxl.load_workbook(SRC)
F = "Arial"
BODY, BOLD = Font(name=F, size=10), Font(name=F, size=10, bold=True)
TITLE = Font(name=F, size=12, bold=True)
NOTE = Font(name=F, size=9, italic=True, color="595959")
HFILL = PatternFill("solid", fgColor="D9E1F2")
THIN = Border(bottom=Side(style="thin", color="808080"))
SCI = "0.00E+00"


def replace_sheet(name):
    i = wb.sheetnames.index(name)
    wb.remove(wb[name])
    return wb.create_sheet(name, i)


def header(ws, row, ncol, col0=1):
    for c in range(col0, col0 + ncol):
        cell = ws.cell(row=row, column=c)
        cell.font, cell.fill, cell.border = BOLD, HFILL, THIN
        cell.alignment = Alignment(vertical="center", wrap_text=False)


def write_df(ws, df, row0=1, sci_cols=()):
    for j, col in enumerate(df.columns, 1):
        ws.cell(row=row0, column=j, value=col)
    header(ws, row0, df.shape[1])
    vals = df.astype(object).where(df.notna(), None).values.tolist()
    for i, r in enumerate(vals, row0 + 1):
        for j, v in enumerate(r, 1):
            if isinstance(v, (np.integer,)):
                v = int(v)
            elif isinstance(v, (np.floating,)):
                v = float(v)
            ws.cell(row=i, column=j, value=v)
    for j, col in enumerate(df.columns, 1):
        if col in sci_cols:
            for i in range(row0 + 1, row0 + 1 + len(df)):
                ws.cell(row=i, column=j).number_format = SCI
    ws.freeze_panes = ws.cell(row=row0 + 1, column=1)
    ws.auto_filter.ref = f"A{row0}:{get_column_letter(df.shape[1])}{row0 + len(df)}"


# ---- S1
ws = replace_sheet("S1_N. of phenotypes")
ws["A1"] = "Number of phenotypes per modality (lncRNA and protein-coding genes)"
ws["A1"].font = TITLE
ws["A2"] = "Vestibular otolith organs, n = 107. Biotypes from Ensembl release 110."
ws["A2"].font = NOTE
ws["A4"], ws["B4"] = "Modality", "N. phenotypes"
header(ws, 4, 2)
n_mod, bt_row0 = len(S1_ORDER), 4 + len(S1_ORDER) + 5          # biotype table header row
bt_first, bt_last = bt_row0 + 1, bt_row0 + len(s1)
for k, m in enumerate(S1_ORDER):
    r = 5 + k
    ws.cell(row=r, column=1, value=S1_LABEL[m])
    ws.cell(row=r, column=2, value=f"=SUMIFS($C${bt_first}:$C${bt_last},$A${bt_first}:$A${bt_last},A{r})")
tot = 5 + n_mod
ws.cell(row=tot, column=1, value="Total").font = BOLD
ws.cell(row=tot, column=2, value=f"=SUM(B5:B{tot - 1})").font = BOLD
ws.cell(row=bt_row0 - 2, column=1, value="Number of phenotypes per modality by gene biotype").font = TITLE
for j, h in enumerate(["Modality", "Gene biotype", "N. phenotypes"], 1):
    ws.cell(row=bt_row0, column=j, value=h)
header(ws, bt_row0, 3)
for i, (m, b, n) in enumerate(s1, bt_first):
    ws.cell(row=i, column=1, value=m); ws.cell(row=i, column=2, value=b); ws.cell(row=i, column=3, value=n)
for r in list(range(5, tot + 1)) + list(range(bt_first, bt_last + 1)):
    for c in (2, 3):
        ws.cell(row=r, column=c).number_format = "#,##0"
ws.column_dimensions["A"].width, ws.column_dimensions["B"].width, ws.column_dimensions["C"].width = 18, 16, 15

# ---- S2
ws = replace_sheet("S2_xQTLs_cross_modality")
write_df(ws, s2, sci_cols=("pval_true_df", "pval_nominal", "pval_perm", "pval_beta"))

# ---- S3
ws = replace_sheet("S3_N._of_sig_xqtls_xgenes")
blocks = [("pairs", "Number of significant independent variant–phenotype pairs per modality (lncRNA and protein-coding genes)"),
          ("genes", "Number of unique genes with a significant QTL per modality (lncRNA and protein-coding genes)"),
          ("variants", "Number of unique significant QTL variants per modality (lncRNA and protein-coding genes)")]
r = 1
for key, title in blocks:
    ws.cell(row=r, column=1, value=title).font = TITLE
    for j, h in enumerate(["Modality", "Separate mapping", "Cross-modality mapping"], 1):
        ws.cell(row=r + 1, column=j, value=h)
    header(ws, r + 1, 3)
    first = r + 2
    for k, m in enumerate(MODS):
        ws.cell(row=first + k, column=1, value=m)
        ws.cell(row=first + k, column=2, value=s3[m][key][0])
        ws.cell(row=first + k, column=3, value=s3[m][key][1])
    last = first + len(MODS) - 1
    ws.cell(row=last + 1, column=1, value="Sum").font = BOLD
    for c, L in ((2, "B"), (3, "C")):
        ws.cell(row=last + 1, column=c, value=f"=SUM({L}{first}:{L}{last})").font = BOLD
    for rr in range(first, last + 2):
        for c in (2, 3):
            ws.cell(row=rr, column=c).number_format = "#,##0"
    r = last + 4
ws.cell(row=r - 1, column=1, value="Significant = independent cis-QTLs from tensorQTL (gene-level q < 0.05, forward-backward conditional mapping).").font = NOTE
ws.column_dimensions["A"].width, ws.column_dimensions["B"].width, ws.column_dimensions["C"].width = 16, 18, 24

# ---- S6
ws = replace_sheet("S6_TWAS_Vestibular")
write_df(ws, s6, sci_cols=("TWAS.P",))

# ---- TWAS dictionary (matches S6 columns)
ws = replace_sheet("TWAS_dictionary_columns")
dic = [("GWAS.TRAIT", "GWAS trait (see S5)"),
       ("RNA.MODALITY", "RNA phenotype modality"),
       ("ID", "RNA phenotype ID (from the TWAS weights file)"),
       ("ENSG", "Ensembl gene ID"),
       ("GENE.NAME", "Gene symbol (Ensembl release 110)"),
       ("GENE.TYPE", "Gene biotype (Ensembl release 110)"),
       ("CHR", "Chromosome"),
       ("BEST.GWAS.ID", "Variant ID (chr_pos_ref_alt, GRCh38) of the most significant GWAS variant in the locus"),
       ("BEST.GWAS.Z", "GWAS Z-score of BEST.GWAS.ID"),
       ("xQTL.ID", "Variant ID (chr_pos_ref_alt, GRCh38) of the best QTL for this RNA phenotype in the weights"),
       ("TWAS.Z", "TWAS Z-score (primary statistic): sign = direction of the association between the genetically "
                  "predicted RNA phenotype and the trait"),
       ("TWAS.P", "TWAS P-value; hits have TWAS.P < 5e-8 / 7 (7 modalities)"),
       ("COLOC.PP4", "Colocalization posterior probability of H4: the RNA phenotype and the trait share one causal variant")]
write_df(ws, pd.DataFrame(dic, columns=["Column", "Description"]))
ws.auto_filter.ref = None
ws.column_dimensions["A"].width, ws.column_dimensions["B"].width = 16, 110

# ---- S7
ws = replace_sheet("S7_TWAS hits_expr_non-expr")
write_df(ws, s7)
n = len(s7) + 3
ws.cell(row=n, column=1, value="Trait x gene pairs with a TWAS hit (lncRNA and protein-coding genes), by modality: "
                               "expression only, other modalities only, or both.").font = NOTE

# ---- restyle untouched sheets (fonts, header) without changing values
for name in ["S4_Sample size vs eGenes", "S5_GWAS_Sumstats", "S8_Fine_map_balance_meta", "S9_DEG_Sac_vs_Utr"]:
    ws = wb[name]
    for row in ws.iter_rows():
        for c in row:
            if c.value is not None:
                c.font = Font(name=F, size=10, bold=c.font.b, italic=c.font.i)
    ncol = max(c.column for c in ws[1] if c.value is not None)
    header(ws, 1, ncol)
    ws.freeze_panes = "A2"
    if name == "S8_Fine_map_balance_meta":   # second table (isoform QTL fine-mapping) further down
        for row in ws.iter_rows():
            if row[0].value == "SNP" and row[0].row > 1:
                header(ws, row[0].row, ncol)
            if isinstance(row[0].value, str) and row[0].value.startswith("Fine mapping"):
                row[0].font = TITLE
    elif name != "S5_GWAS_Sumstats":
        last = max(c.row for c in ws["A"] if c.value is not None)
        ws.auto_filter.ref = f"A1:{get_column_letter(ncol)}{last}"

# ---- column widths for the data sheets
def autowidth(ws, maxw=60, sample=300):
    for col in ws.iter_cols(min_row=1, max_row=min(ws.max_row, sample)):
        w = max((len(str(c.value)) for c in col if c.value is not None and not str(c.value).startswith("=")), default=6)
        ws.column_dimensions[col[0].column_letter].width = min(maxw, max(8, w + 2))


for name in ["S2_xQTLs_cross_modality", "S6_TWAS_Vestibular", "S7_TWAS hits_expr_non-expr", "S4_Sample size vs eGenes",
             "S8_Fine_map_balance_meta", "S9_DEG_Sac_vs_Utr"]:
    autowidth(wb[name])
autowidth(wb["S5_GWAS_Sumstats"], maxw=70)
for c in wb["S5_GWAS_Sumstats"]["D"] + wb["S5_GWAS_Sumstats"]["E"]:
    c.alignment = Alignment(wrap_text=True, vertical="top")

# ---- README (first sheet)
ws = wb.create_sheet("README", 0)
ws["A1"] = "Supplementary Table 2. QTL mapping and TWAS in human vestibular otolith organs"
ws["A1"].font = Font(name=F, size=14, bold=True)
ws["A2"] = "All QTL and TWAS results: vestibular otolith organs, n = 107 individuals; GRCh38; Ensembl release 110."
ws["A2"].font = NOTE
rows = [("Sheet", "Content"),
        ("S1_N. of phenotypes", "Number of RNA phenotypes per modality, overall and by gene biotype (lncRNA and protein-coding genes)."),
        ("S2_xQTLs_cross_modality", "Independent cis-QTLs from cross-modality mapping (tensorQTL), all genes, sorted by pval_beta. "
                                    "phenotype_id = modality:phenotype."),
        ("S3_N._of_sig_xqtls_xgenes", "Number of significant independent QTL pairs, genes and variants per modality, separate vs cross-modality mapping."),
        ("S4_Sample size vs eGenes", "Sample size vs number of significant eGenes: vestibular tissue vs GTEx tissues."),
        ("S5_GWAS_Sumstats", "GWAS summary statistics used for TWAS: trait, file, sample size and source."),
        ("S6_TWAS_Vestibular", "TWAS hits (FUSION) for all traits and modalities, sorted by COLOC.PP4. Columns: see TWAS_dictionary_columns."),
        ("TWAS_dictionary_columns", "Column descriptions for S6."),
        ("S7_TWAS hits_expr_non-expr", "Trait x gene pairs with a TWAS hit, classified as expression only, other modalities only, or both."),
        ("S8_Fine_map_balance_meta", "Fine-mapping (credible sets) at the chr19 ZNF91 locus: balance GWAS meta-analysis and isoform-ratio QTLs."),
        ("S9_DEG_Sac_vs_Utr", "Differential gene expression, saccule vs utricle.")]
for i, (a, b) in enumerate(rows, 4):
    ws.cell(row=i, column=1, value=a); ws.cell(row=i, column=2, value=b)
    ws.cell(row=i, column=1).font = ws.cell(row=i, column=2).font = BODY
header(ws, 4, 2)
r = 4 + len(rows) + 1
log = ["Changes in this version (2026-10-14) — internal note, remove before submission:",
       "• RNA editing phenotypes were re-generated after fixing a bug in Pantry's assemble_bed.py: sites with no coverage in a sample "
       "(0/0) were set to 100% edited instead of missing, so the >40% missing filter and mean imputation never applied.",
       "• Host-gene expression was then regressed out of each editing site (log2 TPM + 1 of all host genes; residuals re-normalized), "
       "as in the RNA editing QTL pipeline / Li et al. 2022, before QTL mapping, heritability, TWAS and colocalization.",
       "• Re-computed with the fixed RNA editing phenotypes (same 107 individuals and pipeline): S1 (RNA editing rows), S2, S3 "
       "(RNA editing rows and all cross-modality counts), S6 (RNA editing hits), S7. S4, S5, S8 and S9 are unchanged.",
       "• TWAS.Z signs: hits from the original weights have enet/top1 signs corrected (Pantry FUSION.compute_weights.R A2-dosage bug); "
       "the new RNA editing weights were computed with the fixed script. Verified: TWAS.Z sign = eQTL direction x GWAS direction for all scored hits.",
       "• Trait 'balance_meta1' renamed 'balance_meta' in S6 to match S5 and S7. Sheets restyled; S3 headers clarified; TWAS dictionary matched to S6 columns."]
for k, line in enumerate(log):
    c = ws.cell(row=r + k, column=1, value=line)
    c.font = Font(name=F, size=10, bold=(k == 0), color="C00000" if k == 0 else "000000")
ws.column_dimensions["A"].width, ws.column_dimensions["B"].width = 30, 120

# ---- Arial everywhere (default style + every written cell), keeping bold/italic/size/color
wb._named_styles["Normal"].font = Font(name=F, size=10)
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if c.value is not None and c.font.name != F:
                f = c.font
                c.font = Font(name=F, size=f.sz if f.sz and f.sz != 11 and f.sz != 12 else 10, bold=f.b, italic=f.i,
                              color=f.color)

wb.active = 0
wb.calculation.fullCalcOnLoad = True  # no LibreOffice here to cache formula results; Excel recalculates on open
wb.save(OUT)
print("wrote", OUT)
print("S1 RNA editing:", [x for x in s1 if x[0] == "RNA editing"])
print("S2 rows:", len(s2), "| S6 rows:", len(s6), s6["RNA.MODALITY"].value_counts().to_dict(), "| S7 rows:", len(s7), s7.modality_hits.value_counts().to_dict())
print("S3:", {m: s3[m] for m in MODS})
