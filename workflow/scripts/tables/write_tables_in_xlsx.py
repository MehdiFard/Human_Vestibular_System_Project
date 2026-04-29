import pandas as pd
import re

def write_sheet(writer, df, sheet_name):
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    wb  = writer.book
    ws  = writer.sheets[sheet_name]
    fmt_hdr      = wb.add_format({"bold": True})
    fmt_pvalue   = wb.add_format({"num_format": "0.0E+00"})
    fmt_position = wb.add_format({"num_format": "#,##0"})
    fmt_2dp      = wb.add_format({"num_format": "0.00"})

    def col_fmt(col_name):
        if re.search(r'Functional p-value|p-value|fp-value|\bFDR\b|FDR', col_name, re.IGNORECASE):
            return fmt_pvalue
        if re.match(r'Position', col_name):
            return fmt_position
        if re.match(r'Beta|Standard error|Genetic correlation estimate|Heritability estimate|Functional prior', col_name):
            return fmt_2dp
        return None

    for col_num, col_name in enumerate(df.columns):
        ws.write(0, col_num, col_name, fmt_hdr)
        fmt = col_fmt(col_name)
        ws.set_column(col_num, col_num, 18, fmt)

    ws.freeze_panes(1, 0)

    ws.set_column(0, len(df.columns) - 1, 18)

sffdr_lead_snps_df = pd.read_csv(snakemake.input.sffdr_lead_snps, sep="\t")
sffdr_lead_snps_df = sffdr_lead_snps_df[["rsID", "CHR", "BP", "A1", "A2", "P", "fp", "fpi0", "flfdr", "missense_gene", "qtl_genes", "nearest_gene_name", "distance_bp", "P.1", "P.2", "P.3", "P.4", "P.5"]]
sffdr_lead_snps_df.rename(columns={"CHR": "Chromosome",
 "BP": "Position",
 "A1": "Effect allele",
 "A2": "Other allele",
 "P": "Dizziness p-value",
 "fp": "Functional p-value",
 "fpi0": "Functional prior probability of the null",
 "flfdr": "Functional local FDR",
 "missense_gene": "Missense gene",
 "qtl_genes": "QTL genes",
 "nearest_gene_name": "Nearest genes",
"P.1": "Tinnitus p-value",
 "P.2": "Hearing loss p-value",
 "P.3": "Anxiety p-value",
 "P.4": "Depression p-value",
 "P.5": "Frailty p-value"},
 inplace=True)
sffdr_lead_snps_df = sffdr_lead_snps_df[["rsID",
 "Chromosome",
 "Position",
 "Effect allele",
 "Other allele",
 "Missense gene",
 "QTL genes",
 "Nearest genes",
 "Dizziness p-value",
 "Functional p-value",
 "Functional prior probability of the null",
 "Functional local FDR",
 "Tinnitus p-value",
 "Hearing loss p-value",
 "Anxiety p-value",
 "Depression p-value",
 "Frailty p-value"]]

with pd.ExcelWriter(snakemake.output[0], engine="xlsxwriter") as writer:
    write_sheet(writer, sffdr_lead_snps_df, "sfFDR lead SNPs")
