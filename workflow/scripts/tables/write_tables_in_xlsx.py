import pandas as pd

def write_sheet(writer, df, sheet_name):
    df = df.to_excel(writer, sheet_name=sheet_name, index=False)
    wb  = writer.book
    ws  = writer.sheets[sheet_name]
    hdr = wb.add_format({"bold": True})

    for col_num, col_name in enumerate(df.columns):
        ws.write(0, col_num, col_name, hdr)

    ws.freeze_panes(1, 0)
    ws.set_column(0, len(df.columns) - 1, 18)

    return df

sffdr_lead_snps_df = pd.read_csv(snakemake.input.sffdr_lead_snps, sep="\t")
sffdr_lead_snps_df = sffdr_lead_snps_df[["rsID", "CHR", "BP", "SNP", "A1", "A2", "P", "fp", "fpi0", "flfdr", "genes", "P.1", "P.2", "P.3", "P.4", "P.5"]]
sffdr_lead_snps_df.rename(columns={"CHR": "Chromosome", "BP": "Position", "A1": "Effect allele", "A2": "Other allele", "P": "Dizziness p-value", "fp": "Functional p-value", "fpi0": "Functional prior probability of the null", "flfdr": "Functional local FDR","P.1": "Tinnitus p-value", "P.2": "Hearing loss p-value", "P.3": "Anxiety p-value", "P.4": "Depression p-value", "P.5": "Frailty p-value"}, inplace=True)
sffdr_lead_snps_df = sffdr_lead_snps_df[["rsID", "Chromosome", "Position", "Effect allele", "Other allele", "Dizziness p-value", "Functional p-value", "Functional prior probability of the null", "Functional local FDR", "Tinnitus p-value", "Hearing loss p-value", "Anxiety p-value", "Depression p-value", "Frailty p-value"]]

with pd.ExcelWriter(snakemake.output[0], engine="xlsxwriter") as writer:
    write_sheet(writer, sffdr_lead_snps_df, "sfFDR lead SNPs")
