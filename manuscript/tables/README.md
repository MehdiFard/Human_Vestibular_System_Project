# manuscript/tables

- `SuppTables/Supplementary_Table1_GWAS.xlsx`: GWAS results (phenotype definitions, gene mapping, MAGMA, LDSC, genomic
  SEM, sfFDR, demographics) provided by the collaborators; cleaned layout, values as provided. Not built by the pipeline.
- `SuppTables/Supplementary_Table2_QTL_TWAS.xlsx`: phenotypes, xQTLs, TWAS (`build_supp_table2.py`)
- `manuscript_numbers/manuscript_numbers.tsv`: every number quoted in the text (`manuscript_numbers.py`)

`Snakefile` rebuilds Supplementary Table 2 and the numbers from the analysis results.
