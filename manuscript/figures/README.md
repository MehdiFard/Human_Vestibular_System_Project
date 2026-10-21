# manuscript/figures

`Snakefile` builds every figure from the analysis results: `snakemake -s manuscript/figures/Snakefile --profile profiles/slurm`.
Each `FigN/`, `SuppFigN/` folder holds the code, `panels/` (one file per panel) and, where the figure is assembled in
code, the final figure at the size used in the manuscript (`Fig3_5.8x5.8in`, `Fig4`, `Fig5_6.45x5.02in`,
`SuppFig3`, `SuppFig4_5.5x5.5in`). Fig 1, Fig 2, Supp Fig 1 and Supp Fig 5 are panels only (assembled in Illustrator).
Supp Fig 2 is a partner figure (MAGMA GTEx tissue enrichment), kept as is; Fig 2 also uses a partner panel in
`Fig2/inputs/` (GSEM path diagram).
In git: the PNG versions and `SuppFig2/SuppFig2_enriched_in_gtex.pdf`. Not in git: the other PDFs (figures, panels,
`Fig2/inputs/`), `data/` (cached inputs) and `.rds` plot objects (intermediates the code recreates).
