# manuscript/figures

`Snakefile` builds every figure from the analysis results: `snakemake -s manuscript/figures/Snakefile --profile profiles/slurm`.
Each `FigN/`, `SuppFigN/` folder holds the code, `panels/` (one file per panel), `source_data/` (data behind each
panel, where made) and the final figure (`.pdf`, `.png`, and the size used in the manuscript, e.g.
`Fig5_6.45x5.02in.pdf`). Fig 2 and Supp Fig 2 also use partner panels in `inputs/` (GSEM, GTEx enrichment).
In git: the PNG versions. Not in git: the PDF versions (figures, panels and the partner panels in `inputs/`), `data/`
(cached inputs), `.rds` plot objects and `render_tmp/` (intermediates the code recreates).
