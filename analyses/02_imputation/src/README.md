# 02_imputation/src (not in git)

Downloaded / prepared by the Snakefile:
- `1kgp_phase3_hg38/`: 1000 Genomes phase 3 phased reference panel, GRCh38 (links to config.yaml `refs.kgp_phase3_vcf_dir`,
  source `refs.kgp_phase3_url`), converted to BCF per chromosome
- `1kgp_phase3_hg38_msav/`: the same panel in Minimac4 `.msav` format
- `genetic_map/genetic_map_hg38_withX.txt.gz` (+ `_chr` version): Eagle genetic map (`refs.eagle_map_url`)
- `grch38_noalt/`: GRCh38 no-alt analysis set FASTA (+ `.fai`) for allele checks (`refs.grch38_noalt_url`)
