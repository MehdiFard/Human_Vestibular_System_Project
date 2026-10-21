# =====================================================================================================================
# Human vestibular (otolith organ) genetics project: main workflow.
# Each analysis in analyses/ has its own Snakefile and can be run alone; this file runs them in order,
# then builds the paper figures (manuscript/figures/Snakefile) and supplementary tables (manuscript/tables/Snakefile).
#   snakemake --profile profiles/slurm -n                         # everything (dry run)
#   snakemake --profile profiles/slurm 08_expression_gtex         # one analysis
#   snakemake -s manuscript/figures/Snakefile --profile profiles/slurm fig3  # one figure
# Shared paths and parameters: config/config.yaml.
# =====================================================================================================================
configfile: "config/config.yaml"
from pathlib import Path
if Path("config/config_lab.yaml").exists():   # lab-specific settings, not in git (overrides config.yaml)
    configfile: "config/config_lab.yaml"

SMK = f"{config['conda_prefix']}/{config['envs']['pheast']}/bin/snakemake"      # Snakemake 9
ANALYSES = ["01_liftover", "02_imputation", "03_sample_matching", "04_gwas_sumstats", "05_pantry", "06_twas",
            "07_fine_mapping", "08_expression_gtex", "09_dge", "10_sffdr"]
AFTER = {"02_imputation": ["01_liftover"], "03_sample_matching": ["02_imputation"], "05_pantry": ["02_imputation"],
         "06_twas": ["04_gwas_sumstats", "05_pantry"], "07_fine_mapping": ["04_gwas_sumstats", "05_pantry"],
         "08_expression_gtex": ["05_pantry"], "09_dge": ["05_pantry"], "10_sffdr": ["04_gwas_sumstats"]}
CMD = {s: f"{SMK} -s {f} --profile {config['profile']}" for s, f in
       [(a, f"analyses/{a}/Snakefile") for a in ANALYSES] + [("figures", "manuscript/figures/Snakefile"), ("tables", "manuscript/tables/Snakefile")]}
# sfFDR (T. Willis) keeps its own layout, profile and container; it runs inside its folder with Snakemake 9.25.
# By default its results/ links to T. Willis's run (read-only, ~650 GB) and is only checked; with sffdr.rerun: true
# the workflow is run into a local results/ (replace the link with a folder first, see analyses/10_sffdr/README.md).
SFFDR_TARGETS = " ".join(config["sffdr"]["targets"])
CMD["10_sffdr"] = (f"cd analyses/10_sffdr && {config['conda_prefix']}/{config['envs']['sffdr_smk']}/bin/snakemake "
                   f"-s workflow/Snakefile --profile " + config["sffdr"]["profile"] + f" {SFFDR_TARGETS}" if config["sffdr"]["rerun"] else
                   "test -s analyses/10_sffdr/results/sffdr/results.tsv.gz")
AFTER.update({"figures": ANALYSES, "tables": ["05_pantry", "06_twas", "figures"]})

localrules: all


rule all:
    input: [f"logs/steps/{s}.done" for s in ANALYSES + ["figures", "tables"]]


for _s in ANALYSES + ["figures", "tables"]:
    rule:
        name: _s
        localrule: True
        input: [f"logs/steps/{d}.done" for d in AFTER.get(_s, [])]
        output: touch(f"logs/steps/{_s}.done")
        params: cmd=CMD[_s]
        shell: "{params.cmd}"
