# =====================================================================================================================
# Human vestibular (otolith organ) project: main workflow.
# Each analysis in analyses/ has its own Snakefile and can be run alone; this file runs them in order,
# then builds the paper figures (manuscript/figures/Snakefile) and supplementary tables (manuscript/tables/Snakefile).
#   snakemake --profile profiles/slurm -n                         # everything (dry run)
#   snakemake --profile profiles/slurm 08_expression_gtex         # one analysis
#   snakemake -s manuscript/figures/Snakefile --profile profiles/slurm fig3  # one figure
# Shared paths and parameters: config/config.yaml.
# =====================================================================================================================
# ---- configuration: shared settings, then lab-specific overrides if present
configfile: "config/config.yaml"
from pathlib import Path
if Path("config/config_lab.yaml").exists():   # lab-specific settings, not in git (overrides config.yaml)
    configfile: "config/config_lab.yaml"

# ---- Snakemake used to run each analysis's own Snakefile
SMK = f"{config['conda_prefix']}/{config['envs']['pheast']}/bin/snakemake"      # Snakemake 9

# ---- steps: the analyses in analyses/, then the paper figures and tables
ANALYSES = ["01_liftover", "02_imputation", "03_sample_matching", "04_gwas_sumstats", "05_pantry", "06_twas",
            "07_fine_mapping", "08_expression_gtex", "09_dge", "10_sffdr"]
STEPS = ANALYSES + ["figures", "tables"]

# ---- run order: a step starts once the steps it depends on have finished (each writes logs/steps/<step>.done)
DEPENDS_ON = {"02_imputation": ["01_liftover"], "03_sample_matching": ["02_imputation"], "05_pantry": ["02_imputation"],
              "06_twas": ["04_gwas_sumstats", "05_pantry"], "07_fine_mapping": ["04_gwas_sumstats", "05_pantry"],
              "08_expression_gtex": ["05_pantry"], "09_dge": ["05_pantry"], "10_sffdr": ["04_gwas_sumstats"],
              "figures": ANALYSES, "tables": ["05_pantry", "06_twas", "figures"]}

# ---- command for each step: run its Snakefile with the Slurm profile
CMD = {s: f"{SMK} -s {f} --profile {config['profile']}" for s, f in
       [(a, f"analyses/{a}/Snakefile") for a in ANALYSES] + [("figures", "manuscript/figures/Snakefile"), ("tables", "manuscript/tables/Snakefile")]}
# sfFDR keeps its own layout, Snakemake profile and container, and runs inside its folder (Snakemake >= 9.25).
# sffdr.rerun: false (default) only checks that its results exist, e.g. from an earlier run; true runs the workflow
# (see analyses/10_sffdr/README.md).
SFFDR_TARGETS = " ".join(config["sffdr"]["targets"])
CMD["10_sffdr"] = (f"cd analyses/10_sffdr && {config['conda_prefix']}/{config['envs']['sffdr_smk']}/bin/snakemake "
                   f"-s workflow/Snakefile --profile " + config["sffdr"]["profile"] + f" {SFFDR_TARGETS}" if config["sffdr"]["rerun"] else
                   "test -s analyses/10_sffdr/results/sffdr/results.tsv.gz")

# ---- these rules only launch the step workflows, so they run in this controller job, not as Slurm jobs
localrules: all


# ---- target: every step finished
rule all:
    input: [f"logs/steps/{s}.done" for s in STEPS]


# ---- one rule per step: wait for the steps it depends on
for _s in STEPS:
    rule:
        name: _s
        localrule: True
        input: [f"logs/steps/{d}.done" for d in DEPENDS_ON.get(_s, [])]
        output: touch(f"logs/steps/{_s}.done")
        params: cmd=CMD[_s]
        shell: "{params.cmd}"


### END