import pandas as pd
import os

# Rule to download TOP-LD reference data for a specific population and chromosome
rule download_topld_reference:
    """
    NB: European-ancestry reference
    """
    output:
        "resources/pred_ld/TOP_LD/EUR/SNV/EUR_{chr}_no_filter_0.2_1000000_LD.csv.gz"
    params:
        url = "http://topld.genetics.unc.edu/downloads/downloads/EUR/SNV/EUR_{chr}_no_filter_0.2_1000000_LD.csv.gz"
    resources:
        runtime = 30
    shell:
        """
        wget -O {output} {params.url}
        """

rule download_topld_annotation:
    output:
        "resources/pred_ld/TOP_LD/EUR/SNV/EUR_{chr}_no_filter_0.2_1000000_info_annotation.csv.gz"
    params:
        url = "http://topld.genetics.unc.edu/downloads/downloads/EUR/SNV/EUR_{chr}_no_filter_0.2_1000000_info_annotation.csv.gz"
    resources:
        runtime = 30
    shell:
        """
        wget -O {output} {params.url}
        """

# Rule to convert harmonised GWAS to PRED-LD input format
rule convert_gwas_to_pred_ld_format:
    input:
        "resources/gwas/{trait_to_impute}_formatted_hg38_ChrPosId.tsv.gz"
    output:
        "results/pred_ld/{trait_to_impute}/input/{chr}.txt"
    params:
        chr_col = config.get('col_names').get('chromosome'),
        bp_col = config.get('col_names').get('base_pair_location'),
        ref_col = config.get('col_names').get('other_allele'),
        alt_col = config.get('col_names').get('effect_allele'),
        rsid_col = config.get('col_names').get('rsid'),
        beta_col = config.get('col_names').get('beta'),
        se_col = config.get('col_names').get('se')
    threads: 4
    resources:
        runtime = 15
    script: script_path("pred_ld/convert_to_pred_ld_format.R")

# Rule to run PRED-LD imputation on a single chromosome
rule run_pred_ld_imputation:
    input:
        gwas = "results/pred_ld/{trait_to_impute}/input/{chr}.txt",
        topld_reference = rules.download_topld_reference.output,
        topld_annotation = rules.download_topld_annotation.output
    output:
        multiext("results/pred_ld/{trait_to_impute}/imputed/EUR/r2_{r2threshold}/maf_0_{maf}/", imputed = "imputation_results_{chr}.txt", info = "LD_info_TOP_LD_{chr}.txt")
    params:
        r2threshold = lambda w: float(w.r2threshold.replace("_", ".")),
        maf = lambda w: float(f"0.{w.maf}"),
        ref_dir = "resources/pred_ld",
        output_dir = "results/pred_ld/{trait_to_impute}/imputed/EUR/r2_{r2threshold}/maf_0_{maf}"
    threads: 8
    resources:
        runtime = 120
    shell:
        """
        python -m pred_ld \
            --file-path {input.gwas} \
            --r2threshold {params.r2threshold} \
            --pop EUR \
            --maf {params.maf} \
            --ref TOP_LD \
            --ref-dir {params.ref_dir} \
            --output-dir {params.output_dir}
        """

# Rule to validate REF2/ALT2 alleles in TOP-LD LD_info output against Ensembl hg38
rule validate_topld_alleles:
    input:
        rules.run_pred_ld_imputation.output.info
    output:
        "results/pred_ld/{trait_to_impute}/imputed/EUR/r2_{r2threshold}/maf_0_{maf}/validate_topld_alleles_{chr}.log"
    params:
        n_snps = 1000
    resources:
        runtime = 15,
        gnomad_api_calls = 1
    script: script_path("pred_ld/validate_topld_alleles.py")

# Rule to concatenate per-chromosome PRED-LD imputation results and return harmonised GWAS columns
rule convert_pred_ld_to_harmonised_format:
    input:
        expand(
            rules.run_pred_ld_imputation.output.imputed,
            chr=[f"chr{i}" for i in range(21, 22)],
            allow_missing=True
        )
    output:
        "results/pred_ld/{trait_to_impute}/imputed/EUR/r2_{r2threshold}/maf_0_{maf}/combined.tsv.gz"
    params:
        chr_col   = config.get('chr_col'),
        bp_col    = config.get('bp_col'),
        ref_col   = config.get('ref_col'),
        alt_col   = config.get('alt_col'),
        rsid_col  = config.get('rsid_col'),
        beta_col  = config.get('beta_col'),
        se_col    = config.get('se_col'),
        p_col     = config.get('p_col')
    threads: 4
    resources:
        runtime = 30
    script: script_path("pred_ld/convert_to_harmonised_format.py")
