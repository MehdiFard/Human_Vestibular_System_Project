rule download_ensembl_db:
    output:
        "resources/ensembl_113_hsapiens_edb.sqlite"
    params:
        search_tokens = ["EnsDb", "Homo sapiens", "Ensembl 113"],
        genome = "GRCh38",
    script: script_path("download_edb.R")

rule distance_clump_gwas:
    input:
        "resources/gwas/{trait}_formatted_hg38_ChrPosId.tsv.gz"
    output:
        "results/gwas/{trait}/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda wildcards: False, #if wildcards.snp_set == 'sans_mhc' else True,
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
    threads: 16
    resources:
        runtime = 5
    script: script_path("distance_clump.R")

rule annotate_lead_snps_with_missense_and_qtl_info:
    input:
        rules.distance_clump_gwas.output
    output:
        "results/gwas/{trait}/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"
    resources:
        runtime = 30
    script: script_path("annotate_with_missense_and_qtl_info.py")

rule annotate_lead_snps_with_nearest_gene:
    input:
        lead = rules.annotate_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/ensembl_113_hsapiens_edb.sqlite"
    output:
        "results/gwas/{trait}/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"
    resources:
        runtime = 30
    script: script_path("annotate_with_nearest_gene.R")

rule finalise_lead_snp_annotations:
    input:
        rules.annotate_lead_snps_with_nearest_gene.output
    output:
        "results/gwas/{trait}/{window_size}_{threshold}_annotated_lead_snps.tsv"
    resources:
        runtime = 30
    script: script_path("finalise_lead_snp_annotations.R")

rule draw_manhattan_with_lead_snp_annotation:
    input:
        gwas = "resources/gwas/{trait}_formatted_hg38_ChrPosId.tsv.gz",
        lead_snps = rules.finalise_lead_snp_annotations.output
    output:
        "results/gwas/{trait}/annotated_manhattan.png"
    params:
        title = '',
        width = 6,
        height = 4,
        p_col = config['col_names']['p_value']
    threads: 16
    resources:
        runtime = 20
    script:
        script_path("plot_gwas_manhattan.R")
