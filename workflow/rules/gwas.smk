rule download_ensembl_db:
    output:
        "resources/ensembl_113_hsapiens_edb.sqlite"
    params:
        search_tokens = ["EnsDb", "Homo sapiens", "Ensembl 113"],
        genome = "GRCh38",
    script: script_path("download_edb.R")

rule draw_manhattan:
    input:
        gwas = "resources/gwas/{trait}_formatted_hg38_ChrPosId.tsv.gz"
    output:
        gwas = "results/gwas/{trait}_manhattan_plot.png"
    params:
        width = 12,
        height = 4,
        p_col = config['col_names']['p_value']
    script: script_path("plot_gwas_manhattan.R")
