rule draw_locuszoomr_plot_for_chr18_associations:
    input:
        gwas = rules.run_sffdr.output.tsv,
        edb = rules.download_ensembl_db.output
    output:
        "results/pub/chr18_associations.png"
    params:
        p_col         = 'fp',
        snps_to_label = ['rs4800996', 'rs4456576']
    script: script_path("figures/draw_DCC_TCF4_lz_plot.R")

rule draw_sffdr_manhattan_plot_for_publication:
    input:
        gwas = rules.run_sffdr.output.tsv,
        lead_snps = "results/sffdr/1000kb_gws_annotated_lead_snps.tsv"
    params:
        p_col = "fp"
    output:
        "results/pub/sffdr_manhattan_plot.png"
    script: script_path("figures/draw_sffdr_manhattan_plot_for_publication.R")

use rule draw_sffdr_manhattan_plot_for_publication as draw_gwas_manhattan_plot_for_publication_with_legend with:
    input:
        gwas = "resources/gwas/{trait}_formatted_hg38_ChrPosId.tsv.gz",
        lead_snps = "results/gwas/{trait}/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/pub/{trait}_manhattan_plot.png"
    params:
        p_col = "P"

rule write_tables_to_xlsx:
    input:
        sffdr_lead_snps = "results/sffdr/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/pub/sffdr_lead_snp_table.xlsx"
    script: script_path("tables/write_tables_in_xlsx.py")
