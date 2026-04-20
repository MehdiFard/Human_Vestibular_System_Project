rule draw_locuszoomr_plot_for_chr18_associations:
    input:
        gwas = rules.run_sffdr.output.tsv,
        edb = rules.download_ensembl_db.output
    output:
        "results/sffdr/locuszoomr/chr18_associations.png"
    params:
        p_col         = 'fp',
        snps_to_label = ['rs4800996', 'rs4456576']
    script: script_path("figures/draw_DCC_TCF4_lz_plot.R")
