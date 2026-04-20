rule merge_sffdr_input_with_1kG:
    input:
        primary = config['sffdr_inputs']['primary'],
        informative = config['sffdr_inputs']['informative'],
        onekgp_panel = "results/1kG/hg38/all/eur/snps_only/005/merged.pvar.zst"
    output:
        "results/sffdr/input_on_1kG_snp_ids.txt"
    threads: 8
    resources:
        runtime = 60
    script: script_path("merge_sffdr_input_with_1kG.R")

rule make_pruned_range_from_sffdr_1kG_join:
    input:
        multiext("results/1kG/hg38/all/eur/snps_only/005/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = rules.merge_sffdr_input_with_1kG.output
    output:
        temp(multiext("results/sffdr/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.prune.in'),
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    resources:
        runtime = 90
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --extract {input.ids} --indep-pairwise {wildcards.window_size} 1 {params.r2} --out {params.out_stem}"

rule generate_sffdr_input:
    input:
        primary = config['sffdr_inputs']['primary'],
        prune = "results/sffdr/1000kb_1_0_2/merged.prune.in",
        informative = config['sffdr_inputs']['informative']
    output:
        "results/sffdr/input.tsv"
    threads: 8
    resources:
        runtime = 60
    script: script_path("generate_sffdr_input.R")

rule run_sffdr:
    input:
        rules.generate_sffdr_input.output
    output:
        tsv = "results/sffdr/results.tsv.gz",
        rds = "results/sffdr/results.rds",
        png = "results/sffdr/results.png"
    params:
        seed = 145,
        use_indep_snps = False
    threads: 16
    resources:
        runtime = 120
    script: script_path("sffdr.R")

rule distance_clump_sffdr:
    input:
        rules.run_sffdr.output
    output:
        "results/sffdr/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = False,
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
        p_col = 'fp'
    threads: 16
    resources:
        runtime = 5
    script: script_path("distance_clump.R")

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_lead_snps_with_missense_and_qtl_info_for_sffdr with:
    input:
        lead_snps = rules.distance_clump_sffdr.output
    output:
        "results/sffdr/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_lead_snps_with_nearest_gene_for_sffdr with:
    input:
        lead = rules.annotate_lead_snps_with_missense_and_qtl_info_for_sffdr.output,
        edb = rules.download_ensembl_db.output
    output:
        "results/sffdr/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"

rule finalise_lead_snp_annotations_for_sffdr:
    input:
        lead_snps = rules.annotate_lead_snps_with_nearest_gene_for_sffdr.output,
        sffdr = rules.run_sffdr.output.tsv
    output:
        "results/sffdr/{window_size}_{threshold}_annotated_lead_snps.tsv"
    threads: 8
    resources:
        runtime = 5
    script: script_path("merge_lead_snps_and_sffdr_results.R")

use rule draw_manhattan_with_lead_snp_annotation as draw_manhattan_with_lead_snp_annotation_for_sffdr with:
    input:
        lead_snps = rules.finalise_lead_snp_annotations_for_sffdr.output,
        gwas = rules.run_sffdr.output.tsv
    output:
        "results/sffdr/{window_size}_{threshold}_manhattan_plot.png"
    params:
        p_col = 'fp'

rule get_pruned_set_from_sffdr_input:
    input:
        "results/sffdr/input.tsv"
    output:
        temp("results/sffdr/pruned_input.tsv")
    shell:
        """
        head -n 1 {input} > {output}
        grep "TRUE" {input}  >> {output}
        """

use rule draw_manhattan_with_lead_snp_annotation as plot_pruned_sffdr_input with:
    input:
        gwas = "results/sffdr/pruned_input.tsv"
    output:
        "results/sffdr/pruned_input.png"

rule draw_locuszoomr_plot_for_sffdr:
    input:
        gwas = rules.run_sffdr.output.tsv,
        edb = rules.download_ensembl_db.output
    output:
        "results/sffdr/locuszoomr/{chr}_{start}_{end}.png"
    params:
        p_col         = 'fp',
        snps_to_label = []
    script: script_path("locuszoomr/plot_locus.R")
