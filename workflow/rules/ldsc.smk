import yaml
import os

# with open("config/ldsc_checksums.yaml") as _f:
#     _ldsc_checksums = yaml.safe_load(_f)

# def lookup_checksum(path):
#     """Return expected MD5 for a file, keyed by basename. Passed to ensure()."""
#     return _ldsc_checksums[os.path.basename(path)]

_ZENODO = "https://zenodo.org/records/10515792/files"

rule download_ldsc_resources:
    """Download hg38 BaselineLD annotation LD scores from Zenodo."""
    output:
        "resources/ldsc/archives/GRCh38.tgz"
    params:
        url = f"{_ZENODO}/GRCh38.tgz?download=1"
    shell:
        "wget -q -O {output} '{params.url}'"

rule download_ldsc_list_of_rsids:
    """NB: Zenodo Gazal record has a file like this, but it doesn't have alleles"""
    output:
        "resources/ldsc/w_hm3.snplist"
    params:
        url = f"https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist"
    shell:
        "wget -q -O {output} '{params.url}'"

rule download_ldsc_frq:
    output:
        "resources/ldsc/1000G_Phase3_frq.tgz"
    params:
        url = f"{_ZENODO}/1000G_Phase3_frq.tgz?download=1"
    shell:
        "wget -q -O {output} '{params.url}'"

rule unzip_gwas_sumstats:
    """Unzip GWAS summary statistics if needed."""
    input:
        "results/pred_ld/{trait}/imputed/EUR/r2_0_2/maf_0_005/combined.tsv.gz"
    output:
        temp("results/pred_ld/{trait}/imputed/EUR/r2_0_2/maf_0_005/combined.tsv")
    resources:
        runtime = 5
    shell:
        "gunzip -c {input} > {output}"

rule munge_sumstats:
    """Convert GWAS summary statistics to LDSC .sumstats.gz format."""
    input:
        gwas     = rules.unzip_gwas_sumstats.output,
        snp_list = "resources/ldsc/w_hm3.snplist"
    output:
        "results/ldsc/sumstats/{trait}.sumstats.gz"
    params:
        rsid_col = config['col_names']['rsid'],
        alt_col  = config['col_names']['effect_allele'],
        ref_col  = config['col_names']['other_allele'],
        beta_col = config['col_names']['beta'],
        samples = lambda w: config.get('gwas_datasets').get(w.trait).get('samples'),
        out_prefix = subpath(output[0], strip_suffix = ".sumstats.gz")
    resources:
        runtime = 30
    shell:
        """
        python /opt/ldsc/munge_sumstats.py \
            --sumstats {input.gwas} \
            --snp {params.rsid_col} \
            --a1 {params.alt_col} \
            --a2 {params.ref_col} \
            --signed-sumstats {params.beta_col},0 \
            --N {params.samples} \
            --merge-alleles {input.snp_list} \
            --out {params.out_prefix}
        """

rule run_ldsc:
    """Run stratified LDSC with the BaselineLD model to estimate per-annotation tau coefficients. Need to do my own MAF filtering as I'm not using the LDSC frq files."""
    input:
        # TODO list the input files properly
        sumstats     = rules.munge_sumstats.output,
        # annot_done   = "resources/ldsc/baselineLD_hg38/.done",
        # weights_done = "resources/ldsc/weights/.done",
    output:
        results = "results/ldsc/{trait}/ldsc.results",
        log     = "results/ldsc/{trait}/ldsc.log"
    params:
        annot_prefix   = "resources/ldsc/GRCh38/baselineLD_v2.2/baselineLD.",
        weights_prefix = "resources/ldsc/GRCh38/weights/weights.hm3_noMHC.",
        frq_prefix = "resources/ldsc/1000G_Phase3_frq/1000G.EUR.QC.",
        out_prefix     = "results/ldsc/{trait}/ldsc",
        case_proportion = lambda w: config.get('gwas_datasets').get(w.trait).get('case_proportion'),
        pop_prevalence = lambda w: config.get('gwas_datasets').get(w.trait).get('pop_prevalence')
    resources:
        runtime = 60,
        mem_mb  = 16000
    shell:
        """
        python /opt/ldsc/ldsc.py \
            --h2 {input.sumstats} \
            --frqfile-chr {params.frq_prefix} \
            --ref-ld-chr {params.annot_prefix} \
            --w-ld-chr {params.weights_prefix} \
            --samp-prev {params.case_proportion} \
            --pop-prev {params.pop_prevalence} \
            --overlap-annot \
            --print-coefficients \
            --out {params.out_prefix}
        """

rule compute_expected_chi2:
    """Compute per-SNP expected chi-squared from LDSC tau coefficients and BaselineLD LD scores."""
    input:
        results    = rules.run_ldsc.output.results,
        sumstats     = rules.munge_sumstats.output[0]
    output:
        "results/ldsc/{trait}/expected_chi2.tsv.gz"
    params:
        annot_prefix = "resources/ldsc/GRCh38/baselineLD_v2.2/baselineLD.",
        chroms       = list(range(1, 23)),
        samples = lambda w: config.get('gwas_datasets').get(w.trait).get('samples'),
        case_proportion = lambda w: config.get('gwas_datasets').get(w.trait).get('case_proportion')
    threads: 1
    resources:
        runtime = 30
    script: script_path("ldsc/compute_expected_chi2.py")
