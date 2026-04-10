# workflow/rules/ldsc.smk
import yaml
import os

with open("config/ldsc_checksums.yaml") as _f:
    _ldsc_checksums = yaml.safe_load(_f)

def lookup_checksum(path):
    """Return expected MD5 for a file, keyed by basename. Passed to ensure()."""
    return _ldsc_checksums[os.path.basename(path)]

_ZENODO = "https://zenodo.org/records/10515792/files"

rule download_ldsc_annotations:
    """Download hg38 BaselineLD annotation LD scores from Zenodo."""
    output:
        ensure("resources/ldsc/archives/GRCh38.tgz", md5=lookup_checksum)
    params:
        url = f"{_ZENODO}/GRCh38.tgz?download=1"
    resources:
        runtime = 120
    shell:
        "wget -q -O {output} '{params.url}'"

rule extract_ldsc_annotations:
    input:
        "resources/ldsc/archives/GRCh38.tgz"
    output:
        touch("resources/ldsc/baselineLD_hg38/.done")
    shell:
        "mkdir -p resources/ldsc/baselineLD_hg38 && tar -xzf {input} -C resources/ldsc/baselineLD_hg38/"

rule download_ldsc_weights:
    """Download 1000G Phase3 hm3 (no MHC) regression weights from Zenodo."""
    output:
        ensure("resources/ldsc/archives/weights.tgz", md5=lookup_checksum)
    params:
        url = f"{_ZENODO}/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1"
    resources:
        runtime = 30
    shell:
        "wget -q -O {output} '{params.url}'"

rule extract_ldsc_weights:
    input:
        "resources/ldsc/archives/weights.tgz"
    output:
        touch("resources/ldsc/weights/.done")
    shell:
        "mkdir -p resources/ldsc/weights && tar -xzf {input} -C resources/ldsc/weights/"

rule download_ldsc_frq:
    """Download 1000G Phase3 EUR allele frequency files from Zenodo."""
    output:
        ensure("resources/ldsc/archives/frq.tgz", md5=lookup_checksum)
    params:
        url = f"{_ZENODO}/1000G_Phase3_frq.tgz?download=1"
    resources:
        runtime = 30
    shell:
        "wget -q -O {output} '{params.url}'"

rule extract_ldsc_frq:
    input:
        "resources/ldsc/archives/frq.tgz"
    output:
        touch("resources/ldsc/frq/.done")
    shell:
        "mkdir -p resources/ldsc/frq && tar -xzf {input} -C resources/ldsc/frq/"

rule download_ldsc_snp_list:
    """Download HapMap3 (no MHC) SNP list from Zenodo."""
    output:
        ensure("resources/ldsc/hm3_no_MHC.list.txt", md5=lookup_checksum)
    params:
        url = f"{_ZENODO}/hm3_no_MHC.list.txt?download=1"
    resources:
        runtime = 10
    shell:
        "wget -q -O {output} '{params.url}'"

rule munge_sumstats:
    """Convert GWAS summary statistics to LDSC .sumstats.gz format."""
    input:
        gwas     = "resources/gwas/{trait}_formatted_hg38_ChrPosId.tsv.gz",
        snp_list = "resources/ldsc/hm3_no_MHC.list.txt"
    output:
        "resources/ldsc/sumstats/{trait}.sumstats.gz"
    params:
        rsid_col = config['col_names']['rsid'],
        alt_col  = config['col_names']['effect_allele'],
        ref_col  = config['col_names']['other_allele'],
        beta_col = config['col_names']['beta'],
        n_col    = config['col_names']['n_col'],
        out_prefix = lambda w, output: str(output[0]).replace('.sumstats.gz', '')
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
            --N-col {params.n_col} \
            --merge-alleles {input.snp_list} \
            --out {params.out_prefix}
        """

rule run_ldsc:
    """Run stratified LDSC with the BaselineLD model to estimate per-annotation tau coefficients."""
    input:
        sumstats     = "resources/ldsc/sumstats/{trait}.sumstats.gz",
        annot_done   = "resources/ldsc/baselineLD_hg38/.done",
        weights_done = "resources/ldsc/weights/.done",
        frq_done     = "resources/ldsc/frq/.done"
    output:
        results = "results/ldsc/{trait}/ldsc.results",
        log     = "results/ldsc/{trait}/ldsc.log"
    params:
        annot_prefix   = "resources/ldsc/baselineLD_hg38/baselineLD.",
        weights_prefix = "resources/ldsc/weights/weights.hm3_noMHC.",
        frq_prefix     = "resources/ldsc/frq/1000G.EUR.",
        out_prefix     = "results/ldsc/{trait}/ldsc"
    resources:
        runtime = 60,
        mem_mb  = 16000
    shell:
        """
        python /opt/ldsc/ldsc.py \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.annot_prefix} \
            --w-ld-chr {params.weights_prefix} \
            --frq-chr {params.frq_prefix} \
            --overlap-annot \
            --print-coefficients \
            --out {params.out_prefix}
        """

rule compute_expected_chi2:
    """Compute per-SNP expected chi-squared from LDSC tau coefficients and BaselineLD LD scores."""
    input:
        results    = "results/ldsc/{trait}/ldsc.results",
        annot_done = "resources/ldsc/baselineLD_hg38/.done"
    output:
        "results/ldsc/{trait}_expected_chi2.tsv.gz"
    params:
        annot_prefix = "resources/ldsc/baselineLD_hg38/baselineLD.",
        chroms       = list(range(1, 23))
    threads: 4
    resources:
        runtime = 30
    script: script_path("ldsc/compute_expected_chi2.py")
