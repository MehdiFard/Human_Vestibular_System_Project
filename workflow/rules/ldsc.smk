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
    resources:
        runtime = 120
    shell:
        "wget -q -O {output} '{_ZENODO}/GRCh38.tgz?download=1'"

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
    resources:
        runtime = 30
    shell:
        "wget -q -O {output} '{_ZENODO}/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1'"

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
    resources:
        runtime = 30
    shell:
        "wget -q -O {output} '{_ZENODO}/1000G_Phase3_frq.tgz?download=1'"

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
    resources:
        runtime = 10
    shell:
        "wget -q -O {output} '{_ZENODO}/hm3_no_MHC.list.txt?download=1'"
