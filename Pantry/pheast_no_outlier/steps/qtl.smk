import os
from pathlib import Path


### Setup helper functions and parameters for nominal QTL mapping steps
"""
This section (lines 5-36) defines some helper functions and parameters for the nominal QTL mapping steps,
which are used in multiple rules. Comment them out if you don't need to run the nominal QTL mapping.
Also comment out these rules: tensorqtl_cis_nominal and tensorqtl_all_signif.
"""
GROUPED_MODALITIES = {
    "alt_polyA",
    "alt_TSS",
    "cross_modality",
    "isoforms",
    "RNA_editing",
    "splicing",
}


DATASET = config.get("dataset") or os.environ.get("DATASET")
if DATASET is None:
    raise ValueError("Dataset must be provided via --config dataset=... or DATASET env var")


def nominal_outdir(wildcards):
    return output_dir / "qtl" / "nominal" / wildcards.modality


def nominal_prefix_name(wildcards):
    return f"{DATASET}_{wildcards.modality}"


def nominal_prefix_path(wildcards):
    return str(nominal_outdir(wildcards) / nominal_prefix_name(wildcards))
###


# Original Pantry setup
def groups_arg(wildcards, input):
    """Pass the phenotype groups file as a Python arg if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return f'--groups {input.groups}'
    else:
        return ''

def groups_input(wildcards):
    """Include the phenotype groups file as an input if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return pheno_dir / f'{wildcards.modality}.phenotype_groups.txt'
    else:
        return []


rule tensorqtl_cis:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno = multiext(str(geno_prefix), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
    output:
        output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = output_dir / 'qtl',
        groups_arg = groups_arg,
    resources:
        mem_mb = 32000,
        runtime = '12h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.groups_arg} \
            --mode cis
        """

rule tensorqtl_cis_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno = multiext(str(geno_prefix), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
        cis = output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    output:
        output_dir / 'qtl' / '{modality}.cis_independent_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        groups_arg = groups_arg,
    resources:
        mem_mb = 32000,
        runtime = '20h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            {params.groups_arg} \
            --mode cis_independent
        """

rule tensorqtl_trans:
    """Map trans-QTL associations without permutation testing.
    Outputs all variant pairs with TSS distance > 5 Mb, MAF > 0.05, and p < 1e-5.
    """
    input:
        geno = multiext(str(geno_prefix), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
    output:
        output_dir / 'qtl' / '{modality}.trans_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = output_dir / 'qtl',
    resources:
        mem_mb = 32000,
        runtime = '12h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode trans
        """


# Nominal mapping
rule tensorqtl_cis_nominal:
    """
    Run TensorQTL cis_nominal to write nominal cis summary stats for all tested SNPs per phenotype.
    TensorQTL writes multiple files per-chromosome using the provided prefix.
    """
    input:
        geno = multiext(str(geno_prefix), ".bed", ".bim", ".fam"),
        bed  = pheno_dir / "{modality}.bed.gz",
        bedi = pheno_dir / "{modality}.bed.gz.tbi",
        covar = interm_dir / "covar" / "{modality}.covar.tsv",
    output:
        nominal_outdir = directory(output_dir / "qtl" / "nominal" / "{modality}"),
        done = output_dir / "qtl" / "nominal" / "{modality}" / ".cis_nominal.done",
    params:
        geno_prefix = geno_prefix,
        outdir = nominal_outdir,
        out_prefix = nominal_prefix_name,
    resources:
        mem_mb = 32000,
        runtime = "12h",
    shell:
        r"""
        mkdir -p {params.outdir}
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {params.out_prefix} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --mode cis_nominal
        touch {output.done}
        """


rule tensorqtl_all_signif:
    input:
        perm = output_dir / "qtl" / "{modality}.cis_qtl.txt.gz",
        nominal_done = output_dir / "qtl" / "nominal" / "{modality}" / ".cis_nominal.done",
    output:
        output_dir / "qtl" / "nominal" / "{modality}.cis_qtl_all_signif.txt.gz",
    params:
        nominal_prefix = nominal_prefix_path,
        grouped = lambda wc: wc.modality in GROUPED_MODALITIES,
        fdr = 0.05,
    shell:
        r"""
        if [[ "{params.grouped}" == "True" ]]; then
        # if [[ "{wildcards.modality}" != "expression" && "{wildcards.modality}" != "stability" ]]; then
            python3 scripts/tensorqtl_all_signif.py \
                {input.perm} \
                {params.nominal_prefix} \
                {output} \
                --fdr {params.fdr} \
                --groups input/{DATASET}/phenotypes/{wildcards.modality}.phenotype_groups.txt
        else
            python3 scripts/tensorqtl_all_signif.py \
                {input.perm} \
                {params.nominal_prefix} \
                {output} \
                --fdr {params.fdr}
        fi
        """

### END