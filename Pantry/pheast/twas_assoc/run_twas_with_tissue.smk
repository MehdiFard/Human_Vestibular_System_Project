# Template for running FUSION TWAS using the weights produced by Pantry Pheast
# This example is based on GWAS data from https://zenodo.org/records/3629742

import pandas as pd

datasets = ['vestibular']
# datasets = ['saccule', 'utricle', "vestibular"]
modalities = ['alt_polyA', 'alt_TSS', 'RNA_editing', 'expression', 'isoforms', 'splicing', 'stability']
# modalities = ['expression']
trait_df = pd.read_csv('../input/gwas_sumstat/gwas_metadata.tsv', sep='\t')
traits = trait_df['trait'].tolist()
gwas_n = {traits[i]: n for i, n in enumerate(trait_df['sample_size'].tolist())}
# chroms = [f'chr{i + 1}' for i in range(22)]
chroms = [f'{i + 1}' for i in range(22)]

wildcard_constraints:
    modality = "[^.]+",
    dataset = "[a-z]+"
    # dataset = "[a-zA-Z0-9_-]+"


rule all:
    input:
        expand("{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv",
            modality=modalities, trait=traits, chrom=chroms, dataset=datasets),
            expand("{dataset}/output/{modality}/twas_hits.tsv", dataset=datasets, modality=modalities),


rule assoc_test_chr:
    input:
        # sumstats = "../input/gwas_sumstat/{trait}.sumstats",
        sumstats = "../../../gwas_sumstats_prep/output/{trait}_hg38_twas_input.sumstats.gz",
        weights = "../intermediate/{dataset}/twas/{modality}.pos",
        ref_ld = lambda wildcards: expand("../input/LDREF/1kg_ldref_hg38_chr{chrom}.{ext}", chrom=[wildcards.chrom], ext=['bed', 'bim', 'fam'])
    output:
        "{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv"
    params:
        weights_dir = "../intermediate/{dataset}/twas/",
        ref_ld_prefix = "../input/LDREF/1kg_ldref_hg38_chr",
        coloc_p = "5e-8", # Will be divided by number of RNA phenotype going to be tested
        gwas_n = lambda w: gwas_n[w.trait]
    resources:
        mem_mb = 64000
    shell:
        """
        mkdir -p {wildcards.dataset}/intermediate/{wildcards.modality}/{wildcards.trait}
        Rscript ../scripts/fusion_twas/FUSION.assoc_test.R \
            --sumstats {input.sumstats} \
            --weights {input.weights} \
            --weights_dir {params.weights_dir} \
            --ref_ld_chr {params.ref_ld_prefix} \
            --chr {wildcards.chrom} \
            --coloc_P {params.coloc_p} \
            --GWASN {params.gwas_n} \
            --out {output}
        """


rule assemble_assoc_tests:
    input:
        expand('{dataset}/intermediate/{{modality}}/{{trait}}/{{modality}}.{{trait}}.{chrom}.tsv',
        dataset=datasets, modality=modalities, trait=traits, chrom=chroms)
    output:
        out = '{dataset}/output/{modality}/{modality}.{trait}.tsv'
    params:
        output_dir = '{dataset}/output/{modality}',
    group: 'twas'
    shell:
        """
        mkdir -p {params.output_dir}
        head -n1 {input[0]} | cut -f3- > {output.out}
        for fname in {input}; do
            tail -n+2 $fname | cut -f3- >> {output.out}
        done
        """

rule assemble_hits_all_traits_modalities:
    input:
        expand('{dataset}/output/{modality}/{modality}.{trait}.tsv',
        dataset=datasets, modality=modalities, trait=traits)
    output:
        out = '{dataset}/output/{modality}/twas_hits.tsv'
    params:
        traits = ' '.join(traits),
        modalities = ' '.join(modalities),
        n_modalities = len(modalities),
    shell:
        """
        head -n1 {input[0]} \
            | sed 's/^/TRAIT\tMODALITY\t/' \
            > {output.out}
        for trait in {params.traits}; do
            for modality in {params.modalities}; do
                tail -n+2 {wildcards.dataset}/output/$modality/$modality.$trait.tsv \
                    | awk '$18 < 5e-8 / {params.n_modalities}' \
                    | sed "s/^/$trait\t$modality\t/" \
                    >> {output.out}
            done
        done
        """

### END