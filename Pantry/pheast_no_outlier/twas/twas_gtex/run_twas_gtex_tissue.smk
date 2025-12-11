### Running FUSION TWAS using the GTEx tissues weights produced by Pantry Pheast 
import pandas as pd


# GTEx tissues (dir names)
gtex_tissue = ['ADPSBQ', 'ARTAORT', 'BREAST', 'BRNCDT', 'BRNCTXA', 'BRNHPT', 'BRNSNG', 'CLNTRN', 'ESPMSL', 'HRTAA', 'LCL', 'MSCLSK', 'PNCREAS',
'SKINNS', 'SNTTRM', 'TESTIS', 'VAGINA', 'ADPVSC', 'ARTCRN', 'BRNACC', 'BRNCHA', 'BRNCTXB', 'BRNNCC', 'BRNSPC', 'ESPGEJ', 'FIBRBLS',
'HRTLV', 'LIVER', 'NERVET', 'PRSTTE', 'SKINS', 'SPLEEN', 'THYROID', 'WHLBLD', 'ADRNLG', 'ARTTBL', 'BRNAMY', 'BRNCHB', 'BRNHPP', 'BRNPTM',
'CLNSGM', 'ESPMCS', 'GEUVADIS', 'KDNCTX', 'LUNG', 'OVARY', 'PTTARY', 'SLVRYG', 'STMACH', 'UTERUS']


# RNA Modalities
modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability'] # , 'RNA_editing'

# Subset for test run
# gtex_tissue = ['BRNAMY']
# modalities = ['expression', 'splicing']


# Read GWAS metadata
trait_df = pd.read_csv('gwas_sumstats/gwas_metadata.tsv', sep='\t')
traits = trait_df['trait'].tolist()
gwas_n = {traits[i]: n for i, n in enumerate(trait_df['sample_size'].tolist())}


# Choose chromosomes format (chr1 or just 1)
chroms = [f'chr{i + 1}' for i in range(22)] # chr1-chr22
# chroms = [f'{i + 1}' for i in range(22)] # 1-22


# chr1, chr2, and chr3 need more memory
def mem_by_chr(wildcards, attempt):
    if wildcards.chrom in ["chr1", "chr2", "chr3"]:
        return 128000   # 128 GB
    else:
        return 64000    # 64 GB for others


wildcard_constraints:
    modality = "[^.]+",


rule all:
    input:
        expand("results/{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv",
            modality=modalities, trait=traits, chrom=chroms, dataset=gtex_tissue),
            # expand("results/{dataset}/output/{modality}/twas_hits.tsv", dataset=gtex_tissue, modality=modalities),
            expand("results/{dataset}/output/twas_hits.tsv", dataset=gtex_tissue)


rule assoc_test_chr:
    input:
        sumstats = "gwas_sumstats/{trait}_hg38_twas_input.chr_b38.sumstats.gz",
        weights = "weights/{dataset}/twas/{modality}.pos",
        ref_ld = lambda wildcards: expand("ld_ref/1kg_ldref_hg38_{chrom}.{ext}", chrom=[wildcards.chrom], ext=['bed', 'bim', 'fam']) # If chr1-chr22
        # ref_ld = lambda wildcards: expand("ld_ref/1kg_ldref_hg38_chr{chrom}.{ext}", chrom=[wildcards.chrom], ext=['bed', 'bim', 'fam']) # If 1-22
    output:
        "results/{dataset}/intermediate/{modality}/{trait}/{modality}.{trait}.{chrom}.tsv"
    params:
        weights_dir = "weights/{dataset}/twas",
        ref_ld_prefix = "ld_ref/1kg_ldref_hg38_", # If chr1-chr22
        # ref_ld_prefix = "ld_ref/1kg_ldref_hg38_chr", # If 1-22
        coloc_p = "5e-8", # Will be divided by number of RNA phenotype going to be tested
        gwas_n = lambda w: gwas_n[w.trait]
    threads: 1
    resources:
        # mem_mb = 64000,
        mem_mb = mem_by_chr,
        runtime = "4h"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/intermediate/{wildcards.modality}/{wildcards.trait}
        Rscript ../../scripts/fusion_twas/FUSION.assoc_test.R \
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
        expand('results/{dataset}/intermediate/{{modality}}/{{trait}}/{{modality}}.{{trait}}.{chrom}.tsv',
        dataset=gtex_tissue, modality=modalities, trait=traits, chrom=chroms)
    output:
        out = 'results/{dataset}/output/{modality}/{modality}.{trait}.tsv'
    params:
        output_dir = 'results/{dataset}/output/{modality}',
    group: 'twas'
    threads: 1
    resources:
        mem_mb = 32000,
        runtime = "4h"
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
        expand('results/{dataset}/output/{modality}/{modality}.{trait}.tsv',
                dataset=gtex_tissue, modality=modalities, trait=traits)
    output:
        out = 'results/{dataset}/output/twas_hits.tsv'
    params:
        traits = ' '.join(traits),
        modalities = ' '.join(modalities),
        n_modalities = len(modalities),
    threads: 1
    resources:
        mem_mb = 32000,
        runtime = "4h"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/output
        head -n1 {input[0]} \
            | sed 's/^/TRAIT\tMODALITY\t/' \
            > {output.out}
        for trait in {params.traits}; do
            for modality in {params.modalities}; do
                tail -n+2 results/{wildcards.dataset}/output/$modality/$modality.$trait.tsv \
                    | awk '$18 < 5e-8 / {params.n_modalities}' \
                    | sed "s/^/$trait\t$modality\t/" \
                    >> {output.out}
            done
        done
        """


### END


