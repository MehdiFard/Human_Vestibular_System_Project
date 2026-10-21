#!/bin/bash
# Prepare pheast inputs for one dataset (run from analyses/05_pantry/pheast, after analyses/05_pantry/phenotyping has finished):
#   1. cross_modality phenotypes in ../phenotyping/output/<dataset> (7 modalities concatenated, "modality:" prefixes)
#   2. phenotypes renamed from sample IDs to individual IDs, restricted to the analysed individuals
#      (input/<dataset>/<dataset>.individuals.txt: 107 individuals for vestibular, 3 PCA outliers removed)
#   3. RNA editing: host-gene expression regressed out (scripts/regress_rna_editing_host_expression.py; Li et al. 2022);
#      the unregressed RNA_editing and cross_modality go to input/<dataset>/phenotypes_unregressed/
#   4. cross_modality rebuilt with the regressed RNA editing rows (scripts/build_cross_modality.py)
#   5. genotypes from the imputed VCF (MAF 0.05, geno 0.01, HWE 1e-6, biallelic SNPs, autosomes) and covariates
# Usage:  bash prepare_inputs.sh [dataset=vestibular] [imputed_vcf] [--phenotypes-only]
set -euo pipefail
# needs the pantry_pheast environment on PATH (plink2, python3 with pandas); ../Snakefile (rule pheast_inputs) sets it
cd "$(dirname "$0")"
DS=${1:-vestibular}
VCF=${2:-../../02_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz}
ONLY=${3:-}
P=../phenotyping
IN=input/$DS

# 1. cross-modality phenotypes (sample IDs)
(cd $P && bash scripts/combine_modalities_updated_060425.sh --path_to_output output/$DS \
    alt_polyA alt_TSS expression RNA_editing isoforms splicing stability)

# 2. individual IDs + analysed individuals
mkdir -p $IN
[ -s $IN/$DS.individuals.txt ] || { echo "missing $IN/$DS.individuals.txt (analysed individuals, one per line)"; exit 1; }
awk '{ print $1 "\tFTL_" $1 }' $P/input/$DS/$DS.samples.tsv > $IN/$DS.sample_map.txt
python3 scripts/prepare_phenotypes_updated072225.py --indir $P/output/$DS --outdir $IN/phenotypes/ \
    --map $IN/$DS.sample_map.txt --individuals $IN/$DS.individuals.txt

# 3-4. RNA editing ~ host-gene expression; cross_modality with the regressed editing
mkdir -p $IN/phenotypes_unregressed
mv -f $IN/phenotypes/RNA_editing.* $IN/phenotypes/cross_modality.* $IN/phenotypes_unregressed/
python3 scripts/regress_rna_editing_host_expression.py
python3 -W ignore scripts/build_cross_modality.py
[ "$ONLY" = "--phenotypes-only" ] && exit 0

# 5. genotypes and covariates
plink2 --vcf $VCF --keep $IN/$DS.individuals.txt \
    --fa $P/input/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --ref-from-fa force \
    --maf 0.05 --geno 0.01 --hwe 1e-6 --snps-only just-acgt --nonfounders --max-alleles 2 --chr 1-22 \
    --set-all-var-ids @_#_'$r'_'$a' --make-bed --out $IN/${DS}_genotype
[ -s $IN/${DS}_covars.tsv ] || cp ../../03_sample_matching/results/sample_sheets/${DS}_covars.tsv $IN/
echo "inputs ready in $IN"
