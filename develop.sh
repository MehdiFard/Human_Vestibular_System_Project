#!/bin/bash
srun -p cpu-gdml-investor -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A cpu-gdml-investor --pty bash -i
mamba activate Env1


# Another source for imputation:
https://www.protocols.io/view/genotype-imputation-workflow-v3-0-e6nvw78dlmkj/v2?step=1

# For liftover:
https://www.protocols.io/view/genotyping-chip-data-lift-over-to-reference-genome-n2bvjmbpvk5w/v2

# Make soft link to genotypes
ln -s /data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/genotype/* data/genotype/

rsync -avh  /data/hps/assoc/private/gdml/user/sesmae/projects/cybertron_data/ref_gtf_indx_vcfs/1k_genomes_ref_vcf_panel_phase3_hg38/ALL.chr*.{gz,gz.tbi} \
            /data/hps/assoc/private/gdml/user/sesmae/refs_gtfs_indx_vcfs/vcfs/1k_genomes_ref_vcf_panel_phase3_hg38/

snakemake -j 16 --printshellcmds -n
snakemake -j 16 --printshellcmds --scheduler greedy -n


samtools view -H data/bam/FTL_834_Aligned.sortedByCoord.out.bam | grep '^@RG' | head


plink \
    --file ../data/genotype/batch12/20240129_hg19 \
    --merge ../data/genotype/batch3/20250515_hg19 \
    --chr 1-22 \
    --maf 0.01 \
    --snps-only just-acgt \
    --nonfounders \
    --make-bed \
    --out data/geno_merged_hg19

plink2 \
    --bfile data/geno_merged_hg19 \
    --set-all-var-ids @_#_\$r_\$a \
    --make-bed \
    --output-chr chrM \
    --out data/geno_merged_hg19_clean


plink2 \
        --bfile data/geno_merged_hg19_clean \
        --fa src/hg19.fa \
        --ref-from-fa force \
        --max-alleles 2 \
        --nonfounders \
        --output-chr chrM \
        --set-all-var-ids @_#_\$r_\$a \
        --recode vcf id-paste=iid bgz \
        --out data/geno_merged_hg19_clean


bcftools query -f '%ID\n' data/geno_merged_hg19_clean.vcf.gz | grep -v '^\.' | head
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' data/geno_merged_hg19_clean.vcf.gz | head


bcftools query -f '%ID\n' results/geno_lifted_hg38.vcf.gz | grep -v '^\.' | head

chr1	886546	chr1_821926_G_C	C	G
chr1	904947	chr1_840327_G_A	G	A

# uniq works, it keeps one ID
# still have duplicated SNP id in final VCF after liftover, check how to unique before phasing


bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' results/geno_lifted_hg38.vcf.gz | \
bcftools norm -d exact -Ou | \
bcftools view -Oz -o results/geno_lifted_hg38_deduped_renamed.vcf

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' results/geno_lifted_hg38_deduped_renamed.vcf.gz | head


EAGEL="/data/hps/assoc/private/gdml/user/sesmae/projects/cybertron_data/tools/Eagle_v2.4.1/eagle"

cp -r /data/hps/assoc/private/gdml/user/sesmae/projects/cybertron_data/tools/Eagle_v2.4.1 /data/hps/assoc/private/gdml/user/sesmae/tools/
cp /data/hps/assoc/private/gdml/user/sesmae/projects/cybertron_data/tools/Eagle_v2.4.1.tar.gz /data/hps/assoc/private/gdml/user/sesmae/tools/

EAGEL="/data/hps/assoc/private/gdml/user/sesmae/tools/Eagle_v2.4.1/eagle"




bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' src/ref_panel/ALL.chr22_GRCh38.genotypes.20170504.bcf | head


bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' src/ref_panel/ALL.chr22_GRCh38.genotypes.20170504.bcf | head

gmap="/data/hps/assoc/private/gdml/user/sesmae/tools/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

zcat $gmap | head -n 10

zcat src/genetic_map_hg38_withX.txt.gz | head -n 10

zcat src/genetic_map_hg38_withX.txt.gz | awk 'NR==1{print $0; next} { $1="chr"$1; print }' OFS='\t' | gzip > src/genetic_map_hg38_withX_chr.txt.gz

zcat src/genetic_map_hg38_withX_chr.txt.gz | head -n 10

zcat results/phased_vcf/geno_lifted_phased_hg38_chr22.vcf.gz > chr22.vcf

MINIMAC="/data/hps/assoc/private/gdml/user/sesmae/tools/minimac4-4.1.6-Linux-x86_64/bin/minimac4"



bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' results/imputed_vcf/geno_lifted_phased_imputed_hg38_chr1.vcf.gz | head

zcat results/imputed_vcf/geno_lifted_phased_imputed_hg38_chr1.vcf.gz > chr1_imputed.vcf


comm -12 <(sort all_rnas.txt) <(sort extra.txt) | wc -l

ln -s /data/hps/assoc/private/gdml/user/sesmae/projects_data/geuvadis_data/geuvadis_rnaseq/* /data/hps/assoc/private/gdml/user/sesmae/projects/cybertron_data/project02_rsem_kallisto_comparison/data/quantification_proj_data/RawData/RNAseqData/

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' phasing_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz | head
6722694


# Rule: Link all RNA-seq FASTQ files to ../data/rna_fq/
ln -sf /data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/rnaseq_fastq/* ../data/rna_fq/




us frick master


sbatch mapping.slurm
2322720

snakemake -j 16 -R map_star_FTL_500_FTL_501 -n

snakemake -j 16 data/bam/FTL_500_501.Aligned.sortedByCoord.out.bam --printshellcmds -n
