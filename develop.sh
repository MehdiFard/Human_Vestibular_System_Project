#!/bin/bash
# sesmae@login-2.hpc.childrens.sea.kids
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
snakemake -j 16 --printshellcmds --scheduler greedy -R summarize_matching -n

snakemake -j 24 --printshellcmds --scheduler greedy -n


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

sbatch mapping.slurm
2322720

snakemake -j 16 -R map_star_FTL_500_FTL_501 -n

snakemake -j 16 data/bam/FTL_500_501.Aligned.sortedByCoord.out.bam --printshellcmds -n



###---------------------------------------------- Pantry phenotypyong
srun -p cpu-gdml-investor -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A cpu-gdml-investor --pty bash -i
mamba activate pantry

## Saccuel
snakemake       --configfile config.yml \
                --config dataset=saccule \
                --jobs 16 \
                --default-resources mem_mb=64000 cpus=16 \
                --printshellcmds \
                --rerun-incomplete \
                --latency-wait 60 -n

# Submit the job to the cluster - saccule
sbatch --export=ALL,DATASET=saccule --job-name=saccule_phenotyping run_phenotyping.slurm
2672069
2840149
2864909


## Utricle
snakemake       --configfile config.yml \
                --config dataset=utricle \
                --jobs 16 \
                --default-resources mem_mb=64000 cpus=16 \
                --printshellcmds \
                --rerun-incomplete \
                --latency-wait 60 -n

# Submit the job to the cluster - utricle
sbatch --export=ALL,DATASET=utricle --job-name=utricle_phenotyping run_phenotyping.slurm
2672072
2840150
2864915


## Vestibular
snakemake       --configfile config.yml \
                --config dataset=vestibular \
                --jobs 16 \
                --default-resources mem_mb=64000 cpus=16 \
                --printshellcmds \
                --rerun-incomplete \
                --latency-wait 60 -n

# Submit the job to the cluster - utricle
sbatch --export=ALL,DATASET=vestibular --job-name=vestibular_phenotyping run_phenotyping.slurm
2704904
2840151
2845105


# ###----------------------- Making GWAS sumstat
# # Copy gwas result
# cp /data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/gwas_results/balance_meta1.tbl.sorted.filtered.gz \
#         Pantry/pheast/input/gwas_sumstat/

# # Downloding LDSC to make the GWAS sumstat
# git clone https://github.com/bulik/ldsc.git /data/hps/assoc/private/gdml/user/sesmae/tools/ldsc
# # Create the environment
# mamba env create --file /data/hps/assoc/private/gdml/user/sesmae/tools/ldsc/environment.yml
# mamba activate ldsc
# LDSC="/data/hps/assoc/private/gdml/user/sesmae/tools/ldsc/munge_sumstats.py"
# $LDSC -h

# $LDSC --sumstats Pantry/pheast/input/gwas_sumstat/balance_meta1.tbl.sorted.filtered.gz \
#                         --out Pantry/pheast/input/gwas_sumstat/balance_meta1 \
#                         --N 420000
# gunzip Pantry/pheast/input/gwas_sumstat/balance_meta1.sumstats.gz

# # Making the gwas metadata file
# awk 'BEGIN { print "trait\tsample_size\nbalance_meta1\t420000" }' > Pantry/pheast/input/gwas_sumstat/gwas_metadata.tsv


###------------------ LIFTOVER GWAS SUMSTAT
# Go to the lift_gwas directory and run the script to lift over the GWAS sumstat to hg38 and prepare the final sumstat file
###------------------ LIFTOVER GWAS SUMSTAT


# ### LDREF - Maybe no need to run if i'm running the previouse step -should check
# cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry
# mkdir -p Pantry/pheast/input/LDREF

# for chr in {1..22}; do
# plink2  --vcf phasing_imputation/src/ref_panel/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz \
#         --maf 0.002 \
#         --geno 0.01 \
#         --hwe 1e-10 \
#         --snps-only just-acgt \
#         --vcf-half-call m \
#         --nonfounders \
#         --max-alleles 2 \
#         --set-all-var-ids @_#_'$r'_'$a' \
#         --make-bed \
#         --out Pantry/pheast/input/LDREF/1kg_ldref_hg38_chr${chr}
# done



###-------------------------- Pheast
srun -p cpu-gdml-investor -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A cpu-gdml-investor --pty bash -i
mamba activate pantry_pheast

###------ utricle
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/utricle alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd ../pheast
mkdir -p input/utricle/

cp ../../src/utricle_sample_ids.tsv input/utricle/utricle.individuals.txt
awk '{ print $1 "\tFTL_" $1 }' ../phenotyping/input/utricle/utricle.samples.tsv > input/utricle/utricle.sample_map.txt

scripts/prepare_phenotypes_updated072225.py \
                --indir ../phenotyping/output/utricle \
                --outdir input/utricle/phenotypes/ \
                --map input/utricle/utricle.sample_map.txt \
                --individuals input/utricle/utricle.individuals.txt


# Prep genotypes
plink2  --vcf ../../phasing_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz \
        --keep input/utricle/utricle.individuals.txt \
        --fa ../phenotyping/input/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --ref-from-fa force \
        --maf 0.05 \
        --geno 0.01 \
        --hwe 1e-6 \
        --snps-only just-acgt \
        --nonfounders \
        --max-alleles 2 \
        --chr 1-22 \
        --set-all-var-ids @_#_'$r'_'$a' \
        --make-bed \
        --out input/utricle/utricle_genotype

# Moving the covar files to pheast/input
cp ../../src/utricle_covars.tsv input/utricle

# Now Run the pheast
# Utricle - Dry run
snakemake --configfile config.yml --config dataset=utricle --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=utricle --job-name=pheast_utricle run_pheast_gpu_wrapper.slurm
2700418
2799432
2879353


###------ Saccule
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/saccule alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd ../pheast
mkdir -p input/saccule/

cp ../../src/saccule_sample_ids.tsv input/saccule/saccule.individuals.txt
awk '{ print $1 "\tFTL_" $1 }' ../phenotyping/input/saccule/saccule.samples.tsv > input/saccule/saccule.sample_map.txt

scripts/prepare_phenotypes_updated072225.py \
        --indir ../phenotyping/output/saccule \
        --outdir input/saccule/phenotypes/ \
        --map input/saccule/saccule.sample_map.txt \
        --individuals input/saccule/saccule.individuals.txt

# Preparing the genotypes
plink2  --vcf ../../phasing_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz \
        --keep input/saccule/saccule.individuals.txt \
        --fa ../phenotyping/input/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --ref-from-fa force \
        --maf 0.05 \
        --geno 0.01 \
        --hwe 1e-6 \
        --snps-only just-acgt \
        --nonfounders \
        --max-alleles 2 \
        --chr 1-22 \
        --set-all-var-ids @_#_'$r'_'$a' \
        --make-bed \
        --out  input/saccule/saccule_genotype

# Moving the covar files to pheast/input
cp ../../src/saccule_covars.tsv input/saccule

# Now Run the pheast
# Saccule - Dry run
snakemake --configfile config.yml --config dataset=saccule --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=saccule --job-name=pheast_saccule run_pheast_gpu_wrapper.slurm
2700410
2870214


###------ Vestibular
# Making cross-modality
cd phenotyping
bash scripts/combine_modalities_updated_060425.sh --path_to_output output/vestibular alt_polyA alt_TSS expression RNA_editing isoforms splicing stability

# Copy phenotypes to pheast input/
cd ../pheast
mkdir -p input/vestibular/

cp ../../src/vestibular_sample_ids.tsv input/vestibular/vestibular.individuals.txt
awk '{ print $1 "\tFTL_" $1 }' ../phenotyping/input/vestibular/vestibular.samples.tsv > input/vestibular/vestibular.sample_map.txt

scripts/prepare_phenotypes_updated072225.py \
        --indir ../phenotyping/output/vestibular \
        --outdir input/vestibular/phenotypes/ \
        --map input/vestibular/vestibular.sample_map.txt \
        --individuals input/vestibular/vestibular.individuals.txt


# Preparing the genotypes
plink2  --vcf ../../phasing_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz \
        --keep input/vestibular/vestibular.individuals.txt \
        --fa ../phenotyping/input/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --ref-from-fa force \
        --maf 0.05 \
        --geno 0.01 \
        --hwe 1e-6 \
        --snps-only just-acgt \
        --nonfounders \
        --max-alleles 2 \
        --chr 1-22 \
        --set-all-var-ids @_#_'$r'_'$a' \
        --make-bed \
        --out  input/vestibular/vestibular_genotype

# Moving the covar files to pheast/input
cp ../../src/vestibular_covars.tsv input/vestibular

# Now Run the pheast
# Vestibular - Dry run
snakemake --configfile config.yml --config dataset=vestibular --cores 16 --printshellcmds -n
# Submit on gpu
sbatch --export=ALL,DATASET=vestibular --job-name=pheast_vestibular run_pheast_gpu_wrapper.slurm
2816964
2897577



###------------ ASE
scp username@remote:/path/to/remote/file /path/to/local/destination

bcftools view -v snps phasing_imputation/results/phased_vcf/geno_lifted_phased_hg38_chr22.vcf.gz | grep -v "##" | head

snakemake --cores 16 --printshellcmds -n


GTF_COLLAPSED="matching_analysis/src/Homo_sapiens.GRCh38.110.chr.gtf"
grep -v '^#' $GTF_COLLAPSED | awk '$3=="exon"' | cut -f1,4,5 | gzip -c > matching_analysis/src/GRCh38_exon_regions.tsv.gz


REF="matching_analysis/src/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
plink2  --vcf phasing_imputation/results/phased_vcf/geno_lifted_phased_hg38_allchrs.vcf.gz \
        --rm-dup \
        --maf 0.05 \
        --fa $REF \
        --ref-from-fa force \
        --not-chr X Y \
        --snps-only 'just-acgt' \
        --recode vcf id-paste=iid bgz \
        --max-alleles 2 \
        --out matching_analysis/data/genotype/geno_lifted_phased_hg38_allchrs_inter_cln

#        --output-chr chrM \
tabix -p vcf matching_analysis/data/genotype/geno_lifted_phased_hg38_allchrs_inter_cln.vcf.gz
# Getting a list of chr-ID-POS of SNPs
bcftools query -f '%CHROM\t%ID\t%POS\n' matching_analysis/data/genotype/geno_lifted_phased_hg38_for_mixup_test.vcf.gz | head


bcftools view matching_analysis/data/genotype/geno_lifted_phased_hg38_allchrs_inter_cln.vcf.gz \
                --regions-file matching_analysis/src/GRCh38_exon_regions.tsv.gz \
                -Ou | bcftools view \
                --min-af 0.1:minor \
                -i 'F_MISSING<0.01 && N_ALT=1' \
                -Oz -o matching_analysis/data/genotype/geno_lifted_phased_hg38_for_mixup_test.vcf.gz

# && TYPE="snp"
tabix -p vcf matching_analysis/data/genotype/geno_lifted_phased_hg38_for_mixup_test.vcf.gz

# Number of SNPs in the VCF
bcftools view -v snps matching_analysis/data/genotype/geno_lifted_phased_hg38_for_mixup_test.vcf.gz | grep -v -c "^#" 

# Make .fai and dictionary index of reference genome
REF="matching_analysis/src/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gatk CreateSequenceDictionary \
        R=$REF \
        O=matching_analysis/src/Homo_sapiens.GRCh38.dna_sm.primary_assembly.dict

# .fai index
samtools faidx $REF

### ASEReadCounter - GATK
# Defining inputs and outputs
cd matching_analysis/
REF="src/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
VCF_IN="data/genotype/geno_lifted_phased_hg38_for_mixup_test.vcf.gz"
OUT="results/ase_count"
PATH_TO_BAM="data/bam/"
mkdir -p $OUT

# Getting a list of BAMs
ls -1 data/bam/*.Aligned.sortedByCoord.out.bam > src/bam_list.txt

# Running ASE
for BAM in src/bam_list.txt; do
lib_id=$(echo $BAM | cut -f3 -d"/" | cut -f1 -d".")
gatk ASEReadCounter \
                -R $REF \
                -I $BAM \
                -V $VCF_IN \
                -O ${OUT}/${lib_id}.readcounts.txt \
                --verbosity ERROR
done

snakemake -j 24 --printshellcmds --scheduler greedy -R rna_to_genotype_similarity -n

sbatch run_workflow.slurm
2798510
2799032


###----------------- TWAS
cd Pantry/pheast/twas_assoc/
mamba activate pantry_pheast

snakemake -s run_twas_with_tissue.smk --jobs 2 --printshellcmds --rerun-incomplete -n
sbatch run_twas_assoc.slurm
2799583
2811504
2839409
