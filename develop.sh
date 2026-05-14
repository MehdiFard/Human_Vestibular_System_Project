#!/bin/bash

# $ sacctmgr show assoc user=$USER format=User,Account,Partition
#       User    Account  Partition 
# ---------- ---------- ---------- 
#     sesmae       core         cpu-core 
#     sesmae       core         cpu-test 
#     sesmae       core         gpu-core 
#     sesmae       core         gpu-test 
#     sesmae       cpu-gdml-+   cpu-gdml-+ 
#     sesmae       gdml         cpu-core 
#     sesmae       gdml         cpu-gdml 
#     sesmae       gdml         cpu-pwb 
#     sesmae       gdml         cpu-pwb_s+ 
#     sesmae       gdml         gpu-core 
#     sesmae       gdml         gpu-gdml 
#     sesmae       gdml         gpu-pwb 
#     sesmae       gpu-gdml-+   gpu-gdml-+ 

# sesmae@login-2.hpc.childrens.sea.kids
# srun -p cpu-gdml-investor -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A cpu-gdml-investor --pty bash -i
srun -p cpu-gdml -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A gdml --pty bash -i
srun -p cpu-core -N 1 -n 1 -c 16 --mem=64g --time=8:00:00 -A core --pty bash -i
mamba activate Env1


salloc --partition=cpu-core-sponsored --cpus-per-task=8 --mem=32G --time=2:00:00
srun --pty bash


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


snakemake assemble_expression_bed   --configfile config.yml   --config dataset=saccule   --jobs 16   --default-resources mem_mb=64000 cpus=16   --printshellcmds   --rerun-incomplete   --latency-wait 60 -n
snakemake assemble_expression_bed   --configfile config.yml   --config dataset=utricle   --jobs 16   --default-resources mem_mb=64000 cpus=16   --printshellcmds   --rerun-incomplete   --latency-wait 60 -n
snakemake assemble_expression_bed   --configfile config.yml   --config dataset=vestibular   --jobs 16   --default-resources mem_mb=64000 cpus=16   --printshellcmds   --rerun-incomplete   --latency-wait 60 -n


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
# Go to the ./Pantry/pheast/lift_gwas directory and run the script to lift over the GWAS sumstat to hg38 and prepare the final sumstat file
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
scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/pca_expression_count_tissue.png .
scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/pca_expression_count_batch.png .

# From local to remote
rsync -avz --partial --progress \
        -e "ssh -o ServerAliveInterval=30 -o ServerAliveCountMax=20" \
        ./balance_meta_all_151.tab \
        sesmae@login-2.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/gwas_sumstats_prep/just_for_liftover/just_for_test/


rsync -avz --partial --progress \
        /Users/sesmae/Downloads/OneDrive_1_3-9-2026/*.gz \
        sesmae@login-2.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/data/gwas_raw_sumstats/Royce_files/


# From remote to local - run on local
FILE="/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/paper/Supp_Fig3_panel_*.png"
rsync -avz --progress \
        sesmae@login-2.hpc.childrens.sea.kids:$FILE .
        
rsync -avz --progress \
        sesmae@login-2.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/paper/figure2e.pdf .


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
4625309



seff 4623743
sacct -j 4623743 --format=JobID,JobName%30,State,Elapsed,ReqMem,MaxRSS,AveRSS,ExitCode


### Path to phased ref panel
/data/hps/assoc/private/gdml/user/twil15/1kGP_pipeline/resources/1kG/hg38


### Downloading GTEx color palette for plotting
wget https://raw.githubusercontent.com/broadinstitute/gtex-viz/master/boxplot/dev/colors.json -O src/gtex_colors.json

scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/pca_expression_gtex_vestibular_tissue.png .
scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/pca_expression_tmm_logcpm_vestibular.png .
scp sesmae@login-2.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/VolcanoPlot_utricle_vs_saccule_batch_fitted.jpeg .


rsync -avz "OneDrive_2_11-5-2025.zip" sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast/lift_gwas/gwas_sumstats/
rsync -avz "OneDrive_3_11-5-2025.zip" sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast/lift_gwas/gwas_sumstats/


### RUN DGE
sbatch Pantry/phenotyping/run_dge.slurm
4737203

# runing differential_gene_expression.R
sbatch Pantry/phenotyping/run_dge.slurm
4737409
4813242

###-------------------------------------------------------- Runing pheast after removing outliers
# Removing this 3 samples: "31947386", "31907773", "32073088"
# Genotype
cd Pantry/pheast_no_outlier/
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


scripts/prepare_phenotypes_updated072225.py \
        --indir ../phenotyping/output/vestibular \
        --outdir input/vestibular/phenotypes/ \
        --map input/vestibular/vestibular.sample_map.txt \
        --individuals input/vestibular/vestibular.individuals.txt


## Now Run the pheast
mamba activate pantry_pheast

#  Dry run
snakemake --configfile config.yml --config dataset=vestibular --cores 24 --printshellcmds -n
snakemake --configfile config.yml --config dataset=vestibular -n

# Submit on gpu
sbatch --export=ALL,DATASET=vestibular --job-name=pheast_vestibular run_pheast_gpu_wrapper.slurm
2816964
2897577


snakemake --profile slurm \
                --config dataset=vestibular \
                --configfile config.yml \
                --printshellcmds \
                --keep-going \
                --rerun-incomplete \
                -j 20 -n

sbatch run_snakemake_profile.slurm
4945407
4950958
6210584

# ~/.config/snakemake/slurm/config.yaml

# executor: slurm
# default-resources:
#   runtime: "24h"
#   mem_mb: 32000
#   slurm_partition: "cpu-core-sponsored"
#   slurm_account: "cpu-core-sponsored"
# set-resources:
#   tensorqtl_cis:
#     slurm_partition: "gpu-core-sponsored"
#     slurm_account: "gpu-core-sponsored"
#     slurm_extra: "'--gres=gpu:1'"
#   tensorqtl_cis_independent:
#     slurm_partition: "gpu-core-sponsored"
#     slurm_account: "gpu-core-sponsored"
#     slurm_extra: "'--gres=gpu:1'"
#   tensorqtl_trans:
#     slurm_partition: "gpu-core-sponsored"
#     slurm_account: "gpu-core-sponsored"
#     slurm_extra: "'--gres=gpu:1'"
# use-conda: true
# latency-wait: 60
# max-jobs-per-timespan: "1/1s"
# max-status-checks-per-second: 5

sacct -j 2816964 --format=JobID,State,Partition,AllocTRES,ReqTRES,AllocCPUS,NodeList


###---------------------- TWAS after removing outliers
cd Pantry/pheast_no_outlier/twas
mamba install -c conda-forge fastparquet

# Dry run
snakemake -s run_twas_with_tissue.smk --jobs 2 --printshellcmds --rerun-incomplete -n

# Submit the job
sbatch run_twas_assoc.slurm
5004182
6210288



/data/hps/assoc/private/gdml/user/dmunro/pantry/pheast/{tissue}/intermediate/twas/

/data/hps/assoc/private/gdml/user/dmunro/pantry/pheast/PRSTTE/intermediate/twas/expression/ENSG00000000460.wgt.RDat
load("/data/hps/assoc/private/gdml/user/dmunro/pantry/pheast/PRSTTE/intermediate/twas/expression/ENSG00000000460.wgt.RDat")

# SNP format in GTEx weight file:
chr1_169162312_A_C_b38 - rs12755385


for tissue in $(ls -1 /data/hps/assoc/private/gdml/user/dmunro/pantry/pheast/); do
        echo "Linking $tissue"
        mkdir -p weights/$tissue/
        ln -s /data/hps/assoc/private/gdml/user/dmunro/pantry/pheast/$tissue/intermediate/twas weights/$tissue
done

# rsync -ah  --info=progress2


###---------------------- TWAS with GTEx tissue weights
cd Pantry/pheast_no_outlier/twas/twas_gtex
mamba activate pantry_pheast

# Dry run
snakemake -s run_twas_gtex_tissue.smk --jobs 10 --printshellcmds --rerun-incomplete -n
snakemake -s run_twas_gtex_tissue.smk -j 1 -p \
        results/BRNAMY/intermediate/expression/anx_eur/expression.anx_eur.chr5.tsv -n


# Dry run - I'm using Slurm profile (loctaed here ~/.config/snakemake/slurm/config.yaml)
snakemake -s run_twas_gtex_tissue.smk \
                --profile slurm \
                --jobs 250 \
                --keep-going \
                --rerun-incomplete \
                -n

# Submit the job
sbatch run_twas_gtex_tissues.slurm
5025548




tmux new -s twas_gtex
Ctrl+b d    → detach
tmux a      → reattach
tmux attach -t twas_gtex

# Run Snakemake
snakemake -s run_twas_gtex_tissue.smk \
                --profile slurm \
                --jobs 250 \
                --keep-going \
                --rerun-incomplete


scontrol show partitions


git config --global user.name "Mehdi Esmaeilifard"
git config --global user.email "mehdi.esmaeilifard@seattlechildrens.org"
git config --global core.editor "code --wait"
git config --global core.autocrlf input
git config --global -e



grep -R --exclude="twas_hits.tsv" -n \
        "ENSG00000167232__ENST00000300619" \
        Pantry/pheast_no_outlier/twas/twas_gtex/results/*/output/ > tmp/search_res.txt

git status -s
git add . # stage files 
git restore --staged . # unstage
git log --oneline # see priviuos commits
git commit -m "Add additional analysis scripts"

# Adding the GitHub remote to your local repo
git remote add origin https://github.com/MehdiFard/Human_Vestibular_System_Project.git

# Check
git remote -v

# Find branch name
git branch -M main # GitHub using main, so let's change it to main
git branch

# Making a key
ssh-keygen -t ed25519 -C "mehdi.esmaeilifard@gmail.com"
git remote set-url origin git@github.com:MehdiFard/Human_Vestibular_System_Project.git

# Pull first, cause I have License and Readme.md file in github but not here
git pull --rebase origin main

# Push
git push -u origin main


###------ Getting some stats
## Lift over
# in bim file after merge, hg19, 871,077 (maf > 0.01)

# No of rejected variants
bcftools view -H liftover/results/rejected_variants_to_lift.vcf | wc -l # 394,795 / for only SNPs add -v snps

# Why rejected
bcftools view -H liftover/results/rejected_variants_to_lift.vcf | cut -f7 | sort | uniq -c
        394342 MismatchedRefAllele
        453 NoTarget


# No of lifted variants
bcftools view -H liftover/results/geno_lifted_hg38.vcf.gz | wc -l # 476282


## Phasing and imputation
vcf="/data/hps/assoc/private/gdml/user/sesmae/refs_gtfs_indx_vcfs/vcfs/1kgp_phase3_hg38_phased/ALL.chr21_GRCh38.genotypes.20170504.vcf.gz"
# To check genome built
bcftools view -H $vcf | head | cut -f1-6

# To get saample siz
bcftools query -l $vcf | wc -l

# To get number of snps after imputation
bcftools view -H phasing_imputation/results/imputed_vcf/geno_lifted_phased_imputed_hg38_allchrs.vcf.gz | wc -l



### R code to investigate vertigo
###--------------------------------------------------------
sumstats_path = "data/gwas_raw_sumstats/vertigo_GRCH38_sorted.tsv.gz"
build = "GRCh38"
out_file = "formatted_sumstats/vertigo_formatted_hg38.tsv.gz"

df = fread("data/gwas_raw_sumstats/vertigo_GRCH38_sorted.tsv.gz")
dim(df)
head(df)
colnames(df)[3] <- "BP"

# 1) Effect size on log scale
df[, BETA := log(OR)]

# 2) Convert SE on OR scale -> SE on log(OR) scale
df[, SE_beta := SE / OR]

# 3) Z-score on log scale
df[, Z := BETA / SE_beta]

# 4) Two-sided P-value
df[, P := 2 * pnorm(-abs(Z))]
range(df$P)

# 5) Clip P to (0,1) for qqman safety
df[, P := pmin(pmax(P, 1e-300), 1 - 1e-15)]


df[, BETA := log(OR)]
df[, Z := BETA / (log(SE))]
df[, P := 2 * pnorm(-abs(Z))]
df[, P := 2 * pnorm(abs(Z), lower.tail = FALSE)]

range(df$P)

table(p_value < 1e-8)

df[, q := p.adjust(P, method = "BH")]
table(df$q < 0.05)

df[, p_bonf := pmin(P * .N, 1)]
table(df$p_bonf < 0.05)

fwrite(df, "data/gwas_raw_sumstats/vertigo_GRCH38_sorted_cln.tsv.gz", sep = "\t")

install.packages("qqman")
library(qqman)
head(df2)
df2 <- df[, .(CHR = CHR, BP = BP, SNP = SNP, P = P)]

png("figures/Manhattan_vertigo.png", width = 1000, height = 500, res = 100)
manhattan(df2, chr = "CHR", bp = "BP", p = "P", snp = "SNP", genomewideline = -log10(5e-08))
dev.off()

range(df2$P)
summary(df2$P)
sum(is.na(df2$P))
sum(!is.finite(df2$P))
sum(df2$P == 0, na.rm = TRUE)

# Replace exact 0 and 1 with safe bounds
df2[, P := pmin(pmax(P, 1e-300), 1 - 1e-15)]
###--------------------------------------------------------



scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/paper/* .

###---------------------- Fine-mapping with SuSiE
mamba create -n susie_env -c conda-forge -c bioconda r-susier r-dplyr r-tidyverse r-data.table r-readr r-tidyr r-stringr r-devtools r-R.utils r-BEDMatrix r-plink2R r-SNPRelate r-bigsnpr
mamba install -c conda-forge r-dplyr r-tidyverse r-data.table r-readr r-tidyr r-stringr r-devtools r-R.utils r-BEDMatrix r-plink2R r-SNPRelate r-bigsnpr
mamba install conda-forge::r-bigsnpr

mamba create -n susie_env2 \
        -c conda-forge -c bioconda \
        r-base \
        r-susier \
        r-data.table \
        r-tidyverse \
        r-R.utils \
        r-BEDMatrix \
        r-plink2R \
        r-SNPRelate \
        r-bigsnpr \
        r-devtools \
        r-matrix \
        r-mvtnorm \
        r-irlba \
        r-rcpp \
        r-rcpparmadillo \
        r-GenomicRanges \
        r-SummarizedExperiment \
        r-VariantAnnotation \
        r-qvalue \
        r-coloc \
        r-ashr \
        r-ggplot2 \
        r-patchwork \
        r-cowplot \
        r-arrow \
        r-fst


mamba create -n susie_env2 \
        -c conda-forge -c bioconda \
        r-base \
        r-susier \
        r-data.table \
        r-tidyverse \
        r-r.utils \
        r-bedmatrix \
        r-bigsnpr \
        r-devtools \
        r-matrix \
        r-mvtnorm \
        r-irlba \
        r-rcpp \
        r-rcpparmadillo \
        bioconductor-genomicranges \
        bioconductor-summarizedexperiment \
        bioconductor-variantannotation \
        bioconductor-qvalue \
        bioconductor-snprelate \
        r-coloc \
        r-ashr \
        r-ggplot2 \
        r-patchwork \
        r-cowplot \
        r-arrow \
        r-fst
mamba install -n susie_env2 -c conda-forge r-bigsnpr


scp sesmae@login-1.hpc.childrens.sea.kids:/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/figures/pca_expression_tmm_logcpm_sac_utr_* .



mkdir -p src/ucsc_recomb/hg38/
cd src/ucsc_recomb/hg38/
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recomb1000GAvg.bw



grep -H "ENSG00000167232__ENST00000300619" \
Pantry/pheast_no_outlier/twas/twas_gtex/results/*/output/twas_hits.tsv \
| sed -E 's|.*/results/([^/]+)/output/twas_hits.tsv:|\1\t|' > src/znf91_hits_in_gtex2.tsv



for f in Pantry/pheast_no_outlier/twas/twas_gtex/results/*/output/twas_hits.tsv; do
        tissue=$(basename "$(dirname "$(dirname "$f")")")
        if grep -q "ENSG00000167232__ENST00000300619" "$f"; then
                grep "ENSG00000167232__ENST00000300619" "$f" | sed "s/^/${tissue}\t/"
        else
                echo -e "${tissue}\tNO_HIT"
        fi
done > src/znf91_hits_in_gtex.tsv


# This one is not much good
grep -R --exclude="twas_hits.tsv" -n \
        "ENSG00000167232__ENST00000300619" \
        Pantry/pheast_no_outlier/twas/twas_gtex/results/*/output/ > tmp/search_res.txt



### Install claude 

npm install -g @anthropic-ai/claude-code
curl -fsSL https://claude.ai/install.sh | bash

curl -fsSL https://claude.ai/install.sh | bash
wget -qO- https://claude.ai/install.sh | bash

conda install -c conda-forge nodejs
node --version
npm --version
npm install -g @anthropic-ai/claude-code

which jq
ldd $(which jq)