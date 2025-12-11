#!/bin/bash
# Script to lift over GWAS summary statistics from hg19 to hg38
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast/lift_gwas

# Activate the environment
mamba activate Env1


# Copy GWAS results
cp /data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/gwas_results/balance_meta1.tbl.sorted.filtered.gz gwas_sumstats/clifford_raw_sumstats/

# Make bed file from gwas sumstat
zcat balance_meta1.tbl.sorted.filtered.gz \
        | awk 'NR>1 {print "chr"$1"\t"$2-1"\t"$2"\t"$3}' \
        > balance_meta1.hg19.bed

# Dl chain file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Dl liftOver tool 
wget gwas_sumstats_prep/src/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x gwas_sumstats_prep/src/liftOver

# Run liftOver to lift gwas sumstat from hg19 to hg38
./gwas_sumstats_prep/src/liftOver balance_meta1.hg19.bed hg19ToHg38.over.chain.gz balance_meta1.hg38.bed unMapped.bed

for chr in {1..22}; do
        echo "Processing chr$chr"
        grep -w "chr${chr}" balance_meta1.hg38.bed > balance_meta1.hg38.chr${chr}.bed

        # Extracting ALT and REF from the reference panel
        REF="../../../phasing_imputation/src/ref_panel/ALL.chr${chr}_GRCh38.genotypes.20170504.bcf"
        bcftools view -R balance_meta1.hg38.chr${chr}.bed $REF -Oz -o balance_meta1.hg38.chr${chr}_with.ref.alt.vcf.gz
        bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\n' balance_meta1.hg38.chr${chr}_with.ref.alt.vcf.gz >> balance_meta1.hg38.ref.alt.bed

        # Remove the temporary files
        rm balance_meta1.hg38.chr${chr}.bed
        rm balance_meta1.hg38.chr${chr}_with.ref.alt.vcf.gz
done

# Updating gwas results positions to hg38
Rscript update_gwas_pos.R \
        balance_meta1.hg38.ref.alt.bed \
        ../input/gwas_sumstat/balance_meta1.tbl.sorted.filtered.gz \
        balance_meta1.tbl.sorted.filtered_hg38.tsv


# Making the final sumstat file
mamba activate ldsc
LDSC="/data/hps/assoc/private/gdml/user/sesmae/tools/ldsc/munge_sumstats.py" # to check run $LDSC -h

# Running the LDSC
$LDSC --sumstats balance_meta1.tbl.sorted.filtered_hg38.tsv \
                        --out balance_meta1.hg38 \
                        --N 420000

# Copy the results to the input directory
cp balance_meta1.hg38.sumstats.gz ../input/gwas_sumstat/

# Unzip the sumstats file
gunzip ../input/gwas_sumstat/balance_meta1.hg38.sumstats.gz

# Making the gwas metadata file
awk 'BEGIN { print "trait\tsample_size\nbalance_meta1.hg38\t420000" }' > ../input/gwas_sumstat/gwas_metadata.tsv


### END

# Make bed file from gwas sumstat
zcat gwas_sumstats/clifford_raw_sumstats/balance_meta_EUR.gz \
        | awk 'NR>1 {print "chr"$1"\t"$2-1"\t"$2"\t"$3}' \
        > Pantry/pheast/lift_gwas/gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg19.bed


./liftOver gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg19.bed hg19ToHg38.over.chain gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.bed gwas_sumstats/formatted_sumstats/unMapped.bed

for chr in {1..22}; do
        echo "Processing chr$chr"
        grep -w "chr${chr}" gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.bed > gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}.bed

        # Extracting ALT and REF from the reference panel
        REF="../../../phasing_imputation/src/ref_panel/ALL.chr${chr}_GRCh38.genotypes.20170504.bcf"
        bcftools view -R gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}.bed $REF -Oz -o gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}_with.ref.alt.vcf.gz
        bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT\n' gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}_with.ref.alt.vcf.gz >> gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.ref.alt.bed

        # Remove the temporary files
        rm gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}.bed
        rm gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.chr${chr}_with.ref.alt.vcf.gz
done


Rscript update_gwas_pos.R \
        gwas_sumstats/formatted_sumstats/balance_meta_EUR.hg38.ref.alt.bed \
        gwas_sumstats/clifford_raw_sumstats/balance_meta_EUR.gz \
        gwas_sumstats/formatted_sumstats/balance_meta_EUR.tbl.sorted.filtered_hg38.tsv

LDSC="/data/hps/assoc/private/gdml/user/sesmae/tools/ldsc/munge_sumstats.py" # to check run $LDSC -h
$LDSC --sumstats  gwas_sumstats/formatted_sumstats/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_hg38_ClnSnpId.tsv.gz \
                        --out gwas_sumstats/formatted_sumstats/GP2_ALL_EUR_ALL_DATASET_HG38_12162024_hg38 \
                        --maf-min 0 \
                        --N 10000


mamba run -n ldsc python "$LDSC" --sumstats formatted_sumstats/gcst90027158_formatted_hg38_ChrPosId.tsv.gz \
                                --out gcst90027158_formatted_hg38_ChrPosId \
                                --maf-min 0 \
                                --N 100000 \
                                --snp SNP \
                                --ignore rsID,HETPVALVARIANT_ALTERNATE_ID,OR


mamba run -n ldsc python "$LDSC" --sumstats formatted_sumstats/tissinketal_commbio2022_cerebralvolume_formatted_hg38_ChrPosId.tsv.gz \
                                --out cerebral_volume_formatted_hg38_ChrPosId \
                                --maf-min 0 \
                                --N 100000 \
                                --snp SNP \
                                --ignore  rsID,HETPVALVARIANT_ALTERNATE_ID,OR,ID,L95,U95,SIGNED_SUMSTAT,UNIQUE_ID \
                                --signed-sumstats BETA,1




wget -P gwas_sumstats/clifford_raw_sumstats/ https://api.kpndataregistry.org/api/d/7j5797
unzip gwas_sumstats/clifford_raw_sumstats/7j5797 -d gwas_sumstats/clifford_raw_sumstats/
rm -r gwas_sumstats/clifford_raw_sumstats/7j5797 gwas_sumstats/clifford_raw_sumstats/__MACOSX

zcat gwas_sumstats/clifford_raw_sumstats/GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz

#!/usr/bin/env Rscript

snakemake -j 24 --printshellcmds --scheduler greedy -n


sbatch run_job.slurm
4570412
4570436
4570519
4570579


sbatch run_job_cluster.slurm
4570540
4570544
4570551
4570563
4570568

trait="XYZ"
if [[ "$trait" =~ ^(cerebellar_volume|cerebral_volume|subcortical_volume)$ ]]; then
        echo "Good"
else
        echo "Not good"
fi

/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/gwas_sumstats_prep

/data/hps/home/sesmae/.config/snakemake/slurm/config.yaml

snakemake -s Snakefile --profile slurm --printshellcmds --scheduler greedy -n
