#!/bin/bash
# Script to lift over GWAS summary statistics from hg19 to hg38
cd /data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry/Pantry/pheast/lift_gwas

# Activate the environment
mamba activate Env1

# Copy GWAS results
cp /data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/gwas_results/balance_meta1.tbl.sorted.filtered.gz .

# Make bed file from gwas sumstat
zcat balance_meta1.tbl.sorted.filtered.gz \
        | awk 'NR>1 {print "chr"$1"\t"$2-1"\t"$2"\t"$3}' \
        > balance_meta1.hg19.bed

# Dl chain file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Dl liftOver tool 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

# Run liftOver to lift gwas sumstat from hg19 to hg38
./liftOver balance_meta1.hg19.bed hg19ToHg38.over.chain.gz balance_meta1.hg38.bed unMapped.bed

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
