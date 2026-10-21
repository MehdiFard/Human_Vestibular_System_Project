### Preparing GWAS sumstats for TWAS analysis

getwd()
# run from analyses/06_twas/gtex_weights
rm(list = ls())

library(data.table)
library(dplyr)


# Path to GWAS sumstats
path_to_sumstats <- "../../04_gwas_sumstats/results/twas_input/"

# Getting list of GWAS sumstats files
sumstat_fn <- list.files(path_to_sumstats, pattern = "sumstats.gz")

# Making input and output file paths
input_files <- paste0(path_to_sumstats, sumstat_fn)
output_files <- paste0("gwas_sumstats/", sub(".sumstats.gz", "", sumstat_fn), ".chr_b38.sumstats.gz")

# Loop through each file, modifying SNP IDs and writing output
for (i in seq_along(input_files)) {
    print(paste("Processing file:", input_files[i]))

    sumstat <- fread(input_files[i])[
    , SNP := paste0("chr", SNP, "_b38")
]
    fwrite(sumstat, output_files[i], sep = "\t")
    
    print(paste("Written output to:", output_files[i]))
}


### Preparing LD reference .bim file
# Listing LD reference .bim files
ld_ref_fn <- list.files("../../input/LDREF/", pattern = ".bim")

# Reading, modifying, and writing the .bim file
for(i in seq_along(ld_ref_fn)){
    print(paste("Processing LD reference file:", ld_ref_fn[i]))
    bim_file <- fread(paste0("../../input/LDREF/", ld_ref_fn[i])) %>%
        mutate(
            V2 = paste0("chr", V2, "_b38"),
            V1 = paste0("chr", V1)
        )
    # Writing modified .bim file
    fwrite(bim_file, paste0("ld_ref/", ld_ref_fn[i]), col.names = F, sep = "\t")
}


# Copying .bed and .fam files
system("cp ../../input/LDREF/1kg_ldref_hg38_chr{1..22}.{bed,fam} ld_ref/")


### END