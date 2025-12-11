# Making metadata for genotype and RNAseq data and a list of fastq files for mapping

rm(list = ls())

library(data.table)
library(dplyr)
library(readxl)
library(janitor)
library(lubridate)
library(stringr)
library(tidyr)
library(readr)


### Making genotype metadata
# Batch 1 and 2
b12 <- read_excel("data/genotype/genotype_metadata_GDA_human_samples_b1b2b3.xlsx", sheet = 1) %>%
        clean_names()

# Batch 3
b3 <- read_excel("data/genotype/genotype_metadata_GDA_human_samples_b1b2b3.xlsx", sheet = 2) %>%
        clean_names()

# Fix date problem for four samples
b3[which(b3$sample_id == "FTL_31475704"), "age"] <- "12/5/1959"
b3[which(b3$sample_id == "FTL_32517665"), "age"] <- "8/2/1998"
b3[which(b3$sample_id == "FTL_32510348"), "age"] <- "6/10/1962"
b3[which(b3$sample_id == "FTL_32559131"), "age"] <- "5/10/1972"

# Make it as date
b3$age <- as.Date(b3$age, format = "%d/%m/%Y", origin = "1899-12-30")

# Calculate age
birthdate <- b3$age
target_date <- as.Date("2023-01-12") # Had to pick this date to make the age consistent through all batches (there is no birthdate for batch1 and 2)
age <- interval(birthdate, target_date) / years(1)
b3$age <- floor(age)  # round down to full years

# Adding human_id
b3$human_id <- sub("FTL_", "", b3$sample_id)

# Rbind to get the full metadata for all batches
geno_meta <- b12 %>%
    select(sample_id = library_id, human_id, age, sex) %>%
    mutate(human_id = as.character(human_id)) %>%
    bind_rows(b3 %>% select(sample_id, human_id, age, sex)) %>%
    mutate(sex = str_replace_all(sex, c("Male" = "M", "Female" = "F"))) 

# Write out the metadata
write_tsv(geno_meta, "src/metadata_geno.tsv")


### Making RNA metadata
# Making metadata for batch 3
b1_rna <- read_excel("data/rna_fq/rnaseq_metadata_human_b1b2b3_gs.xlsx", sheet = 1) %>%
    clean_names() %>%
    select(library_id = rna_seq_library_id, human_id,
        tissue, age, sex, library_type, blood) %>%
    mutate(tissue = str_replace_all(tissue, c("Saccule" = "saccule", "saccule 1" = "saccule", "saccule 2" = "saccule", "Utricle" = "utricle")),
    batch = "b1",
    read_length = 101) 
# Note: FTL_500 and FTL_501, both are from saccule and from one sample (31661958).


# Making metadata for batch 2
b2_rna <- read_excel("data/rna_fq/rnaseq_metadata_human_b1b2b3_gs.xlsx", sheet = 2) %>%
    clean_names() %>%
    select(library_id = rna_seq_library_id,human_id,tissue,
        age,sex,library_type,blood) %>%
    mutate(tissue = case_when(
        str_detect(tissue, regex("^sac$", ignore_case = TRUE)) ~ "saccule",
        str_detect(tissue, regex("^utricle$", ignore_case = TRUE)) ~ "utricle",
        TRUE ~ tissue),
        batch = "b2",
        read_length = 151)

# Fixing typo in tissue name for library FTL_840
b2_rna[which(b2_rna$tissue == "utrilcle"), "tissue"] <- "utricle"


# Need to rename batch 3 fastq files to make all batch data consistent
rename_b3 <- readxl::read_excel("data/rna_fq/rnaseq_metadata_human_b1b2b3_gs.xlsx", sheet = 3) %>%
    clean_names() %>%
    select(library_id = human_id, sample_id = rna_seq_library_id ) %>%
    mutate(library_id = paste0("FTL_", library_id))

# Getting fastqs paths and renaming them
b3_oldpath_newpath <- list.files("data/rna_fq/batch3", pattern = "FTL_.*\\.fastq\\.gz$", full.names = TRUE) %>%
    tibble(path = .) %>%
    mutate(library_id = str_extract(path, "FTL_\\d+_(Saccule|Utricle)")) %>%
    left_join(rename_b3, by = c("library_id" = "sample_id")) %>%
    mutate(new_path = str_replace(path, library_id, library_id.y))

# Actually rename the files
# file.rename(from = b3_oldpath_newpath$path, to = b3_oldpath_newpath$new_path)

# Making metadata for batch 3
b3_rna <- read_excel("data/rna_fq/rnaseq_metadata_human_b1b2b3_gs.xlsx", sheet = 3) %>%
    clean_names() %>%
    mutate(library_id = paste0("FTL_", human_id),
            human_id = str_extract(rna_seq_library_id, "(?<=FTL_)\\d+")) %>%
    select(library_id, human_id, tissue, age, sex, library_type, blood) %>%
    mutate(tissue = str_replace_all(tissue, c("Saccule" = "saccule", "Utricle" = "utricle")),
    sex = str_replace_all(sex, c("Male" = "M", "Female" = "F")),
    batch = "b3",
    read_length = 151)

# Fix date problem for four samples
b3_rna[which(b3_rna$human_id == "31475704"), "age"] <- "12/5/1959"
b3_rna[which(b3_rna$human_id == "32517665"), "age"] <- "8/2/1998"
b3_rna[which(b3_rna$human_id == "32510348"), "age"] <- "6/10/1962"
b3_rna[which(b3_rna$human_id == "32559131"), "age"] <- "5/10/1972"

# Make it as date
b3_rna$age <- as.Date(b3_rna$age, format = "%d/%m/%Y", origin = "1899-12-30")

# Calculate age
birthdate <- b3_rna$age
target_date <- as.Date("2023-01-12") # Had to pick this date to make the age consistent through all batches (there is no birthdate for batch1 and 2)
age <- interval(birthdate, target_date) / years(1)
b3_rna$age <- floor(age)  # round down to full years

# Making metadata for gold standard
gs_rna <- readxl::read_excel("data/rna_fq/rnaseq_metadata_human_b1b2b3_gs.xlsx", sheet = 4) %>%
    clean_names() %>%
    select(library_id = rna_seq_library_id, human_id,
        tissue, age, , sex, library_type, blood) %>%
    mutate(library_id = str_remove(library_id, "_RNAseq_human"),
    batch = rep("gs", length(human_id)),
    read_length = 101)

# Row bind all batches
meta_rna <- rbind(b1_rna, b2_rna, b3_rna, gs_rna) %>%
    mutate(blood = str_replace_all(blood, c("YES" = "yes", "Yes" = "yes", "NO" = "no")),
    library_type = str_replace_all(library_type, c("Ribodepleted" = "ribodepleted")))


# Merging RNA metadata with genotype metadata
meta_merged <- meta_rna %>%
    left_join(geno_meta, by = "human_id") %>%
    select(library_id, sample_id_geno = sample_id, human_id, tissue,
    age = age.x, sex = sex.y, batch, read_length, blood, library_type) %>%
    mutate(sample_id_geno = replace_na(sample_id_geno, "no_genotype"))


###---------------------------------------------------------------------------- Adding matching results - Beginning of the block.
# Note: comment this block if you don't want to use matching results

# Read in the mix-up table
mix_up <- read.table("matching_analysis/results/rna_to_genotype_summary.tsv", header = T, sep = "\t")

# Extract row for FTL_500_501 and duplicate it
indx_n <- which(mix_up$RNA_ID == "FTL_500_501")
new_rows <- mix_up[indx_n, ] %>%
    slice(rep(1, 2)) %>%
    mutate(RNA_ID = c("FTL_500", "FTL_501"))
# Replace row 46 with the two new rows
mix_up <- bind_rows(mix_up[-indx_n, ], new_rows)

# Merging metadata and mix-up table
mix_meta_for_check <- mix_up %>%
    # select(-top_variant_count) %>%
    inner_join(meta_merged[, c("library_id", "sample_id_geno", "tissue")],
        by = c("RNA_ID" = "library_id"))

# Check the number of matched and unmatched IDs
table(mix_meta_for_check$top_simil_ID == mix_meta_for_check$sample_id_geno)

# Save the results (this one is only for looking at in Excel)
write_tsv(mix_meta_for_check[, c(1,2,8,3,4,6,7,9)], file = "src/mixup_table_forcheck.tsv")

### Now making RNA metadata with matching results
mix_meta_final <- meta_merged %>%
        inner_join(mix_up[, c("RNA_ID", "top_simil_ID", "top_simil")],
            by = c("library_id" = "RNA_ID"))  %>%
        mutate(top_simil_ID = if_else(sample_id_geno == "no_genotype", "no_genotype", top_simil_ID), # If there is no genotype, then top_simil_ID becomes also no_genotype
        matched = if_else(sample_id_geno == top_simil_ID, "Yes", "No")) %>%
        filter(sample_id_geno != "no_genotype") %>%
        mutate(human_id_top_simil = str_extract(top_simil_ID, "(?<=FTL_)\\d+")) %>%
        dplyr::rename(sample_id_geno_org = sample_id_geno,
            human_id_org = human_id,
            sample_id_geno = top_simil_ID,
            human_id = human_id_top_simil)

# I'm overwriting the meta_merged to keep the rest of the script consistent
meta_merged <- mix_meta_final

###---------------------------------------------------------------------------- End of the block

# Write out the metadata
write_tsv(meta_merged, "src/metadata_rna.tsv")

### Subsetting ids and covar files for each tissue + vestibular
# Utricle
samples_utricle <- meta_merged %>%
    filter(tissue == "utricle") %>%
    filter(sample_id_geno != "no_genotype") %>%
    filter(!(matched == "No" & top_simil < 0.81)) %>% # Note: Comment this if you don't want to use matching results
    select(sample_id_geno,age,sex) %>%
    mutate(sex = recode(sex, "M" = 1, "F" = 2)) %>%
    distinct(sample_id_geno, .keep_all = TRUE)

# Write out
write.table(samples_utricle$sample_id_geno, "src/utricle_sample_ids.tsv", quote = F, col.names = F, row.names = F)
write_tsv(samples_utricle, "src/utricle_covars.tsv")

# Saccule
samples_saccule <- meta_merged %>%
    filter(tissue == "saccule") %>%
    filter(sample_id_geno != "no_genotype") %>%
    filter(!(matched == "No" & top_simil < 0.81)) %>% # Note: Comment this if you don't want to use matching results
    select(sample_id_geno,age,sex) %>%
    mutate(sex = recode(sex, "M" = 1, "F" = 2)) %>%
    distinct(sample_id_geno, .keep_all = TRUE)

# Write out
write.table(samples_saccule$sample_id_geno, "src/saccule_sample_ids.tsv", quote = F, col.names = F, row.names = F)
write_tsv(samples_saccule, "src/saccule_covars.tsv")

# Vestibular
samples_vestibular <- meta_merged %>%
    filter(sample_id_geno != "no_genotype") %>%
    filter(!(matched == "No" & top_simil < 0.81)) %>% # Note: Comment this if you don't want to use matching results
    select(sample_id_geno,age,sex) %>%
    mutate(sex = recode(sex, "M" = 1, "F" = 2)) %>%
    distinct(sample_id_geno, .keep_all = TRUE)

# Write out
write.table(samples_vestibular$sample_id_geno, "src/vestibular_sample_ids.tsv", quote = F, col.names = F, row.names = F)
write_tsv(samples_vestibular, "src/vestibular_covars.tsv")

### Making a list of R1 and R2 reads for mapping
# Batch 1
b1_path_fqs <- list.files("data/rna_fq/batch1", pattern = "FTL_.*\\.fastq\\.gz$", full.names = TRUE) %>%
    tibble(path = .) %>%
    mutate(library_id = str_extract(path, "FTL_\\d+"),
        batch = "b1",
        read_length = 101,
        read = str_extract(path, "R[12]")
    ) %>%
    pivot_wider(names_from = read, values_from = path) %>%
    select(R1, R2, library_id, read_length, batch)

# Batch 2
b2_path_fqs <- list.files("data/rna_fq/batch2", pattern = "FTL_.*\\.fastq\\.gz$", full.names = TRUE) %>%
    tibble(path = .) %>%
    mutate(library_id = str_extract(path, "FTL_\\d+"),
        batch = "b2",
        read_length = 151,
        read = str_extract(path, "R[12]")
    ) %>%
    pivot_wider(names_from = read, values_from = path) %>%
    select(R1, R2, library_id, read_length, batch)

# Batch 3
b3_path_fqs <- list.files("data/rna_fq/batch3", pattern = "FTL_.*\\.fastq\\.gz$", full.names = TRUE) %>%
    tibble(path = .) %>%
    mutate(library_id = str_extract(path, "FTL_\\d+"),
        batch = "b3",
        read_length = 151,
        read = str_extract(path, "R[12]")
    ) %>%
    pivot_wider(names_from = read, values_from = path) %>%
    select(R1, R2, library_id, read_length, batch)

# Batch gold standard
gs_path_fqs <- list.files("data/rna_fq/gold_standard", pattern = "FTL_.*\\.fastq\\.gz$", full.names = TRUE) %>%
    tibble(path = .) %>%
    mutate(library_id = str_extract(path, "FTL_\\d+"),
        batch = "gs",
        read_length = 101,
        read = str_extract(path, "R[12]")
    ) %>%
    pivot_wider(names_from = read, values_from = path) %>%
    select(R1, R2, library_id, read_length, batch)

# Row bind all batches
all_path_fqs <- rbind(b1_path_fqs,
                    b2_path_fqs,
                    b3_path_fqs,
                    gs_path_fqs)

# Add human_id
rna_meta_with_path <- meta_merged %>%
    select(library_id,human_id, tissue) %>%
    left_join(all_path_fqs, by = "library_id") %>%
    select(R1, R2, read_length, library_id, human_id, batch, tissue)

# Write out the fastq files list
write_tsv(rna_meta_with_path, "src/path_to_fqs.tsv")

# # Function to get read length from gzipped FASTQ (first read only)
# get_read_length <- function(fq) {
#     con <- gzfile(fq, "rt")
#     readLines(con, n = 2) -> lines
#     close(con)
#     nchar(lines[2])
# }
# # mutate(read_length = get_read_length(R1))


### END