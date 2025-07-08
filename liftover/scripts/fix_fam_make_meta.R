### Fixing the fam IID and making genotype metadata
rm(list = ls())

library(data.table)
library(dplyr)
library(readxl)
library(janitor)
library(lubridate)
library(stringr)

# Fix problem with IID and FID
fam_cln <- read.table("data/genotype/geno_merged_hg19.fam") %>%
        mutate(V2 = sub("^FTL[ _]*([0-9]+)$", "FTL_\\1", V2),
                V1 = seq(1, length(V1)))
write.table(fam_cln, "data/genotype/geno_merged_hg19_fixed.fam", col.names = F, row.names = F, quote = F, sep = "\t")


## Making Metadata
# Batch 1 and 2
b12 <- readxl::read_excel("../data/genotype/genotype_metadata_GDA_human_samples_b1b2b3.xlsx", sheet = 1)
b12 <- janitor::clean_names(b12)

# Batch 3
b3 <- readxl::read_excel("../data/genotype/genotype_metadata_GDA_human_samples_b1b2b3.xlsx", sheet = 2)
b3 <- janitor::clean_names(b3)

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
    mutate(sex = str_replace_all(sex, c("Male" = "M", "Female" = "F"))) %>%
    slice(match(fam_cln$V2, sample_id))

# Write out the metadata
write.table(geno_meta, "results/geno_metadata.tsv", col.names = T, row.names = F, quote = F, sep = "\t")


### END