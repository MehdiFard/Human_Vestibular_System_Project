### Fixing the fam file IID and making metadata for genotype data 
rm(list = ls())

library(data.table)
library(dplyr)
library(readxl)
library(janitor)
library(lubridate)
library(stringr)
library(tidyr)
library(readr)


fq_list <- fread("src/path_to_fqs.tsv")
View(fq_list)

utricle_fqs <- fq_list %>% 
    filter(tissue == "utricle") %>%
    select(R1, R2)

saccule_fqs <- fq_list %>% 
    filter(tissue == "saccule") %>%
    select(R1, R2)
View(saccule_fqs)

library(fs)
library(purrr)
library(data.table)

### Making symlinks for fastq files
# Function to create a symlink with resolved real path
make_symlink <- function(source_path, dest_dir) {
  # Resolve symlink to real file
  real_path <- normalizePath(source_path, mustWork = TRUE)
  link_name <- file.path(dest_dir, basename(real_path))

  # Create symlink
  file.symlink(from = real_path, to = link_name)
}

# Path for utricle
dest_dir_utr <- "Pantry/phenotyping/input/utricle/fastq/"
# Symlinking utricle fastq files
walk(utricle_fqs$R1, make_symlink, dest_dir = dest_dir_utr)
walk(utricle_fqs$R2, make_symlink, dest_dir = dest_dir_utr)


# Path for saccule
dest_dir_sac <- "Pantry/phenotyping/input/saccule/fastq/"
# Symlinking utricle fastq files
walk(saccule_fqs$R1, make_symlink, dest_dir = dest_dir_sac)
walk(saccule_fqs$R2, make_symlink, dest_dir = dest_dir_sac)


# Dl edit sites ref
download.file("https://github.com/vargasliqin/mpileup/raw/refs/heads/master/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz",
              destfile = "src/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz",
              method = "wget")

# reading edit sites reference file
edit_ref <- fread("src/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz",
              header = FALSE)

# Removing 'chr' prefix to match te reference genome
edit_ref_nochr <- edit_ref %>%
  mutate(V1 = str_remove(V1, "^chr"))

# Writing the modified reference file
write_tsv(edit_ref_nochr, "Pantry/phenotyping/input/ref/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.nochr.bed", 
          col_names = FALSE)


# Making fastq_map file and list of samples
# utricle
fastq_map_utr <- fq_list %>%
  filter(tissue == "utricle") %>%
  select(R1, R2, human_id, read_length) %>%
  mutate(
    R1 = str_replace(R1, "data/rna_fq/(batch\\d+|gold_standard)", "input/utricle/fastq"),
    R2 = str_replace(R2, "data/rna_fq/(batch\\d+|gold_standard)", "input/utricle/fastq")
  )

# Writing the fastq_map file for utricle
write_tsv(fastq_map_utr, "Pantry/phenotyping/input/utricle/utricle.fastq_map.tsv", col_names = F)

# Write sample ids to include in the analysis
write.table(fastq_map_utr$human_id, "Pantry/phenotyping/input/utricle/utricle.samples.tsv", quote = F, col.names = F, row.names = F)

# Saccule
fastq_map_sac <- fq_list %>%
  filter(tissue == "saccule") %>%
  select(R1, R2, human_id, read_length) %>%
  mutate(
    R1 = str_replace(R1, "data/rna_fq/(batch\\d+|gold_standard)", "input/saccule/fastq"),
    R2 = str_replace(R2, "data/rna_fq/(batch\\d+|gold_standard)", "input/saccule/fastq")
  )

# Writing the fastq_map file for utricle
write_tsv(fastq_map_sac, "Pantry/phenotyping/input/saccule/saccule.fastq_map.tsv", col_names = F)

# Write sample ids to include in the analysis
write.table(fastq_map_sac$human_id, "Pantry/phenotyping/input/saccule/saccule.samples.tsv", quote = F, col.names = F, row.names = F)

