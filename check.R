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
library(fs)
library(purrr)

# set wd
setwd("/data/hps/assoc/private/gdml/user/sesmae/projects/vestibulart_project_human_pantry")

fq_list <- fread("src/path_to_fqs.tsv")
View(fq_list)

# Utricle
utricle_fqs <- fq_list %>% 
    filter(tissue == "utricle") %>%
    select(R1, R2)

# Saccule
saccule_fqs <- fq_list %>% 
    filter(tissue == "saccule") %>%
    select(R1, R2)

# Both tissues as vesibular system
vestibular_fqs <- fq_list %>%
    select(R1, R2)


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

# Path for vestibular system
dest_dir_vest <- "Pantry/phenotyping/input/vestibular/fastq"
dir.create(dest_dir_vest, recursive = TRUE)
# Symlinking vestibular fastq files
walk(vestibular_fqs$R1, make_symlink, dest_dir = dest_dir_vest)
walk(vestibular_fqs$R2, make_symlink, dest_dir = dest_dir_vest)


### Making fastq_map file and list of samples
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
write.table(unique(fastq_map_utr$human_id), "Pantry/phenotyping/input/utricle/utricle.samples.tsv", quote = F, col.names = F, row.names = F)


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
write.table(unique(fastq_map_sac$human_id), "Pantry/phenotyping/input/saccule/saccule.samples.tsv", quote = F, col.names = F, row.names = F)


# vestibular
fastq_map_vestibular <- fq_list %>%
  select(R1, R2, human_id, read_length) %>%
  mutate(
    R1 = str_replace(R1, "data/rna_fq/(batch\\d+|gold_standard)", "input/vestibular/fastq"),
    R2 = str_replace(R2, "data/rna_fq/(batch\\d+|gold_standard)", "input/vestibular/fastq")
  )

# Writing the fastq_map file for vestibular
write_tsv(fastq_map_vestibular, "Pantry/phenotyping/input/vestibular/vestibular.fastq_map.tsv", col_names = F)
# Write sample ids to include in the analysis
write.table(unique(fastq_map_vestibular$human_id), "Pantry/phenotyping/input/vestibular/vestibular.samples.tsv", quote = F, col.names = F, row.names = F)


### Dl edit sites ref
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


### END
