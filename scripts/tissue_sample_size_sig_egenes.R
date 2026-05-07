library(data.table)

pheast_dir <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/pheast"
phenos_dir <- "/data/hps/assoc/private/gdml/user/dmunro/pantry/GTEx/pheast"

tissues <- basename(list.dirs(pheast_dir, recursive = FALSE, full.names = TRUE))

tissue_summary <- rbindlist(lapply(tissues, function(tissue) {
  qtl_file <- file.path(pheast_dir, tissue, "output", "qtl", "expression.cis_independent_qtl.txt.gz")
  pheno_file <- file.path(phenos_dir, tissue, "input", "phenotypes", "expression.bed.gz")

  sig_egenes <- NA_integer_
  sample_size <- NA_integer_

  if (file.exists(qtl_file)) {
    qtl <- fread(cmd = sprintf("zcat %s", shQuote(qtl_file)), select = "phenotype_id")
    sig_egenes <- uniqueN(qtl$phenotype_id)
  }

  if (file.exists(pheno_file)) {
    phen <- fread(cmd = sprintf("zcat %s | head -n 1", shQuote(pheno_file)), nrows = 0)
    sample_size <- ncol(phen) - 4L
  }

  data.table(
    tissue = tissue,
    sample_size = sample_size,
    sig_egenes = sig_egenes
  )
}), fill = TRUE)

setorder(tissue_summary, tissue)

out_file <- "tmp/tissue_sample_size_sig_egenes_rerun.tsv"
fwrite(tissue_summary, out_file, sep = "\t")

cat("Saved:", out_file, "\n")
print(head(tissue_summary))
