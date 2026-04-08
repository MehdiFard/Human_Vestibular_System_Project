library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)
library(tomics)

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]
ref_col <- snakemake@config[["col_names"]][["other_allele"]]
alt_col <- snakemake@config[["col_names"]][["effect_allele"]]
beta_col <- snakemake@config[["col_names"]][["beta"]]
se_col <- snakemake@config[["col_names"]][["se"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]

if(!is.null(snakemake@params$p_col)) {
  p_col <- snakemake@params$p_col
} else {
  p_col <- snakemake@config[["col_names"]][["p_value"]]
}

distance_window <- snakemake@params[['distance_window']]
index_threshold <- snakemake@params[['index_threshold']]

if(p_col != snakemake@config[["col_names"]][["p_value"]]) {
  dat <- fread(snakemake@input[[1]], sep = '\t', select = c(chr_col, bp_col, ref_col, alt_col, rsid_col, p_col, snakemake@config[["col_names"]][["p_value"]], beta_col, se_col))
} else {
  dat <- fread(snakemake@input[[1]], sep = '\t', select = c(chr_col, bp_col, ref_col, alt_col, rsid_col, p_col, beta_col, se_col))
}

if(!snakemake@params$mhc) {
  dat <- dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

dat[, `:=` (chr = as.character(chr), pos = as.integer(pos), p = as.numeric(p)), env = list(chr = chr_col, pos = bp_col, p = p_col)]

dat <- distance_clump(dat, chr_col = chr_col, bp_col = bp_col, p_col = p_col, distance_window = distance_window, index_threshold = index_threshold)

dat <- dat[order(get(chr_col), get(bp_col))]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
