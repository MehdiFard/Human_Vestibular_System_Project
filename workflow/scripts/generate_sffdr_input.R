library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config[["col_names"]][["chromosome"]]
bp_col <- snakemake@config[["col_names"]][["base_pair_location"]]
ref_col <- snakemake@config[["col_names"]][["other_allele"]]
alt_col <- snakemake@config[["col_names"]][["effect_allele"]]
## beta_col <- snakemake@config[["col_names"]][["beta"]]
## se_col <- snakemake@config[["col_names"]][["se"]]
p_col <- snakemake@config[["col_names"]][["p_value"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)
stat_cols <- c(p_col)

cols_to_keep <- c(coord_cols, stat_cols)

primary <- fread(snakemake@input$primary, sep = '\t', header = T, select = c(coord_cols, stat_cols, rsid_col))

primary[, (chr_col) := as.character(get(chr_col))]
primary[, (bp_col) := as.integer(get(bp_col))]

prune <- fread(snakemake@input$prune, header = F, col.names = 'ID')

prune[, c(chr_col, bp_col, alt_col, ref_col) := tstrsplit(ID, split = ':', keep = 1:4)]
prune[, (chr_col) := as.character(get(chr_col))]
prune[, (bp_col) := as.integer(get(bp_col))]

prune[, prune_in := TRUE]

primary <- merge(primary, prune[, -c('ID')], all.x = T, by = coord_cols)
primary[is.na(prune_in), prune_in := FALSE]

for(i in seq_along(snakemake@input$informative)) {
  informative <- fread(snakemake@input$informative[[i]])

  informative[, (chr_col) := as.character(get(chr_col))]
  informative[, (bp_col) := as.integer(get(bp_col))]

  if("chi2_expected" %in% names(informative)) {
    setnames(informative, "chi2_expected", p_col)
  }

  informative <- informative[!is.na(get(p_col))]

  primary <- merge(primary, informative[, ..cols_to_keep], by = coord_cols, suffixes = c('', sprintf(".%d", i)))
}

# Drop a very extended definition of the MHC
primary <- primary[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]

primary <- unique(primary, by = c(chr_col, bp_col, ref_col, alt_col))

fwrite(primary, file = snakemake@output[[1]], sep = '\t')
