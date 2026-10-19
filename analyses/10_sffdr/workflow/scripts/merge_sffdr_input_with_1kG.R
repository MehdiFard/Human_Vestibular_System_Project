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

primary <- fread(snakemake@input$primary, sep = '\t', header = T, select = c(coord_cols, rsid_col))

primary[, (chr_col) := as.character(get(chr_col))]
primary[, (bp_col) := as.integer(get(bp_col))]

onekgp <- fread(cmd = sprintf("zstdcat %s | grep -v \"^##\"", snakemake@input$onekgp_panel), header = T, skip = "CHROM", select = c("ID", "#CHROM", "POS", "REF", "ALT"))

names(onekgp) <- c("ID", chr_col, bp_col, ref_col, alt_col)
onekgp[, (chr_col) := as.character(get(chr_col))]
onekgp[, (bp_col) := as.integer(get(bp_col))]

primary <- merge(primary, onekgp, by = coord_cols)

for(i in seq_along(snakemake@input$informative)) {
  informative <- fread(snakemake@input$informative[[i]])

  informative[, (chr_col) := as.character(get(chr_col))]
  informative[, (bp_col) := as.integer(get(bp_col))]

  if("chi2_expected" %in% names(informative)) {
    setnames(informative, "chi2_expected", p_col)
  }

  informative <- informative[!is.na(get(p_col))]

  primary <- merge(primary, informative[, ..coord_cols], by = coord_cols)
}

# Drop a very extended definition of the MHC
primary <- primary[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]

primary <- unique(primary, by = c(chr_col, bp_col, ref_col, alt_col))

fwrite(primary[, .(ID)], file = snakemake@output[[1]], sep = '\t')
