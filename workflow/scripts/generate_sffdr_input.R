library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col
rsid_col <- snakemake@config$rsid_col

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)
stat_cols <- c(beta_col, se_col, p_col)

primary <- fread(snakemake@input$primary, sep = '\t', header = T, select = c(coord_cols, stat_cols, rsid_col))

primary[, (chr_col) := as.character(get(chr_col))]
primary[, (bp_col) := as.integer(get(bp_col))]

prune <- fread(snakemake@input$prune, header = F, col.names = 'ID')

prune[, c('chromosome', 'base_pair_location', 'other_allele', 'effect_allele') := tstrsplit(ID, split = ':', keep = 1:4)]
prune[, (chr_col) := as.character(get(chr_col))]
prune[, (bp_col) := as.integer(get(bp_col))]

prune[, prune_in := TRUE]

primary <- merge(primary, prune[, -c('ID')], all.x = T, by = coord_cols)
primary[is.na(prune_in), prune_in := FALSE]

for(i in seq_along(snakemake@input$informative)) {
  informative <- fread(snakemake@input$informative[[i]], select = c(coord_cols, stat_cols))

  informative[, (chr_col) := as.character(get(chr_col))]
  informative[, (bp_col) := as.integer(get(bp_col))]

  informative <- informative[!is.na(get(p_col))]

  primary <- merge(primary, informative, by = coord_cols, suffixes = c('', sprintf(".%d", i)))
}

# Drop a very extended definition of the MHC
primary <- primary[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]

primary <- unique(primary, by = c(chr_col, bp_col, ref_col, alt_col))

fwrite(primary, file = snakemake@output[[1]], sep = '\t')
