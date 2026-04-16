library(data.table)
library(locuszoomr)
library(ensembldb)

chr_col  <- snakemake@config[["col_names"]][["chromosome"]]
bp_col   <- snakemake@config[["col_names"]][["base_pair_location"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]
p_col    <- snakemake@params$p_col

chrom    <- snakemake@wildcards$chrom
start_bp <- as.integer(snakemake@wildcards$start)
end_bp   <- as.integer(snakemake@wildcards$end)

dat <- fread(snakemake@input$gwas, sep = '\t', header = TRUE,
             select = c(chr_col, bp_col, rsid_col, p_col))

dat[[chr_col]] <- as.character(dat[[chr_col]])

dat <- dat[get(chr_col) == chrom & get(bp_col) %between% c(start_bp, end_bp)]

dat <- na.omit(dat, cols = c(p_col, bp_col))

edb <- EnsDb(snakemake@input$edb)

loc <- locus(
    data     = dat,
    chrom    = chr_col,
    pos      = bp_col,
    p        = p_col,
    labs     = rsid_col,
    seqname  = chrom,
    xrange   = c(start_bp, end_bp),
    ens_db   = edb
)

ggsave(gg_scatter(loc, showLD = FALSE),
       file = snakemake@output[[1]], width = 6, height = 4)
