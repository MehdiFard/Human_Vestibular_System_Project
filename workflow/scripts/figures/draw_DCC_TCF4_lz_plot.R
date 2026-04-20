library(data.table)
library(locuszoomr)
library(ensembldb)
library(patchwork)
library(ggplot2)
theme_set(theme_bw())

chr_col  <- snakemake@config[["col_names"]][["chromosome"]]
bp_col   <- snakemake@config[["col_names"]][["base_pair_location"]]
rsid_col <- snakemake@config[["col_names"]][["rsid"]]
p_col    <- snakemake@params$p_col

chr    <- "18"
start_bp <- 52e6
end_bp   <- 57e6

snps_to_label <- snakemake@params$snps_to_label
if (length(snps_to_label) == 0) snps_to_label <- NULL

dat <- fread(snakemake@input$gwas, sep = '\t', header = TRUE,
             select = c(chr_col, bp_col, rsid_col, p_col))

dat[, (chr_col) := as.character(get(chr_col))]

dat <- dat[get(chr_col) == chr & get(bp_col) %between% c(start_bp, end_bp)]

dat <- na.omit(dat, cols = c(p_col, bp_col))

edb <- EnsDb(snakemake@input$edb)

loc <- locus(
    data       = dat,
    chrom      = chr_col,
    pos        = bp_col,
    p          = p_col,
    labs       = rsid_col,
    seqname    = chr,
    xrange     = c(start_bp, end_bp),
    index_snp  = snps_to_label[1],
    ens_db     = edb
)

p_scatter <- gg_scatter(loc, showLD = FALSE, labels = snps_to_label,
                        min.segment.length = 0, nudge_y = 1)+ylab("-log10 functional\np-value")+theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ylim(0, 11)

ggsave(p_scatter / gg_genetracks(loc, filter_gene_biotype = 'protein_coding') + plot_layout(heights = c(3, 1)),
       file = snakemake@output[[1]], width = 6, height = 6)
