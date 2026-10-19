library(data.table)
library(locuszoomr)
library(ensembldb)
library(patchwork)
library(ggplot2)
theme_set(theme_bw())
save.image('dcc.RData')
stop()

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
    ens_db     = edb,
  )

loc <- link_LD(loc, token = snakemake@config$apis$ldlink$token, pop = "EUR")
loc <- link_recomb(loc)

# NB: necessary to customise y axes.
rscale <- max(loc$recomb$value, na.rm = TRUE) / max(loc$data$logP, na.rm = TRUE)

p_scatter <- gg_scatter(loc, showLD = TRUE, labels = snps_to_label,
                        min.segment.length = 0, nudge_y = 1)+ylim(0, 11)+scale_y_continuous(
     name    = expression(-log[10](plain("functional p-value"))),
     sec.axis = sec_axis(
       ~ . * rscale,
       name   = "Recombination rate (cM/Mb)",
       breaks = seq(0, 100, 25)
     )
   )+theme(legend.background = element_rect(fill = "transparent"))

ggsave(p_scatter / gg_genetracks(loc, filter_gene_biotype = 'protein_coding') + plot_layout(heights = c(3, 1)),
       file = snakemake@output[[1]], width = 6, height = 6)
