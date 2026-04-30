
rm(list=ls())

library(locuszoomr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(ggrepel)
library(ggplot2)
library(AnnotationHub)
library(data.table)
library(dplyr)


# Test data
data(SLE_gwas_sub)  ## limited subset of data from SLE GWAS
head(SLE_gwas_sub)

# extracting a subset of the data 
loc <- locus(data = SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5,
              ens_db = "EnsDb.Hsapiens.v75")



### My data - balance_meta1
gwas_meta <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_meta)

# Clean GWAS data
gwas_meta_cln <- gwas_meta %>%
                  select(SNP, rsID, chr = CHR, pos = BP, A1, A2, beta = BETA, se = SE, p = P) %>%
                  arrange(p)
head(gwas_meta_cln)

# Top snp in gwas
# 19_23421337_G_A   rs612285

# Finding lead snp rsIDs
# 19_23421111_A_G: colocalized gwas snp (in twas) 
gwas_meta[which(gwas_meta_cln$SNP == "19_23421111_A_G"), ] # rs695017

# 19_23307401_A_G: colocalized isoform qtl (in twas)
gwas_meta[which(gwas_meta$SNP == "19_23307401_A_G"), ] # rs167154


# Using AnnotationHub to get EnsDb objects
ah <- AnnotationHub()
query(ah, c("EnsDb", "Homo sapiens"))
ensDb_v113 <- ah[["AH119325"]]
summary(ensDb_v113@tables)
ensDb_v113@tables

# Extracting a subset of the data around gene or region of interest
loc_meta <- locus(data = gwas_meta_cln,
                  gene = 'ZNF91',
                  # index_snp = "rs695017",
                  flank = 2.5e5,
                  # flank = 1e6,
                  ens_db = ensDb_v113
                  )
summary(loc_meta)

# Link LD information
library(LDlinkR)
loc_meta_ld <- link_LD(loc_meta,
                        token = "7d78d9d5e1ed",
                        pop = "EUR",
                        genome_build = "grch38"
                        )

# Link recombination rates
# options(timeout = 300)



library(rtracklayer)

recomb_gr <- import.bw("src/ucsc_recomb/hg38/recomb1000GAvg.bw")

recomb_gr
class(recomb_gr)

loc_meta_ld_recom <- link_recomb(
  loc_meta_ld,
  genome = "hg38",
  recomb = recomb_gr
)



library(httr)

options(httr.verbose = TRUE)
httr::set_config(httr::verbose())
loc_meta_ld_recom <- link_recomb(loc_meta_ld,
                                genome = "hg38",
                                # table = "recomb1000GAvg"
                                )



# Adding eQTL data
loc_meta_ld_recom_eqtls <- link_eqtl(loc_meta_ld_recom,
                                      token = "7d78d9d5e1ed",
                                      genome_build = "grch38"
                                      )
# Check eQTL data
head(loc_meta_ld_recom_eqtls$LDexp)
str(loc_meta_ld_recom_eqtls$LDexp)
table(loc_meta_ld_recom_eqtls$LDexp$Gene_Symbol)
table(loc_meta_ld_recom_eqtls$LDexp$Tissue)

# From fine-mapping of isoform mapping results: rs426564 and rs10415081
# GWAS plot

# To see gene biotypes
ensembldb::listGenebiotypes(ensDb_v113)


# Size of the plot
w_mm <- 180
h_mm <- 140
pdf("figures/locuszoom_gwas2.pdf",
    width = w_mm / 25.4,
    height = h_mm / 25.4)

locus_plot(loc_meta_ld_recom_eqtls, 
          labels = c("rs695017", "rs167154"), # rs612285: gwas top snp (= "index", ), - rs695017: gwas snp coloced in twas, - rs167154: Isoform qtl coloced in twas
          label_x = c(-1, -5),
          highlight = c("ZNF724", "ZNF91", "ZNF675"),
          # filter_gene_biotype = "protein_coding",
          # gene_col = 'grey',
          # exon_col = 'orange',
          # exon_border = 'darkgrey'
          )

dev.off()





# Ploting
pdf("figures/locuszoom_isoform.pdf",
    width = 8, height = 2)


# Extracting a subset of the data around gene or region of interest
loc_meta_scatter <- locus(data = gwas_meta_cln,
                  gene = 'ZNF91',
                  flank = 1e6,
                  # index_snp = "rs695017",
                  ens_db = ensDb_v113
                  )

pdf("figures/locuszoom_gwas_scatter.pdf",
    width = 8, height = 2)

scatter_plot(loc_meta_scatter,
              labels = c("rs695017", "rs167154"),
              label_x = c(-3, -5),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "black", # outline
              bg = "black",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              # ylim = c(0, 20)
              )

dev.off()


###-------------------------- Gene track
pdf("figures/locuszoom_qtls_gene_track.pdf",
    width = 8, height = 2)

genetracks(loc_meta_isoform,
            highlight = "ZNF91",
            filter_gene_biotype = "protein_coding",
            
            )
dev.off()


###-------------------------- Isoform
library(arrow)
library(tidyr)

qtl_isoform <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/isoforms/vestibular_isoforms.cis_qtl_pairs.19.parquet")
dim(qtl_isoform)
head(qtl_isoform)


# Check if our TWAS top SNP is present
which(qtl_isoform$variant_id == "19_23307401_A_G")


qtl_cln_isoform <- qtl_isoform %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_isoform)
head(qtl_cln_isoform)
which(qtl_cln_isoform$rsID == "rs167154")
length(unique(qtl_cln_isoform$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_isoform <- qtl_cln_isoform %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_isoform)


# Extracting a subset of the data around gene or region of interest
loc_meta_isoform <- locus(data = qtl_cln_filt_isoform,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_isoform)

# Ploting
pdf("figures/locuszoom_isoform.pdf",
    width = 8, height = 2)

scatter_plot(loc_meta_isoform,
              labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#377eb8", # outline
              bg = "#377eb8",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()


# locus_plot(
#   loc_meta_isoform,
#   labels = c("rs167154"),
#   # label_x = c(2, -5),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   # gene_col = 'grey',
#   # exon_col = 'orange',
#   # exon_border = 'darkgrey',
#   cex.text = .5, # control label size
#   col = "#377eb8", # outline
#   bg = "#377eb8",      # dot fill color
#   pcutoff = NULL,
# )


###-------------------------- expression
qtl_expr <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/expression/vestibular_expression.cis_qtl_pairs.19.parquet")
dim(qtl_expr)
View(qtl_expr[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_expr$variant_id == "19_23307401_A_G")


qtl_cln_expr <- qtl_expr %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_expr)
head(qtl_cln_expr)
which(qtl_cln_expr$rsID == "rs167154")
length(unique(qtl_cln_expr$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_expr <- qtl_cln_expr %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_expr)


# Extracting a subset of the data around gene or region of interest
loc_meta_expr <- locus(data = qtl_cln_filt_expr,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_expr)

# Ploting
pdf("figures/locuszoom_expression.pdf",
    width = 8, height = 2)


# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_expr,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#e41a1c", # outline
              bg = "#e41a1c",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()

# locus_plot(
#   loc_meta_expr,
#   # labels = c("rs167154"),
#   # label_x = c(2, -5),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   # gene_col = 'grey',
#   # exon_col = 'orange',
#   # exon_border = 'darkgrey',
#   cex.text = .5, # control label size
#   col = "#377eb8", # outline
#   bg = "#377eb8",      # dot fill color
#   # maxrows = 0,
# )


###-------------------------- alternative polyA QTLs
qtl_alt_polyA <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/alt_polyA/vestibular_alt_polyA.cis_qtl_pairs.19.parquet")
dim(qtl_alt_polyA)
View(qtl_alt_polyA[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_alt_polyA$variant_id == "19_23307401_A_G")


qtl_cln_alt_polyA <- qtl_alt_polyA %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_alt_polyA)
head(qtl_cln_alt_polyA)
which(qtl_cln_alt_polyA$rsID == "rs167154")
length(unique(qtl_cln_alt_polyA$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_alt_polyA <- qtl_cln_alt_polyA %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_alt_polyA)


# Extracting a subset of the data around gene or region of interest
loc_meta_alt_polyA <- locus(data = qtl_cln_filt_alt_polyA,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_alt_polyA)

# Ploting
pdf("figures/locuszoom_alt_polyA.pdf",
    width = 8, height = 2)

# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_alt_polyA,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#ff7f00", # outline
              bg = "#ff7f00",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()

# locus_plot(
#   loc_meta_alt_polyA,
#   # labels = c("rs167154"),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   cex.text = .5, # control label size
#   col = "#ff7f00", # outline
#   bg = "#ff7f00",      # dot fill color
# )


###-------------------------- alternative TSS QTLs
qtl_alt_TSS <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/alt_TSS/vestibular_alt_TSS.cis_qtl_pairs.19.parquet")
dim(qtl_alt_TSS)
View(qtl_alt_TSS[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_alt_TSS$variant_id == "19_23307401_A_G")


qtl_cln_alt_TSS <- qtl_alt_TSS %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_alt_TSS)
head(qtl_cln_alt_TSS)
which(qtl_cln_alt_TSS$rsID == "rs167154")
length(unique(qtl_cln_alt_TSS$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_alt_TSS <- qtl_cln_alt_TSS %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_alt_TSS)


# Extracting a subset of the data around gene or region of interest
loc_meta_alt_TSS <- locus(data = qtl_cln_filt_alt_TSS,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_alt_TSS)

# Ploting
pdf("figures/locuszoom_alt_TSS.pdf",
    width = 8, height = 2)

# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_alt_TSS,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#984ea3", # outline
              bg = "#984ea3",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()

# locus_plot(
#   loc_meta_alt_TSS,
#   # labels = c("rs167154"),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   cex.text = .5, # control label size
#   col = "#984ea3", # outline
#   bg = "#984ea3",      # dot fill color
# )



###-------------------------- RNA editing QTLs
qtl_rna_editing <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/RNA_editing/vestibular_RNA_editing.cis_qtl_pairs.19.parquet")
dim(qtl_rna_editing)
View(qtl_rna_editing[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_rna_editing$variant_id == "19_23307401_A_G")


qtl_cln_rna_editing <- qtl_rna_editing %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_rna_editing)
head(qtl_cln_rna_editing)
which(qtl_cln_rna_editing$rsID == "rs167154")
length(unique(qtl_cln_rna_editing$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_rna_editing <- qtl_cln_rna_editing %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_rna_editing)


# Extracting a subset of the data around gene or region of interest
loc_meta_rna_editing <- locus(data = qtl_cln_filt_rna_editing,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_rna_editing)

# Ploting
pdf("figures/locuszoom_rna_editing.pdf",
    width = 8, height = 2)

# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_rna_editing,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#984210", # outline
              bg = "#984210",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()

# locus_plot(
#   loc_meta_rna_editing,
#   # labels = c("rs167154"),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   cex.text = .5, # control label size
#   col = "#984210", # outline
#   bg = "#984210",      # dot fill color
# )



###-------------------------- splicing QTLs
qtl_splicing <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/splicing/vestibular_splicing.cis_qtl_pairs.19.parquet")
dim(qtl_splicing)
View(qtl_splicing[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_splicing$variant_id == "19_23307401_A_G")


qtl_cln_splicing <- qtl_splicing %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_splicing)
head(qtl_cln_splicing)
which(qtl_cln_splicing$rsID == "rs167154")
length(unique(qtl_cln_splicing$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_splicing <- qtl_cln_splicing %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_splicing)


# Extracting a subset of the data around gene or region of interest
loc_meta_splicing <- locus(data = qtl_cln_filt_splicing,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_splicing)

# Ploting
pdf("figures/locuszoom_splicing.pdf",
    width = 8, height = 2)

# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_splicing,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#4daf4a", # outline
              bg = "#4daf4a",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()

# locus_plot(
#   loc_meta_splicing,
#   # labels = c("rs167154"),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   cex.text = .5, # control label size
#   col = "#4daf4a", # outline
#   bg = "#4daf4a",      # dot fill color
# )



###-------------------------- stability QTLs
qtl_stability <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/stability/vestibular_stability.cis_qtl_pairs.19.parquet")
dim(qtl_stability)
View(qtl_stability[1:10, ])


# Check if our TWAS top SNP is present
which(qtl_stability$variant_id == "19_23307401_A_G")


qtl_cln_stability <- qtl_stability %>%
  # select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = pval_nominal) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln_stability)
head(qtl_cln_stability)
which(qtl_cln_stability$rsID == "rs167154")
length(unique(qtl_cln_stability$rsID))


# Filtering to keep only the most significant association for each SNP
qtl_cln_filt_stability <- qtl_cln_stability %>%
    group_by(rsID) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
dim(qtl_cln_filt_stability)


# Extracting a subset of the data around gene or region of interest
loc_meta_stability <- locus(data = qtl_cln_filt_stability,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta_stability)

# Ploting
pdf("figures/locuszoom_stability.pdf",
    width = 8, height = 2)

# Alternative: scatter plot without gene tracks
scatter_plot(loc_meta_stability,
              # labels = c("rs167154"),
              ylab = expression(-log[10]~pvalue),
              cex = 0.5,           # control dot size
              col = "#fdc11c", # outline
              bg = "#fdc11c",      # dot fill color
              pcutoff = NULL,
              xticks = FALSE,
              ylim = c(0, 8)
              )

dev.off()


# locus_plot(
#   loc_meta_stability,
#   # labels = c("rs167154"),
#   ylab = expression(-log[10]~pvalue),
#   highlight = "ZNF91",   # highlight = c("ZNF724", "ZNF91", "ZNF675"),
#   cex = 0.5,           # control dot size
#   filter_gene_biotype = "protein_coding",
#   cex.text = .5, # control label size
#   col = "#fdc11c", # outline
#   bg = "#fdc11c",      # dot fill color
# )

















#########################

# SNP from TWAS
# 19_23307401_A_G, rs167154


# FROM isoform QTL mapping
# index_snp
[1] "rs426564"
19:23314211 (GRCh38)
ENSG00000167232

###---
bim <- fread("Pantry/pheast_no_outlier/input/vestibular/vestibular_genotype.bim")
dim(bim)
head(bim)

# IsoformQTL identified by twas (19_23307401_A_G, rs167154), check the input genotype file
which(bim$V2 == "19_23307401_A_G")
bim[bim$V2 == "19_23307401_A_G" , ]
# It's present in the genotype file

# But its not present in our isoform QTL results
which(qtl$variant_id == "19_23307401_A_G")
# integer(0)

# But there is another SNP nearby that is significant
which(bim$V2 == "19_23314211_G_C") # Present in genotype file
which(qtl$variant_id == "19_23314211_G_C") # Present in isoform QTL results
qtl[qtl$variant_id == "19_23314211_G_C", ]
which(qtl_cln$rsID == "rs426564")
qtl_cln[which(qtl_cln$rsID == "rs426564"), ]



###------------------------------------ Checking expression

eqtl <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/expression.cis_qtl.txt.gz")  # change to your file
head(eqtl)
dim(eqtl)


library(tidyr)
eqtl_cln <- eqtl %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  # select(phenotype_id, group_id, snp_chr_pos = variant_id, p = qval) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(eqtl_cln)
head(eqtl_cln)


# Extracting a subset of the data around gene or region of interest
loc_meta <- locus(data = eqtl_cln,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  # flank = 2.5e5,
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta)

# Ploting
png("figures/locuszoom_expression.png",
    width = 12, height = 10, units = "in", res = 300)

locus_plot(loc_meta,
                  # labels = c("rs426564"),
                  # label_x = c(2, -5),
                  ylab = expression(-log[10]~qvalue),
                  highlight = "ZNF91",
                  )

dev.off()


###------------------------------------ Checking RNA Editing QTLs
rna_editing_qtl <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/RNA_editing.cis_qtl.txt.gz")  # change to your file
head(rna_editing_qtl)
dim(rna_editing_qtl)


library(tidyr)
rna_editing_qtl_cln <- rna_editing_qtl %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  # select(phenotype_id, group_id, snp_chr_pos = variant_id, p = qval) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(rna_editing_qtl_cln)
head(rna_editing_qtl_cln)


# Extracting a subset of the data around gene or region of interest
loc_meta <- locus(data = rna_editing_qtl_cln,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  # flank = 2.5e5,
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta)

# Ploting
png("figures/locuszoom_RNA_Editing.png",
    width = 12, height = 10, units = "in", res = 300)

locus_plot(loc_meta,
                  # labels = c("rs426564"),
                  # label_x = c(2, -5),
                  ylab = expression(-log[10]~qvalue),
                  highlight = "ZNF91",
                  )

dev.off()


### END



library(arrow)
chr1 <- read_parquet("Pantry/pheast_no_outlier/output/vestibular/qtl/nominal/expression/vestibular_expression.cis_qtl_pairs.1.parquet")
dim(chr1)
View(chr1[1:10, ])
head(chr1)




# Overlay GWAS and eQTL plot
png("figures/locuszoom_gwas_gtex_eqtl_Brain_Cortex.png",
    width = 12, height = 10, units = "in", res = 300)

# overlay eqtl plot
overlay_plot(loc_meta_ld_recom_eqtls, tissue = "Brain - Cortex",  eqtl_gene = "ZNF91",
          labels = c("index", "rs167154", "rs428549"),
          label_x = c(2, -5),
          highlight = "ZNF91",
          # gene_col = 'grey',
          # exon_col = 'orange',
          # exon_border = 'darkgrey'
          )
dev.off()

###------


### Stacked plots butt not working with overlay_plot 
pdf("figures/myplot_test.pdf", width = 6, height = 8)
oldpar <- set_layers(2)  # number of stacked panels
scatter_plot(loc_meta_ld_recom_eqtls, xticks = FALSE,
          labels = c("index", "rs167154", "rs428549"),
          label_x = c(2, -5))

line_plot(loc_meta_ld_recom_eqtls, col = "orange", xticks = FALSE)

genetracks(loc_meta_ld_recom_eqtls,
          highlight = c("ZNF724", "ZNF91", "ZNF675"))  # automatically goes to bottom track

par(oldpar)      # IMPORTANT: restore original par()
dev.off()



### -------------------------- SFFDR results - regional locuszoom plot for chr18:52-57Mb

sff <- fread("/data/hps/assoc/private/gdml/user/twil15/sffdr_for_dizziness/results/sffdr/results.tsv.gz") %>%
  select(SNP = rsID, Chromosome = CHR, Position = BP, SFFDR_P = fp) %>%
  mutate(
    SNP = as.character(SNP),
    Chromosome = as.integer(Chromosome),
    Position = as.numeric(Position),
    SFFDR_P = as.numeric(SFFDR_P)
  ) %>%
  filter(
    !is.na(SNP),
    SNP != "",
    !is.na(Chromosome),
    Chromosome >= 1,
    Chromosome <= 22,
    !is.na(Position),
    !is.na(SFFDR_P),
    SFFDR_P > 0
  ) %>%
  arrange(Chromosome, Position)



# Extracting the region of interest for chromosome 18: 52-57Mb
sff_chr18_reg <- sff %>%
  filter(
    Chromosome == 18,
    Position >= 52000000,
    Position <= 57000000
  ) %>%
  transmute(
    rsid = SNP,
    chr = Chromosome,
    pos = Position,
    p = SFFDR_P
  ) %>%
  arrange(p)

dim(sff_chr18_reg)
head(sff_chr18_reg)

stopifnot(nrow(sff_chr18_reg) > 0)

index_snp_chr18 <- sff_chr18_reg$rsid[which.min(sff_chr18_reg$p)]
index_snp_chr18

# Create locus object for the region of interest
loc_sff_chr18 <- locus(
  data = as.data.frame(sff_chr18_reg),
  seqname = 18,
  xrange = c(52000000, 57000000),
  ens_db = ensDb_v113,
  chrom = "chr",
  pos = "pos",
  p = "p",
  labs = "rsid",
  index_snp = index_snp_chr18
)

# Add LD information using 1000 Genomes EUR reference panel
loc_sff_chr18_ld <- link_LD(
  loc_sff_chr18,
  token = "7d78d9d5e1ed",
  pop = "EUR",
  genome_build = "grch38"
)

# Add recombination rate information using UCSC recombination map for hg38
recomb_gr <- rtracklayer::import.bw("src/ucsc_recomb/hg38/recomb1000GAvg.bw")
loc_sff_chr18_ld_recomb <- link_recomb(
  loc_sff_chr18_ld,
  genome = "hg38",
  recomb = recomb_gr
)

# Plotting the locuszoom plot for the region of interest
jpeg(
  "figures/locuszoom_sffdr_chr18_52_57Mb.jpg",
  width = 140,
  height =140,
  units = "mm", 
  res = 300
)

# Optional gene-track controls
# Set genes_keep to show only specific genes in the track (preferred)
# Set genes_drop to remove specific genes from the track
genes_keep <- c("DCC", "TCF4", "MBD2", "CCDC68", "TXNL1", "WDR7", "POLI")
genes_drop <- NULL

gene_filter_for_plot <- genes_keep

if (is.null(gene_filter_for_plot) && !is.null(genes_drop)) {
  tx_col <- intersect(c("gene_name", "symbol", "SYMBOL", "gene"), names(loc_sff_chr18_ld_recomb$TX))[1]
  if (!is.na(tx_col)) {
    available_genes <- unique(as.character(loc_sff_chr18_ld_recomb$TX[[tx_col]]))
    gene_filter_for_plot <- setdiff(available_genes, genes_drop)
  }
}

# Plot the locuszoom plot with gene tracks, highlighting specific genes and using the filtered gene list
locus_plot(
  loc_sff_chr18_ld_recomb,
  labels = c("index", "rs4800996"),
  # filter_gene_name = gene_filter_for_plot,
  ylab = expression(-log[10](SFFDR~P)),
  filter_gene_biotype = "protein_coding",
  highlight = c("DCC", "TCF4")
)

dev.off()

### END
