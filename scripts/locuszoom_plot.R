
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
summary(loc)
locus_plot(loc)



### My data - balance_meta1
gwas_meta <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_meta)

# Clean GWAS data
gwas_meta_cln <- gwas_meta %>%
                  select(rsID, chr = CHR, pos = BP, A1, A2, beta = BETA, se = SE, p = P)
head(gwas_meta_cln)


# Finding lead snp rsIDs
gwas_meta[which(gwas_meta$SNP == "19_23421111_A_G"), ]
# rs695017
gwas_meta[which(gwas_meta$SNP == "19_23307401_A_G"), ]
# rs167154 - isoform QTL
gwas_meta[which(gwas_meta$SNP == "19_23359488_A_G"), ]
# rs428549 - RNA_Editing QTL


# Using AnnotationHub to get EnsDb objects
ah <- AnnotationHub()
query(ah, c("EnsDb", "Homo sapiens"))
ensDb_v113 <- ah[["AH119325"]]
summary(ensDb_v113@tables)
ensDb_v113@tables

# Extracting a subset of the data around gene or region of interest
loc_meta <- locus(data = gwas_meta_cln,
                  # gene = 'ZNF91',
                  # flank = 2.5e5,
                  index_snp = "rs695017",
                  flank = 2.5e5,
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
loc_meta_ld_recom <- link_recomb(loc_meta_ld)

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
png("figures/locuszoom_gwas.png",
    width = 12, height = 10, units = "in", res = 300)

locus_plot(loc_meta_ld_recom_eqtls, 
          labels = c("index", "rs167154"), # rs612285, rs426564
          label_x = c(2, -5),
          highlight = c("ZNF724", "ZNF91", "ZNF675"),
          # filter_gene_biotype = "protein_coding",
          gene_col = 'grey',
          exon_col = 'orange',
          exon_border = 'darkgrey'
          )

dev.off()


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


#------- Our isoformQTLs
qtl <- fread("Pantry/pheast_no_outlier/output/vestibular/qtl/isoforms.cis_qtl.txt.gz")  # change to your file
head(qtl)
dim(qtl)


# Check if our TWAS top SNP is present
which(qtl$variant_id == "19_23307401_A_G")
# Not present in isoform QTLs


library(tidyr)
qtl_cln <- qtl %>%
  select(SNP = variant_id, beta = slope, se = slope_se, p = qval) %>%
  # select(phenotype_id, group_id, snp_chr_pos = variant_id, p = qval) %>%
  separate(col = SNP, into  = c("chr", "pos", "A1", "A2"), sep   = "_", remove = FALSE) %>%
  mutate(
    pos = as.integer(pos)) %>%
  inner_join(
    gwas_meta[, c("SNP", "rsID")], by = "SNP") %>%
  select(rsID, chr, pos, A1, A2, beta, se, p)   # enforce final order

dim(qtl_cln)
head(qtl_cln)
which(qtl_cln$rsID == "rs167154")


# Extracting a subset of the data around gene or region of interest
loc_meta <- locus(data = qtl_cln,
                  gene = 'ZNF91',
                  # index_snp = "rs167154",
                  # flank = 2.5e5,
                  flank = 1e6,
                  ens_db = ensDb_v113,
                  )
summary(loc_meta)

png("figures/locuszoom_isoform.png",
    width = 12, height = 10, units = "in", res = 300)

locus_plot(loc_meta,
                  labels = c("rs426564"),
                  label_x = c(2, -5),
                  ylab = expression(-log[10]~qvalue),
                  highlight = "ZNF91",
                  )

dev.off()

# SNP from TWAS
# 19_23307401_A_G, rs167154


# FROM isoform QTL mapping
$index_snp
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