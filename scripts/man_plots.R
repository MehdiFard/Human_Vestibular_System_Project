### Making manhattan plots

rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(qqman)


### Manhattan plot - hg38 result
gwas_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta1_formatted_hg38_ChrPosId.tsv.gz")
dim(gwas_hg38)

gwas_cln_hg38 <- gwas_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(CHR %in% 1:22) %>%
    # filter(P < 0.01) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_cln_hg38)
unique(gwas_cln_hg38$CHR)

# Threshold
fdr_threshold <- 0.05/dim(gwas_hg38)[1]

# Manhattan plot
w_mm <- 188
h_mm <- 94

pdf("figures/man_plot_meta_balance_hg38.pdf", height = h_mm/25.4, width = w_mm/25.4) # 4H * 8W

manhattan(
    gwas_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = FALSE # -log10(1e-5)
)

dev.off()

gwas_cln_hg38_sig <- gwas_cln_hg38 %>%
    dplyr::filter(CHR %in% c(17, 18, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_cln_hg38_sig)
head(gwas_cln_hg38_sig)
fwrite(gwas_cln_hg38_sig, "figures/data_tables/balance_meta_top_hits_hg38.tsv", sep = "\t")




### QQ plot
w_mm <- 70
h_mm <- 70

p <- gwas_cln_hg38$P
p <- p[is.finite(p) & p > 0 & p <= 1]

p_sorted <- sort(p, decreasing = FALSE)
n_p <- length(p_sorted)
max_points <- 50000

if (n_p > max_points) {
    idx <- unique(c(
        as.integer(seq(1, n_p, length.out = max_points)),
        which(p_sorted < 1e-4)
    ))
} else {
    idx <- seq_len(n_p)
}

qq_df <- data.frame(
    expected = -log10(ppoints(n_p)[idx]),
    observed = -log10(p_sorted[idx])
)

qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.4, alpha = 0.5, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
        # title = "Balance Meta-(Lambda = 1.09)",
        x = expression(Expected~~-log[10](italic(P))),
        y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
        axis.title = element_text(size = rel(0.6)),
        axis.text  = element_text(size = rel(0.6))
    )

ggsave(
    filename = "figures/qq_plot_meta_balance_lambda_1.12.pdf",
    plot = qq_plot,
    width = w_mm,
    height = h_mm,
    units = "mm",
    device = "pdf"
)

lambda <- median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda


###----------- Manhattan plot - hg19 result
gwas_hg19 <- fread("data/gwas_raw_sumstats/balance_meta1.tbl.sorted.filtered.gz")

# Clean and filter the GWAS data
gwas_cln_hg19 <- gwas_hg19 %>%
    select(CHR = Chromosome, BP = Position, P = "P-value", SNP = MarkerName) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))

# Threshold
fdr_threshold <- 0.05/dim(gwas_hg19)[1]

w_mm <- 180
h_mm <- 90

pdf("figures/man_plot_meta_balance_hg19.pdf", height = h_mm/25.4, width = w_mm/25.4) # 4H * 8W

manhattan(
    gwas_cln_hg19, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    cex = 0.4,
    cex.axis = 0.6,
    main = "Balance Meta-All",
    cex.main = 0.6,
    cex.lab = 0.7,
    genomewideline = FALSE, # -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = FALSE
)
abline(h = -log10(fdr_threshold), col = "red", lty = 2) # , lwd = 0.8
dev.off()


### QQ plot
gwas_cln_hg19_qq <- gwas_hg19 %>%
    select(CHR = Chromosome, BP = Position, P = "P-value", SNP = MarkerName) %>%
    filter(CHR %in% 1:22) %>%
    mutate(CHR = as.integer(CHR))

p <- gwas_cln_hg19_qq$P
lambda <- median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda


### QQ plot
w_mm <- 60
h_mm <- 60

p <- gwas_cln_hg19_qq$P
p <- p[is.finite(p) & p > 0 & p <= 1]

p_sorted <- sort(p, decreasing = FALSE)
n_p <- length(p_sorted)
max_points <- 50000

if (n_p > max_points) {
    idx <- unique(c(
        as.integer(seq(1, n_p, length.out = max_points)),
        which(p_sorted < 1e-4)
    ))
} else {
    idx <- seq_len(n_p)
}

qq_df <- data.frame(
    expected = -log10(ppoints(n_p)[idx]),
    observed = -log10(p_sorted[idx])
)

qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.4, alpha = 0.5, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
        # title = "Balance Meta-(Lambda = 1.09)",
        x = expression(Expected~~-log[10](italic(P))),
        y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
        axis.title = element_text(size = rel(0.6)),
        axis.text  = element_text(size = rel(0.6))
    )

ggsave(
    filename = "figures/qq_plot_meta_balance_all_lambda_1.12.pdf",
    plot = qq_plot,
    width = w_mm,
    height = h_mm,
    units = "mm",
    device = "pdf"
)


###----------- Manhattan plot - meta_Eur result - hg38
gwas_meta_eur_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/balance_meta_eur_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_meta_eur_hg38)
dim(gwas_meta_eur_hg38)

# Clean and filter the GWAS data
gwas_meta_eur_cln_hg38 <- gwas_meta_eur_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(CHR %in% 1:22) %>%
    filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_meta_eur_cln_hg38)

# Save significant hits
gwas_meta_eur_sig_hg38 <- gwas_meta_eur_cln_hg38 %>%
    dplyr::filter(CHR %in% c(17, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

# Write to file
fwrite(gwas_meta_eur_sig_hg38, "figures/data_tables/balance_meta_eur_top_hits_hg38.tsv", sep = "\t")

# Threshold
fdr_threshold <- 0.05/dim(gwas_meta_eur_hg38)[1]

# Manhattan plot
w_mm <- 180
h_mm <- 90

pdf("figures/man_plot_balance_meta_eur_hg38.pdf", height = h_mm/25.4, width = w_mm/25.4) # 4H * 8W

manhattan(
    gwas_meta_eur_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    cex = 0.4,
    cex.axis = 0.6,
    main = "Balance Meta-EUR",
    cex.main = 0.6,
    cex.lab = 0.7,
    genomewideline = FALSE, # -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = FALSE
)
abline(h = -log10(fdr_threshold), col = "red", lty = 2) # , lwd = 0.8
dev.off()


### QQ plot
w_mm <- 60
h_mm <- 60

p <- gwas_meta_eur_cln_hg38$P
p <- p[is.finite(p) & p > 0 & p <= 1]

p_sorted <- sort(p, decreasing = FALSE)
n_p <- length(p_sorted)
max_points <- 50000

if (n_p > max_points) {
    idx <- unique(c(
        as.integer(seq(1, n_p, length.out = max_points)),
        which(p_sorted < 1e-4)
    ))
} else {
    idx <- seq_len(n_p)
}

qq_df <- data.frame(
    expected = -log10(ppoints(n_p)[idx]),
    observed = -log10(p_sorted[idx])
)

qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.4, alpha = 0.5, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
        # title = "Balance Meta-(Lambda = 1.09)",
        x = expression(Expected~~-log[10](italic(P))),
        y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
        axis.title = element_text(size = rel(0.6)),
        axis.text  = element_text(size = rel(0.6))
    )

ggsave(
    filename = "figures/qq_plot_meta_balance_eur_lambda_1.11.pdf",
    plot = qq_plot,
    width = w_mm,
    height = h_mm,
    units = "mm",
    device = "pdf"
)

lambda <- median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda



###----------- Manhattan plot - UKBB result - hg19
gwas_ukb <- fread("data/gwas_raw_sumstats/Royce_files/UKB_EUR_balance.txt.gz")
head(gwas_ukb)

# Clean and filter the GWAS data
gwas_ukb_cln_hg19 <- gwas_ukb %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(P < 1e-1) %>%
    mutate(CHR = as.integer(CHR))

# Save significant hits
gwas_ukb_cln_sig <- gwas_ukb_cln_hg19 %>%
    dplyr::filter(CHR %in% c(3, 13, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)
# Write to file
fwrite(gwas_ukb_cln_sig, "figures/data_tables/balance_ukb_top_hits.tsv", sep = "\t")


# Threshold
fdr_threshold <- 0.05/dim(gwas_ukb)[1]

# Manhattan plot
w_mm <- 180
h_mm <- 90

pdf("figures/man_plot_balance_ukb_hg19.pdf", height = h_mm/25.4, width = w_mm/25.4) # 4H * 8W

manhattan(
    gwas_ukb_cln_hg19, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    cex = 0.4,
    cex.axis = 0.6,
    main = "Balance UKBB",
    cex.main = 0.6,
    genomewideline = FALSE, # -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = FALSE
)
abline(h = -log10(fdr_threshold), col = "red", lty = 2) # , lwd = 0.8
dev.off()


# Clean for qq plot and lambda calculation
gwas_ukb_cln_hg19_qq <- gwas_ukb %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    filter(CHR %in% 1:22) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_ukb_cln_hg19_qq)


### QQ plot
w_mm <- 60
h_mm <- 60

p <- gwas_ukb_cln_hg19_qq$P
p <- p[is.finite(p) & p > 0 & p <= 1]

p_sorted <- sort(p, decreasing = FALSE)
n_p <- length(p_sorted)
max_points <- 50000

if (n_p > max_points) {
    idx <- unique(c(
        as.integer(seq(1, n_p, length.out = max_points)),
        which(p_sorted < 1e-4)
    ))
} else {
    idx <- seq_len(n_p)
}

qq_df <- data.frame(
    expected = -log10(ppoints(n_p)[idx]),
    observed = -log10(p_sorted[idx])
)

qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.4, alpha = 0.5, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
        # title = "Balance Meta-(Lambda = 1.09)",
        x = expression(Expected~~-log[10](italic(P))),
        y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
        axis.title = element_text(size = rel(0.6)),
        axis.text  = element_text(size = rel(0.6))
    )

ggsave(
    filename = "figures/qq_plot_UKB_eur_lambda_1.09.pdf",
    plot = qq_plot,
    width = w_mm,
    height = h_mm,
    units = "mm",
    device = "pdf"
)

lambda <- median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda


###----------- Manhattan plot - UKBB result - hg38
gwas_ukb_hg38 <- fread("gwas_sumstats_prep/formatted_sumstats/UKB_EUR_balance_formatted_hg38_ChrPosId.tsv.gz")
head(gwas_ukb_hg38)

# Clean and filter the GWAS data
gwas_ukb_cln_hg38 <- gwas_ukb_hg38 %>%
    select(CHR = CHR, BP = BP, P = P, SNP = SNP) %>%
    # filter(P < 1e-2) %>%
    mutate(CHR = as.integer(CHR))
dim(gwas_ukb_cln_hg38)

# Threshold
fdr_threshold <- 0.05/dim(gwas_ukb_hg38)[1]

# Manhattan plot
pdf("figures/man_plot_balance_ukb_hg38.pdf", height = 6, width = 10)

manhattan(
    gwas_ukb_cln_hg38, # subset for plotting
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray60", "lightpink3"),
    genomewideline = -log10(fdr_threshold), # -log10(5e-8),
    suggestiveline = -log10(1e-5)
)

dev.off()

gwas_ukb_cln_hg38_sig <- gwas_ukb_cln_hg38 %>%
    dplyr::filter(CHR %in% c(13, 19)) %>%
    dplyr::filter(P < 1e-5) %>%
    dplyr::arrange(CHR, P)

dim(gwas_ukb_cln_hg38_sig)
head(gwas_ukb_cln_hg38_sig)
fwrite(gwas_ukb_cln_hg38_sig, "figures/data_tables/balance_ukb_top_hits_hg38.tsv", sep = "\t")


## QQ
qq(gwas_ukb_cln_hg38$P, main = "QQ plot - balance UKBB hg38")

p <- gwas_ukb_cln_hg38$P
p <- p[!is.na(p) & p > 0 & p <= 1]
chisq <- qchisq(1 - p, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
lambda

lambda <- median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda

n_case = 46178
n_control = 317856
lambda1000 <- 1 + (lambda - 1) *
  ((1 / n_case + 1 / n_control) / (1 / 1000 + 1 / 1000))
lambda1000


### END

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Input
infile <- "data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz"

# Read MAGMA gene-level results
magma <- fread(infile)

# Check required columns
required_cols <- c("GENE", "CHR", "START", "STOP", "P", "SYMBOL", "set")
stopifnot(all(required_cols %in% colnames(magma)))

# Optional: choose one result set
# Your file has a column called "set", for example UKB, meta_EUR, meta_all.
# Change this if you want a different subset.
plot_set <- "UKB"

gwas_data <- magma %>%
#   filter(set == plot_set) %>%
  mutate(
    CHR = as.character(CHR),
    CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
    CHR = ifelse(CHR == "X", "23", CHR),
    CHR = as.numeric(CHR),
    BP = as.numeric(START),
    P = as.numeric(P),
    SYMBOL = as.character(SYMBOL)
  ) %>%
  filter(
    !is.na(CHR),
    CHR >= 1,
    CHR <= 22,
    !is.na(BP),
    !is.na(P),
    P > 0
  ) %>%
  arrange(CHR, BP)

# Cumulative genomic position
data_cum <- gwas_data %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHR, bp_add)

gwas_data <- gwas_data %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

# Chromosome axis centers
axis_set <- gwas_data %>%
  group_by(CHR) %>%
  summarise(center = mean(bp_cum, na.rm = TRUE), .groups = "drop")

# Bonferroni threshold
# Use number of tested genes in this subset rather than hard-coding 19042
n_genes <- nrow(gwas_data)
sig <- 0.05 / n_genes

# Significant genes
sig_genes <- gwas_data %>%
  filter(P <= sig) %>%
  arrange(P)

# Optional: label only the most significant genes if there are too many
# For example, label all significant genes, but cap at top 25.
label_genes <- sig_genes %>%
  slice_head(n = 25)


# Plot
magma_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
  geom_point(
    aes(color = as.factor(CHR %% 2)),
    alpha = 0.55,
    size = 0.45
  ) +
  scale_color_manual(values = c("grey60", "lightpink3")) +
  geom_hline(
    yintercept = -log10(sig),
    color = "red",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    data = sig_genes,
    aes(x = bp_cum, y = -log10(P)),
    color = "black",
    size = 0.7
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = SYMBOL),
    size = 2.6,
    box.padding = 0.30,
    point.padding = 0.18,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "red",
    segment.size = 0.28,
    segment.alpha = 0.8,
    segment.linetype = "dashed",
    seed = 123
  ) +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = axis_set$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 7),
    plot.title = element_text(hjust = 0.5),

    axis.title.x = element_text(size = 8, face = "plain"),
    axis.title.y = element_text(size = 8, face = "plain"),

    # axis.line.x = element_line(color = "black", linewidth = 0.3),
    # axis.line.y = element_line(color = "black", linewidth = 0.3)
  )
magma_plot

ggsave(
    filename = "figures/MAGMA_gene_manhattan.pdf",
    plot = magma_plot,
    width = w_mm,
    height = h_mm,
    units = "mm",
    device = "pdf",
    dpi = 300
)



#### By qqman

library(data.table)
library(dplyr)
library(qqman)

magma <- fread("data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz")

magma_qq <- magma %>%
    mutate(
        CHR = as.character(CHR),
        CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
        CHR = ifelse(CHR == "X", "23", CHR),
        CHR = as.numeric(CHR),
        BP = as.numeric(START),
        P = as.numeric(P),
        SNP = as.character(SYMBOL)
    ) %>%
    filter(
        !is.na(CHR),
        CHR >= 1,
        CHR <= 22,
        !is.na(BP),
        !is.na(P),
        P > 0,
        !is.na(SNP)
    ) %>%
    arrange(CHR, BP) %>%
    select(SNP, CHR, BP, P)

sig <- 0.05 / nrow(magma_qq)


w_mm <- 180
h_mm <- 90

pdf("figures/MAGMA_qqman_all_sets1.pdf", height = h_mm/25.4, width = w_mm/25.4) # 4H * 8W

manhattan(
    magma_qq,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("grey60", "lightpink3"),
    # cex = 0.4,
    cex.axis = 0.6,
    genomewideline = FALSE, # -log10(sig),
    suggestiveline = FALSE,
    annotatePval = sig,
    annotateTop = FALSE,
    cex.axis = 0.7,
    cex.lab = 0.8,
    ylab = expression(-log[10](P)),
    xlab = "Chromosome",
    main = ""
)
abline(h = -log10(sig), col = "red", lty = 2) # , lwd = 0.8
dev.off()


### Trying to make the lables smaller and less crowded, but qqman doesn't have as much flexibility as ggplot2 for this.
library(data.table)
library(dplyr)
library(qqman)

#-----------------------------
# Input / output
#-----------------------------
infile <- "data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

outfile <- "figures/MAGMA_qqman_all_sets.pdf"

w_mm <- 180
h_mm <- 90

#-----------------------------
# Read data
#-----------------------------
magma <- fread(infile)

#-----------------------------
# Prepare data for qqman
# No filtering by set
#-----------------------------
magma_qq <- magma %>%
  mutate(
    CHR = as.character(CHR),
    CHR = gsub("^chr", "", CHR, ignore.case = TRUE),
    CHR = ifelse(CHR == "X", "23", CHR),
    CHR = as.numeric(CHR),
    BP  = as.numeric(START),
    P   = as.numeric(P),
    SNP = as.character(SYMBOL)
  ) %>%
  filter(
    !is.na(CHR),
    CHR >= 1,
    CHR <= 22,
    !is.na(BP),
    !is.na(P),
    P > 0,
    !is.na(SNP),
    SNP != ""
  ) %>%
  arrange(CHR, BP) %>%
  select(SNP, CHR, BP, P)

# Bonferroni threshold across all rows being plotted
sig <- 0.05 / nrow(magma_qq)

#-----------------------------
# Recreate cumulative positions for manual labels
#-----------------------------
chr_info <- magma_qq %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(bp_add = lag(cumsum(chr_len), default = 0))

plot_df <- magma_qq %>%
  left_join(chr_info, by = "CHR") %>%
  mutate(
    pos  = BP + bp_add,
    logp = -log10(P)
  ) %>%
  arrange(CHR, BP)

#-----------------------------
# Significant genes to label
#-----------------------------
label_df <- plot_df %>%
  filter(P <= sig) %>%
  arrange(CHR, BP, P)

# Optional cap if too many labels:
# label_df <- plot_df %>%
#   filter(P <= sig) %>%
#   arrange(P) %>%
#   slice_head(n = 30)

#-----------------------------
# Label placement
#-----------------------------
x_min   <- min(plot_df$pos, na.rm = TRUE)
x_max   <- max(plot_df$pos, na.rm = TRUE)
x_range <- x_max - x_min

n_tiers       <- 6
x_offset_base <- 0.008 * x_range
x_offset_step <- 0.004 * x_range
y_offset_base <- 0.7
y_offset_step <- 0.45
x_buffer      <- 0.01 * x_range

if (nrow(label_df) > 0) {
  label_df <- label_df %>%
    group_by(CHR) %>%
    mutate(
      idx_chr = row_number(),
      tier    = ((idx_chr - 1) %% n_tiers) + 1,
      side    = ifelse(idx_chr %% 2 == 1, 1, -1),
      label_x = pos + side * (x_offset_base + (tier - 1) * x_offset_step),
      label_y = logp + y_offset_base + (tier - 1) * y_offset_step
    ) %>%
    ungroup() %>%
    mutate(
      label_x = pmax(x_min + x_buffer, pmin(label_x, x_max - x_buffer))
    )
}

ylim_top <- max(
  c(plot_df$logp, if (nrow(label_df) > 0) label_df$label_y else 0),
  na.rm = TRUE
) + 0.5

#-----------------------------
# Plot
#-----------------------------
pdf(outfile, height = h_mm / 25.4, width = w_mm / 25.4)

par(
  mar = c(3.2, 3.2, 0.8, 0.6),
  mgp = c(1.8, 0.45, 0),
  tcl = -0.2,
  xpd = NA
)

manhattan(
  magma_qq,
  chr = "CHR",
  bp  = "BP",
  p   = "P",
  snp = "SNP",
  col = c("grey60", "lightpink3"),
  cex = 0.4,
  cex.axis = 0.65,
  cex.lab  = 0.8,
  genomewideline = FALSE,
  suggestiveline = FALSE,
  annotatePval = NULL,
  ylim = c(0, ylim_top),
  ylab = expression(-log[10](P)),
  xlab = "Chromosome",
  main = ""
)

abline(h = -log10(sig), col = "red", lty = 2, lwd = 0.8)

if (nrow(label_df) > 0) {
  segments(
    x0 = label_df$pos,
    y0 = label_df$logp,
    x1 = label_df$label_x,
    y1 = label_df$label_y - 0.12,
    col = "red",
    lty = 2,
    lwd = 0.5
  )

  text(
    x = label_df$label_x,
    y = label_df$label_y,
    labels = label_df$SNP,
    cex = 0.55,
    font = 1,
    col = "black"
  )
}

dev.off()



###--------------------------- CMPlots
# install.packages("CMplot")
library("CMplot")

data(pig60K)   #calculated p-values by MLM
data(cattle50K)   #calculated SNP effects by rrblup
head(pig60K)


meta_all <- fread("data/gwas_raw_sumstats/balance_meta1.tbl.sorted.filtered.gz") %>%
    select(SNP = MarkerName, Chromosome = Chromosome, Position = Position, Mate_All = "P-value")
head(meta_all)
dim(meta_all)

meta_eur <- fread("data/gwas_raw_sumstats/Royce_files/balance_meta_EUR.tbl.sorted.tab.gz") %>%
    select(SNP = SNP, Chromosome = Chromosome, Position = Position, Mate_Eur = "P-value")
head(meta_eur)
dim(meta_eur)


gwas_ukb <- fread("data/gwas_raw_sumstats/Royce_files/UKB_EUR_balance.txt.gz") %>%
    select(SNP = SNP, Chromosome = CHR, Position = BP, GWAS_UKB = P)
head(gwas_ukb)
dim(gwas_ukb)

table(meta_all$SNP %in% meta_eur$SNP)


library(dplyr)
library(data.table)
library(tidyr)

meta_all <- fread("data/gwas_raw_sumstats/balance_meta1.tbl.sorted.filtered.gz") %>%
  select(SNP = MarkerName, Chromosome = Chromosome, Position = Position, Mate_All = "P-value") %>%
  mutate(
    SNP = as.character(SNP),
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position)
  ) %>%
  distinct(SNP, Chromosome, Position, .keep_all = TRUE)


meta_eur <- fread("data/gwas_raw_sumstats/Royce_files/balance_meta_EUR.tbl.sorted.tab.gz") %>%
  select(SNP = SNP, Chromosome = Chromosome, Position = Position, Mate_Eur = "P-value") %>%
  mutate(
    SNP = as.character(SNP),
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position)
  ) %>%
  distinct(SNP, Chromosome, Position, .keep_all = TRUE)


gwas_ukb <- fread("data/gwas_raw_sumstats/Royce_files/UKB_EUR_balance.txt.gz") %>%
  select(SNP = SNP, Chromosome = CHR, Position = BP, GWAS_UKB = P) %>%
  mutate(
    SNP = as.character(SNP),
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position)
  ) %>%
  distinct(SNP, Chromosome, Position, .keep_all = TRUE)



merged_df <- meta_eur %>%
    left_join(meta_all, by = c("SNP", "Chromosome", "Position")) %>%
    left_join(gwas_ukb, by = c("SNP", "Chromosome", "Position")) %>%
    mutate(
        Mate_All = replace_na(Mate_All, 0.90),
        Mate_Eur = replace_na(Mate_Eur, 0.90),
        GWAS_UKB = replace_na(GWAS_UKB, 0.90)
    ) %>%
    select(SNP, Chromosome, Position, Mate_All, Mate_Eur, GWAS_UKB)


unique(merged_df$Chromosome)
dim(merged_df)
head(merged_df)
0.05 / nrow(merged_df)

# Circular Manhattan plot
CMplot(merged_df,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:22),sep=""),
      threshold=c(4.5e-9,1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
      "blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
      bin.size=1e6,outward=FALSE,file="jpg",file.name=NULL,dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)



### Manplot

thr_meta_all <- 0.05 / nrow(meta_all)
thr_meta_eur <- 0.05 / nrow(meta_eur)
thr_gwas_ukb <- 0.05 / nrow(gwas_ukb)


thr = 5e-08

SNPs_meta_all <- merged_df$SNP[merged_df$Mate_All < thr]
SNPs_meta_eur <- merged_df$SNP[merged_df$Mate_Eur < thr]
SNPs_gwas_ukb <- merged_df$SNP[merged_df$GWAS_UKB < thr]

highlight_snps <- list(
    SNPs_meta_all,
    SNPs_meta_eur,
    SNPs_gwas_ukb
)

CMplot(
    merged_df,
    type = "p",
    plot.type = "m",
    LOG10 = TRUE,
    highlight = highlight_snps,
    highlight.type = "l",
    threshold = 5e-08,
    threshold.col = "red",   # red threshold line
    threshold.lty = 2,       # dashed
    col = c("grey60", "#4197d8"),
    signal.cex = 0.6,
    signal.col = "green4", # chocolate2
    highlight.col = "gray40", # highlighted SNPs in gray
    highlight.cex = 0.5,
    file = "jpg",
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    multracks = TRUE
)


# Multiple tracks Q-Q plot - not good as we have a lot of missing values during the merge which replaced by 0.9

set.seed(123)
merged_df_qq <- merged_df
merged_df_qq$Mate_All[sample(seq_len(nrow(merged_df_qq)), round(nrow(merged_df_qq) * 0.80))] <- NA
merged_df_qq$Mate_Eur[sample(seq_len(nrow(merged_df_qq)), round(nrow(merged_df_qq) * 0.25))] <- NA

CMplot(
    merged_df_qq,
    plot.type = "q",
    col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),
    multracks = TRUE,
    threshold = thr,
    ylab.pos = 2,
    signal.pch = c(19, 6, 4),
    signal.cex = 1.2,
    signal.col = "red",
    conf.int = TRUE,
    box = FALSE,
    axis.cex = 2,
    file = "jpg",
    file.name = NULL,
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    ylim = c(0, 8),
    width = 5,
    height = 5
)

###
CMplot(
    meta_all,
    plot.type = "q",
    box = FALSE,
    file = "jpg",
    dpi = 300,
    conf.int = TRUE,
    conf.int.col = NULL,
    threshold.col = "red",
    threshold.lty = 2,
    file.output = TRUE,
    verbose = TRUE,
    width = 5,
    height = 5
)

CMplot(
    meta_eur,
    plot.type = "q",
    box = FALSE,
    file = "jpg",
    dpi = 300,
    conf.int = TRUE,
    conf.int.col = NULL,
    threshold.col = "red",
    threshold.lty = 2,
    file.output = TRUE,
    verbose = TRUE,
    width = 5,
    height = 5
)

CMplot(
    gwas_ukb,
    plot.type = "q",
    box = FALSE,
    file = "jpg",
    dpi = 300,
    conf.int = TRUE,
    conf.int.col = NULL,
    threshold.col = "red",
    threshold.lty = 2,
    file.output = TRUE,
    verbose = TRUE,
    width = 5,
    height = 5
)



####




SNPs <- pig60K[pig60K[,5] < (0.05 / nrow(pig60K)), 1]
genes <- paste("GENE", 1:length(SNPs), sep = "_")
set.seed(666666)

CMplot(
    pig60K[, c(1:3, 5)],
    plot.type = "m",
    LOG10 = TRUE,
    col = c("grey30", "grey60"),
    highlight = SNPs,
    highlight.col = rep(c("green", "blue"), length.out = length(SNPs)),
    highlight.cex = 1,
    highlight.text = genes,
    highlight.text.col = rep("red", length(SNPs)),
    threshold = 0.05 / nrow(pig60K),
    threshold.lty = 2,
    amplify = FALSE,
    file = "jpg",
    file.name = "pig60K_manhattan",
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    width = 14,
    height = 6
)



  ### CMplot version: MAGMA gene-based Manhattan (minimal)
  magma_cm <- fread("data/gwas_raw_sumstats/Royce_files/magma_all_sets.csv.gz") %>%
    transmute(
      GENE = as.character(SYMBOL),
      SNP = make.unique(as.character(SYMBOL)),
      Chromosome = as.integer(gsub("^chr", "", ifelse(as.character(CHR) == "X", "23", as.character(CHR)), ignore.case = TRUE)),
      Position = as.numeric(START),
      MAGMA_P = as.numeric(P)
    ) %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(Chromosome),
      Chromosome >= 1,
      Chromosome <= 22,
      !is.na(Position),
      !is.na(MAGMA_P),
      MAGMA_P > 0
    ) %>%
    arrange(Chromosome, Position)

  sig <- 0.05 / nrow(magma_cm)

  sig_hits <- magma_cm %>%
    filter(MAGMA_P <= sig) %>%
    arrange(MAGMA_P) # %>%
    # slice_head(n = 20)

  CMplot(
    magma_cm[, c("SNP", "Chromosome", "Position", "MAGMA_P")],
    type = "p",
    plot.type = "m",
    LOG10 = TRUE,
    cex = 0.5,
    col = c("grey60", "#4197d8"), # col = c("grey60", "lightpink3"),
    threshold = sig,
    threshold.col = "red",
    threshold.lty = 2,
    amplify = FALSE,
    highlight = sig_hits$SNP,
    highlight.col = "red",
    highlight.cex = 0.6,
    highlight.text = sig_hits$GENE,
    highlight.text.col = "black",
    highlight.text.cex = 0.8,
    file = "jpg",
    file.name = "MAGMA_CMplot_gene_manhattan",
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    width = 14,
    height = 6
  )



### CMplot version: SFFDR results
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
head(sff)


snp_gene_map <- data.table(
  SNP = c("rs7098100", "rs9536418", "rs74564337", "rs4456576", "rs4800996", "rs12980498"),
  Nearest_genes = c("MLLT10", "OLFM4", "C17orf58", "DCC", "TCF4", "ZNF675")
)

sff_cm <- sff %>%
  left_join(snp_gene_map, by = "SNP")

sig_sff <- 0.05 / nrow(sff_cm)

label_hits_sff <- sff_cm %>%
  filter(!is.na(Nearest_genes))

CMplot(
  sff_cm[, c("SNP", "Chromosome", "Position", "SFFDR_P")],
  type = "p",
  plot.type = "m",
  LOG10 = TRUE,
  cex = 0.5,
  col = c("grey60", "#4197d8"),
  threshold = 5e-8, # sig_sff,
  threshold.col = "red",
  threshold.lty = 2,
  amplify = FALSE,
  highlight = label_hits_sff$SNP,
  highlight.col = "red",
  highlight.cex = 0.6,
  highlight.text = label_hits_sff$Nearest_genes,
  highlight.text.col = "black",
  highlight.text.cex = 0.8,
  file = "jpg",
  file.name = "SFFDR_CMplot_gene_manhattan",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 14,
  height = 6
)


### CMplot version: MAGMA + SFFDR in one multi-track figure
magma_sff_multitrack <- full_join(
  magma_cm[, c("SNP", "Chromosome", "Position", "MAGMA_P")],
  sff_cm[, c("SNP", "Chromosome", "Position", "SFFDR_P")],
  by = c("SNP", "Chromosome", "Position")
) %>%
  mutate(
    MAGMA = dplyr::coalesce(MAGMA_P, 0.90),
    SFFDR = dplyr::coalesce(SFFDR_P, 0.90)
  ) %>%
  select(-MAGMA_P, -SFFDR_P) %>%
  arrange(Chromosome, Position)

head(magma_sff_multitrack)

highlight_tracks <- list(
  sig_hits$SNP,
  label_hits_sff$SNP
)

highlight_labels <- list(
  sig_hits$GENE,
  label_hits_sff$Nearest_genes
)


mm_to_in <- function(x) x / 25.4

w_mm <- 250
h_mm <- 90

CMplot(
  magma_sff_multitrack,
  type = "p",
  plot.type = "m",
  LOG10 = TRUE,
  col = c("grey60", "#4197d8"),
  threshold = list(sig, 5e-8), # 5e-8, # c(sig, 5e-8),
  cex = 0.5,
  threshold.col = "red",
  threshold.lty = 2,
  amplify = FALSE,
  highlight = highlight_tracks,
  highlight.col = "red",
  highlight.cex = 0.6,
  highlight.text = highlight_labels,
  highlight.text.col = "black",
  highlight.text.cex = 0.8,
  file = "jpg",
  file.name = "MAGMA_SFFDR_CMplot_multitrack_manhattan",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  multracks = TRUE,
  width = 14, # mm_to_in(w_mm),
  height = 4.724 # 6 # mm_to_in(h_mm)
)


### End