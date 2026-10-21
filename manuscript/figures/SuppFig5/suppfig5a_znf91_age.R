## Supplementary Figure 5a: ZNF91 expression (TPM) vs donor age in vestibular samples, with linear fit.
## Run from the project root:  Rscript manuscript/figures/SuppFig5/suppfig5a_znf91_age.R
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2) })
dir.create("manuscript/figures/SuppFig5", showWarnings = FALSE, recursive = TRUE)
expr <- fread("analyses/05_pantry/phenotyping/output/vestibular/unnorm/expression.bed")
expr_znf91 <- expr %>%
  as.data.frame() %>%
  dplyr::filter(phenotype_id == "ENSG00000167232") %>%
  dplyr::select(-c(1:3)) %>%
  data.table::transpose(keep.names = "human_id")

expr_znf91 <- expr_znf91[-1, ]
colnames(expr_znf91) <- c("human_id", "ZNF91")
head(expr_znf91)


mata_date <- fread("analyses/03_sample_matching/results/sample_sheets/metadata_geno.tsv")
mata_date$human_id <- as.character(mata_date$human_id)
dim(mata_date)
head(mata_date)
length(unique(mata_date$human_id))

expr_znf91_age <- merge(mata_date[, -1], expr_znf91, by = "human_id")
expr_znf91_age$ZNF91 <- as.numeric(expr_znf91_age$ZNF91)

head(expr_znf91_age)
expr_znf91_age <- expr_znf91_age[order(expr_znf91_age$ZNF91, decreasing = TRUE), ]

# Fit linear regression: ZNF91 ~ age
lm_znf91 <- lm(as.numeric(ZNF91) ~ age, data = expr_znf91_age)
lm_sum   <- summary(lm_znf91)
lm_beta  <- round(coef(lm_znf91)["age"], 4)
lm_pval  <- round(coef(lm_sum)["age", "Pr(>|t|)"], 4)
lm_r2    <- round(lm_sum$r.squared, 3)
lm_label <- paste0("beta = ", lm_beta, "\np = ", lm_pval, "\nR2 = ", lm_r2)
# Plot ZNF91 expression vs age
ggplot(expr_znf91_age, aes(x = age, y = as.numeric(ZNF91))) + # [-1, ]
  geom_point(color = "deepskyblue4", size = 1, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = lm_label,
            hjust = 1.1, vjust = 1.5, size = 2.5, color = "black") +
  labs(
    title = "ZNF91 Expression vs Age in Vestibular Samples",
    x = "Age",
    y = "ZNF91 Expression (tpm)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    plot.title = element_blank(),
  )

ggsave(
  last_plot(),
  filename = "manuscript/figures/SuppFig5/panels/znf91_expression_vs_age_vestibular.pdf",
  units = "mm",
  width = 80,
  height = 80,
  dpi = 300
)
ggsave(last_plot(), filename = "manuscript/figures/SuppFig5/panels/znf91_expression_vs_age_vestibular.png", units = "mm", width = 80, height = 80, dpi = 300)
cat(lm_label, "\n")


### END