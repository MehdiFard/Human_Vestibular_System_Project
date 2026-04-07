library(readr)
library(sffdr)
library(tibble)
library(dplyr)
library(ggplot2)

sumstats <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

informative_cols <- names(sumstats)[grep("p_value.\\d+", names(sumstats))]

p_cols <- c("p_value", informative_cols)

sumstats <- sumstats %>% filter(if_all(all_of(p_cols), ~ !is.na(.)))

z <- as.matrix(select(sumstats, all_of(informative_cols)))

mpi0 <- pi0_model(
  z = z,
  indep_snps = sumstats$prune_in
)

fpi0 <- fpi0est(
  sumstats$p_value,
  mpi0,
  indep_snps = sumstats$prune_in
)

sffdr_out <- sffdr(
  p.value = sumstats$p_value,
  fpi0 = fpi0$fpi0,
)

saveRDS(sffdr_out, file = snakemake@output$rds)

sumstats <- sumstats %>% mutate(
  fp = sffdr_out$fpvalues,
  fq = sffdr_out$fqvalues,
  flfdr = sffdr_out$flfdr,
  fpi0 = sffdr_out$fpi0
)

write_tsv(sumstats, file = snakemake@output$tsv)

theme_set(theme_minimal())

ggsave(
  plot(sffdr_out, rng = c(1e-15, 1)),
  file = snakemake@output$png,
  height = 8,
  width = 8
)
