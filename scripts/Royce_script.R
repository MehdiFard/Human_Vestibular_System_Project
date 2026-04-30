##Manhattan Plot
library(tidyverse)
library(ggtext)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/Owner/Desktop/Balance/")

UKB <- read_table("C:/Users/Owner/Desktop/Balance/UKB_balance_gwas/balance.imputed.stats/balance_UKB.tab.gz", 
               col_types = cols(CHISQ_LINREG = col_skip(), 
               P_LINREG = col_skip(), CHISQ_BOLT_LMM_INF = col_skip(), 
               P_BOLT_LMM_INF = col_skip(), CHISQ_BOLT_LMM = col_skip(),
               GENPOS = col_skip()))
colnames(UKB)[10] <- "P"

sig_data <- UKB |>
  subset(P < 0.05)
notsig_data <- UKB |>
  subset(P >= 0.05) |>
  group_by(CHR) |>
  sample_frac(0.05)

gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- 5e-08
snpsOfInterest <- c("rs2806949", 'rs73014202')
UKB_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), 
  color = as_factor(CHR), size = -log10(P) )) +
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") +
  geom_point(alpha = 0.50, size = 1.3) +
  geom_point(data = snpsOfInterest, color = 'orange', size = 2) +
  geom_label_repel(data=snpsOfInterest, aes(label=SNP), size=2) +
  scale_x_continuous(
    label = axis_set$CHR,
    breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(
    c("orchid", "lightskyblue"),
    unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = NULL, y = "-log<sub>10</sub>(P)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
     
 ggsave("C:/Users/Owner/Desktop/Balance/UKB.png", plot = UKB_plot, width = 6, height = 4,
        units = "in", dpi = 300)   
#######################################
 
gwas_data$highlight_snp[gwas_data$CHR == 19 & gwas_data$BP > 23468509 &
                 gwas_data$BP < 23727940] <- TRUE
  
UKB_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
   # Create the main plot
   geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = .4) +
   scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
   scale_x_continuous(
     label = unique(gwas_data$CHR),
     breaks = gwas_data %>%
       group_by(CHR) %>%
       summarize(center = mean(bp_cum)) %>%
       pull(center)) + 
  
   geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
   labs(
     x = "Chromosome",
     y = "-log10(p-value)"
   ) +
   theme_minimal() +
   theme(
     legend.position = "none", 
     panel.grid.major.x = element_blank(),
     panel.grid.minor.x = element_blank(), 
     axis.text.x = element_text(size = 6)
   ) +
  geom_point(data = subset(gwas_data, highlight_snp == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .4)
   

ggsave("C:/Users/Owner/Desktop/Balance/UKB.png", plot = UKB_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300)   
#########################################################################
###META-EUR
EUR_meta <- read_table("C:/Users/Owner/Desktop/Balance/balance_meta_EUR.tbl.sorted.filter.gz",  
col_types = cols(Allele1 = col_skip(), 
                 Allele2 = col_skip(), Freq1 = col_skip(), 
                 FreqSE = col_skip(), MinFreq = col_skip(),
                 MaxFreq = col_skip(), Effect = col_skip(), 
                 StdErr = col_skip(), Direction = col_skip()))               
colnames(EUR_meta)[4] <- "P"
colnames(EUR_meta)[1] <- "CHR"
colnames(EUR_meta)[2] <- 'BP'
colnames(EUR_meta)[3] <- 'SNP'

sig_data <- EUR_meta |>
  subset(P < 0.05)
notsig_data <- EUR_meta |>
  subset(P >= 0.05) |>
  group_by(CHR) |>
  sample_frac(0.05)

gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- 5e-08

gwas_data$highlight_snp <- ifelse(
  (gwas_data$CHR == 3 & gwas_data$BP > 131447659 & gwas_data$BP < 131510204) |
  (gwas_data$CHR == 17 & gwas_data$BP >  65763870 & gwas_data$BP < 66098155) |
  (gwas_data$CHR == 19 & gwas_data$BP >  23468509 & gwas_data$BP < 23741858),
                          TRUE, FALSE)

UKB_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
  # Create the main plot
  geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = .4) +
  scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
  scale_x_continuous(
    label = unique(gwas_data$CHR),
    breaks = gwas_data %>%
      group_by(CHR) %>%
      summarize(center = mean(bp_cum)) %>%
      pull(center)) + 
  
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
  labs(
    x = "Chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(size = 6)
  ) +
  geom_point(data = subset(gwas_data, highlight_snp == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .4)
ggsave("C:/Users/Owner/Desktop/Balance/meta_EUR.png", plot = UKB_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300)   

#########################################################################
###META
META <- read_table("Meta_analysis_balance/balance_meta1.tbl.sorted.filtered.gz", 
           col_types = cols(Allele1 = col_skip(), 
           Allele2 = col_skip(), Freq1 = col_skip(), 
           FreqSE = col_skip(), MinFreq = col_skip(), 
           MaxFreq = col_skip(), Effect = col_skip(), 
           StdErr = col_skip(), Direction = col_skip(), 
            HetISq = col_skip(), HetChiSq = col_skip()))   
META <- META[META$HetPVal > .05,1:4]
colnames(META)[4] <- "P"
colnames(META)[1] <- "CHR"
colnames(META)[2] <- 'BP'
colnames(META)[3] <- 'SNP'

sig_data <- META |>
  subset(P < 0.05)
notsig_data <- META |>
  subset(P >= 0.05) |>
  group_by(CHR) |>
  sample_frac(0.05)

gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- 5e-08

gwas_data$highlight_snp <- ifelse(
  (gwas_data$CHR == 17 & gwas_data$BP >  65814381 & gwas_data$BP < 66098155) |
    (gwas_data$CHR == 18 & gwas_data$BP > 50555930 & gwas_data$BP <	50641697) |
    (gwas_data$CHR == 19 & gwas_data$BP > 23468509 & gwas_data$BP < 23772486),
              TRUE, FALSE)

UKB_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
  # Create the main plot
  geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = .4) +
  scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
  scale_x_continuous(
    label = unique(gwas_data$CHR),
    breaks = gwas_data %>%
      group_by(CHR) %>%
      summarize(center = mean(bp_cum)) %>%
      pull(center)) + 
  
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
  labs(
    x = "Chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(size = 6)
  ) +
  geom_point(data = subset(gwas_data, highlight_snp == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .4)


ggsave("C:/Users/Owner/Desktop/Balance/META.png", plot = UKB_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300) 
############################################################
#########################     Genes  ##########################################
###UKB

UKB_genes <- read_excel("Meta_analysis_balance.xlsm", 
              sheet = "A3. UKB_magma genes", col_types = c("skip", 
              "numeric", "numeric", "numeric", 
              "numeric", "numeric", "skip", "skip", 
              "numeric", "text", "skip"))
colnames(UKB_genes)[2] <- 'BP'
colnames(UKB_genes)[7] <- 'GENE'

gwas_data <- UKB_genes %>% 
  arrange(CHR, BP)

data_cum <- gwas_data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- .05/19042

gwas_data$highlight_gene <- ifelse(gwas_data$P <= sig, TRUE, FALSE)
highlight_gene <- gwas_data$GENE[gwas_data$highlight_gene]
  

UKB_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
  # Create the main plot
  geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = .4) +
  scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
  scale_x_continuous(
    label = unique(gwas_data$CHR),
    breaks = gwas_data %>%
      group_by(CHR) %>%
      summarize(center = mean(bp_cum)) %>%
      pull(center)) + 
  geom_point(data = subset(gwas_data, highlight_gene == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .6)  +
  geom_text_repel(data = highlight_gene, aes(label = GENE),
                  box.padding = 0.5, point.padding = 0.5,
                  segment.color = 'black') +
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
  labs(
    x = "Chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(size = 6)
  ) 

ggsave("C:/Users/Owner/Desktop/Balance/UKB_genes.png", plot = UKB_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300) 
#######META_EUR
EUR_genes <- read_excel("Meta-analysis-balance.xlsm", 
                         sheet = "B3. EURmagma genes", col_types = c("skip", 
                  "numeric", "numeric", "numeric", 
                  "numeric", "numeric", "skip", "skip", 
                    "numeric", "text", "skip", "skip"))
colnames(EUR_genes)[2] <- 'BP'
colnames(EUR_genes)[7] <- 'GENE'
gwas_data <- EUR_genes %>% 
  arrange(CHR, BP)

data_cum <- gwas_data |>
  group_by(CHR) |>
  arrange(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- .05/18410

gwas_data$highlight_gene <- ifelse(gwas_data$P <= sig, TRUE, FALSE)
highlight_gene <- gwas_data$GENE[gwas_data$highlight_gene]

EUR_gene_plot <- ggplot(data = gwas_data, aes(x = bp_cum, y = -log10(P))) +
  # Create the main plot
  geom_point(gwas_data, aes(color = as.factor(CHR)
                            .), alpha = 0.5, size = .4) +
  scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
  scale_x_continuous(
    label = unique(CHR),
    breaks = gwas_data %>%
      group_by(CHR) %>%
      summarize(center = mean(bp_cum)) %>%
      pull(center)) + 
  geom_point(data = subset(gwas_data, highlight_gene == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .6)  +
  geom_text_repel(data = highlight_gene, aes(label = GENE),
                  box.padding = 0.5, point.padding = 0.5,
                  segment.color = 'black') +
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
  labs(
    x = "Chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(size = 6)
  ) 
#####################################################


ggsave("C:/Users/Owner/Desktop/Balance/EUR_genes.png", plot = EUR_gene_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300) 

#########META_genes
META_genes <- read_excel("Meta_analysis_balance.xlsm", 
                         sheet = "C3. meta_magma genes", col_types = c("skip", 
                          "numeric", "numeric", "numeric", 
                          "numeric", "numeric", "skip", "skip", 
                          "numeric", "text", "skip", "skip"))
colnames(META_genes)[2] <- 'BP'
colnames(META_genes)[7] <- 'GENE'

gwas_data <- META_genes %>%
  arrange(CHR, BP)

data_cum <- gwas_data |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(CHR, bp_add)

gwas_data <- gwas_data |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gwas_data |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- .05/18977

gwas_data$highlight_gene <- ifelse(gwas_data$P <= sig, TRUE, FALSE)
highlight_gene <- gwas_data$GENE[gwas_data$highlight_gene]

META_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
  # Create the main plot
  geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = .4) +
  scale_color_manual(values = rep(c("orchid", "skyblue1"), 22)) +
  scale_x_continuous(
    label = unique(gwas_data$CHR),
    breaks = gwas_data %>%
      group_by(CHR) %>%
      summarize(center = mean(bp_cum)) %>%
      pull(center)) + 
geom_point(data = subset(gwas_data, highlight_gene == TRUE),
             aes(x = bp_cum, y = -log10(P)),
             color = "brown", size = .6)  +
  geom_text_repel(data = highlight_gene, aes(label = GENE),
                  box.padding = 0.5, point.padding = 0.5,
                  segment.color = 'black') +
  geom_hline(yintercept = -log10(sig), color = "red",linetype = "dashed") + 
  labs(
    x = "Chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(size = 6)
  ) 

ggsave("C:/Users/Owner/Desktop/Balance/META_genes.png", plot = META_plot, 
       width = 6, height = 4,
       units = "in", dpi = 300) 


