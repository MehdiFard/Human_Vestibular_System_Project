## ------------------------------------------------------------------
## Figure 1 — Cohort demographics
## Each panel: 8 cm x 8 cm (square, identical canvas)
## Font 7 pt; donut % labels sit at wedge radial midpoint
## ------------------------------------------------------------------
rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

OUT_DIR <- "figures/demographics/"
dir.create(OUT_DIR, showWarnings = FALSE)

SIDE_MM <- 80

save_panel <- function(p, name) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = SIDE_MM, height = SIDE_MM,
         units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = SIDE_MM, height = SIDE_MM,
         units = "mm", dpi = 600)
}

## ---- Data ---------------------------------------------------------
demo <- tibble::tribble(
  ~cohort,  ~ancestry, ~sex,     ~cases, ~controls,
  "UKB",    "EUR",     "Male",   17364,  147507,
  "UKB",    "EUR",     "Female", 28814,  170349,
  "MVP",    "EUR",     "Male",   33605,  247696,
  "MVP",    "EUR",     "Female",  2968,   19474,
  "MVP",    "AFR",     "Male",    8801,   59280,
  "MVP",    "AFR",     "Female",  1334,    9973,
  "MVP",    "HIS",     "Male",    3282,   27304,
  "MVP",    "HIS",     "Female",   349,    3173
) |>
  mutate(
    total    = cases + controls,
    case_pct = 100 * cases / total,
    stratum  = paste(cohort, ancestry, sep = "-")
  )

stratum_levels <- demo |>
  group_by(stratum) |>
  summarise(N = sum(total), .groups = "drop") |>
  arrange(desc(N)) |>
  pull(stratum)

demo <- demo |>
  mutate(stratum = factor(stratum, levels = stratum_levels),
         sex     = factor(sex, levels = c("Male", "Female")))

## ---- Theme & palette ---------------------------------------------
BASE_SIZE <- 7

nature_theme <- theme_classic(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    axis.title          = element_text(size = BASE_SIZE),
    axis.title.x        = element_text(margin = margin(t = 3)),
    axis.text           = element_text(size = BASE_SIZE, colour = "black"),
    axis.line           = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks          = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length   = unit(0.08, "cm"),
    plot.title          = element_blank(),
    legend.title        = element_text(size = BASE_SIZE),
    legend.text         = element_text(size = BASE_SIZE),
    legend.key.size     = unit(0.28, "cm"),
    legend.margin       = margin(0, 0, 0, 0),
    legend.box.spacing  = unit(0.05, "cm"),
    legend.position     = "top",
    legend.justification = "left",
    plot.margin         = margin(4, 6, 4, 4)
  )

col_case    <- "#C44536"
col_control <- "#9CADBF"
col_male    <- "#3B6E8F"
col_female  <- "#E08E45"

strata_pal <- c(
  "UKB-EUR" = "#1F4E79",
  "MVP-EUR" = "#5B9BD5",
  "MVP-AFR" = "#C44536",
  "MVP-HIS" = "#E0A458"
)

## ====================================================================
## Panel A — Donut (8 x 8 cm)
## ====================================================================
## Approach: build the donut from polygon coordinates so we can place
## labels exactly at the radial midpoint of each wedge. This avoids
## ggplot's coord_polar weirdness with text positioning.
##
## Donut geometry (in arbitrary plot units):
##   inner radius = 1.0
##   outer radius = 1.8
##   label radius = 1.4 (radial midpoint of the wedge band)
## ------------------------------------------------------------------

R_IN  <- 1.0
R_OUT <- 1.8
R_LAB <- (R_IN + R_OUT) / 2

donut_df <- demo |>
  group_by(stratum) |>
  summarise(N = sum(total), .groups = "drop") |>
  mutate(
    stratum = factor(stratum, levels = stratum_levels),
    pct     = 100 * N / sum(N),
    pct_lab = sprintf("%.1f%%", pct)
  ) |>
  arrange(stratum) |>
  mutate(
    frac      = N / sum(N),
    end_ang   = cumsum(frac) * 2 * pi,
    start_ang = lag(end_ang, default = 0),
    mid_ang   = (start_ang + end_ang) / 2,
    ## label x/y (Cartesian; angle measured from 12 o'clock, clockwise,
    ## to match coord_polar with start = 0)
    lab_x = R_LAB * sin(mid_ang),
    lab_y = R_LAB * cos(mid_ang)
  )

## Build wedge polygons explicitly
wedge_poly <- function(start_ang, end_ang, r_in, r_out, n = 60) {
  a <- seq(start_ang, end_ang, length.out = n)
  data.frame(
    x = c(r_out * sin(a), rev(r_in * sin(a))),
    y = c(r_out * cos(a), rev(r_in * cos(a)))
  )
}

poly_df <- donut_df |>
  rowwise() |>
  mutate(poly = list(wedge_poly(start_ang, end_ang, R_IN, R_OUT))) |>
  tidyr::unnest(poly) |>
  ungroup()

## Choose label colour: white inside wedge if good contrast, dark for tan
## All four palette colors are dark enough that white text reads well.
donut_df <- donut_df |>
  mutate(label_col = ifelse(stratum == "MVP-HIS", "grey15", "white"))

pA <- ggplot() +
  geom_polygon(data = poly_df,
               aes(x = x, y = y, group = stratum, fill = stratum),
               colour = "white", linewidth = 0.5) +
  geom_text(data = donut_df,
            aes(x = lab_x, y = lab_y, label = pct_lab,
                colour = label_col),
            size = BASE_SIZE / .pt, fontface = "bold",
            show.legend = FALSE) +
  annotate("text", x = 0, y = 0.22,
           label = sprintf("n = %s", scales::comma(sum(demo$total))),
           size = BASE_SIZE / .pt, fontface = "bold") +
  annotate("text", x = 0, y = -0.05,
           label = sprintf("Cases: %s (%.0f%%)",
                           scales::comma(sum(demo$cases)),
                           100 * sum(demo$cases) / sum(demo$total)),
           size = 5.5 / .pt) +
  annotate("text", x = 0, y = -0.25,
           label = sprintf("Controls: %s (%.0f%%)",
                           scales::comma(sum(demo$controls)),
                           100 * sum(demo$controls) / sum(demo$total)),
           size = 5.5 / .pt) +
  scale_fill_manual(values = strata_pal, name = NULL) +
  scale_colour_identity() +
  coord_fixed(xlim = c(-2.0, 2.0), ylim = c(-2.0, 2.0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(colour = NA))) +
  theme_void(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    legend.position     = "bottom",
    legend.title        = element_blank(),
    legend.text         = element_text(size = BASE_SIZE),
    legend.key.size     = unit(0.32, "cm"),
    legend.margin       = margin(0, 0, 0, 0),
    legend.box.spacing  = unit(0.15, "cm"),
    plot.margin         = margin(4, 4, 4, 4)
  )
pA
save_panel(pA, "panel_A_donut")

## ====================================================================
## Panel B — Cases vs controls per stratum (8 x 8 cm)
## ====================================================================
plotB_df <- demo |>
  pivot_longer(c(cases, controls), names_to = "status", values_to = "n") |>
  mutate(status  = factor(ifelse(status == "cases", "Cases", "Controls"),
                          levels = c("Controls", "Cases")),
         row_lab = paste0(stratum, "  ", as.character(sex)))

row_order_B <- demo |>
  arrange(desc(stratum), desc(sex)) |>
  mutate(row_lab = paste0(stratum, "  ", as.character(sex))) |>
  pull(row_lab)

plotB_df <- plotB_df |>
  mutate(row_lab = factor(row_lab, levels = row_order_B))

totalsB <- demo |>
  mutate(row_lab = factor(paste0(stratum, "  ", as.character(sex)),
                          levels = row_order_B))

pB <- ggplot(plotB_df, aes(x = n, y = row_lab, fill = status)) +
  geom_col(width = 0.7, colour = NA) +
  geom_text(data = totalsB,
            aes(x = total, y = row_lab,
                label = scales::comma(total)),
            inherit.aes = FALSE,
            hjust = -0.12, size = BASE_SIZE / .pt, colour = "grey20") +
  scale_fill_manual(values = c("Controls" = col_control,
                               "Cases"    = col_case),
                    breaks = c("Cases", "Controls"),
                    name = NULL) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-3,
                                                   suffix = "K",
                                                   accuracy = 1),
                     limits = c(0, 360000),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Participants (n)", y = NULL) +
  nature_theme +
  theme(panel.grid.major.x = element_line(colour = "grey92",
                                          linewidth = 0.25),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(family = "mono", size = BASE_SIZE))
pB
save_panel(pB, "panel_B_stratum_n")

## ====================================================================
## Panel C — Sex composition per stratum (8 x 8 cm)
## ====================================================================
sex_df <- demo |>
  group_by(stratum, sex) |>
  summarise(n = sum(total), .groups = "drop") |>
  group_by(stratum) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup() |>
  mutate(sex = factor(sex, levels = c("Female", "Male")))

pC <- ggplot(sex_df, aes(x = pct, y = stratum, fill = sex)) +
  geom_col(width = 0.65, colour = NA) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = BASE_SIZE / .pt, fontface = "bold", colour = "white") +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("Female" = col_female,
                               "Male"   = col_male),
                    breaks = c("Male", "Female"),
                    name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Sex composition", y = NULL) +
  nature_theme +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank())
pC
save_panel(pC, "panel_C_sex_balance")

cat("done — three panels at 8 x 8 cm each\n")

## ====================================================================
## Panel D — Age donut (8 x 8 cm)
## Decorative donut with 4 equal quarters (one per cohort × status).
## ====================================================================

## ---- Source data --------------------------------------------------
## Order around the donut, starting at 12 o'clock, going clockwise:
##   1. UKB Cases         (top-right)
##   2. MVP Cases         (bottom-right)
##   3. MVP Controls      (bottom-left)
##   4. UKB Controls      (top-left)
age_df <- tibble::tribble(
  ~order, ~cohort, ~status,    ~mean_age, ~sd_age,
  1,      "UKB",   "Cases",     56.9,      NA_real_,
  2,      "MVP",   "Cases",     67.0,      11.73,
  3,      "MVP",   "Controls",  61.0,      14.62,
  4,      "UKB",   "Controls",  56.6,      NA_real_
) |>
  mutate(
    status  = factor(status, levels = c("Cases", "Controls")),
    age_lab = sprintf("%.1f y", mean_age),
    sd_lab  = ifelse(is.na(sd_age), NA_character_,
                     sprintf("(SD %.2f)", sd_age))
  )

## ---- Geometry ----------------------------------------------------
R_IN  <- 0.85
R_OUT <- 1.75
R_LAB <- (R_IN + R_OUT) / 2

n_wedges <- nrow(age_df)
age_df <- age_df |>
  mutate(
    start_ang = (order - 1) / n_wedges * 2 * pi,
    end_ang   = order        / n_wedges * 2 * pi,
    mid_ang   = (start_ang + end_ang) / 2,
    lab_x = R_LAB * sin(mid_ang),
    lab_y = R_LAB * cos(mid_ang),
    lab_y_age = ifelse(is.na(sd_lab), lab_y, lab_y + 0.07),
    lab_y_sd  = lab_y - 0.09
  )

poly_df <- age_df |>
  rowwise() |>
  mutate(verts = list(wedge_poly(start_ang, end_ang, R_IN, R_OUT))) |>
  unnest(verts) |>
  ungroup() |>
  mutate(group_id = paste(cohort, status))

## Cohort tags above/below the donut
cohort_tags <- tibble::tibble(
  label = c("UKB", "MVP"),
  x     = c(0, 0),
  y     = c(R_OUT + 0.18, -(R_OUT + 0.18))
)

## ---- Plot --------------------------------------------------------
pD <- ggplot() +
  geom_polygon(data = poly_df,
               aes(x = x, y = y, group = group_id, fill = status),
               colour = "white", linewidth = 0.6) +
  geom_text(data = age_df,
            aes(x = lab_x, y = lab_y_age, label = age_lab),
            size = BASE_SIZE / .pt, fontface = "bold",
            colour = "white") +
  geom_text(data = dplyr::filter(age_df, !is.na(sd_lab)),
            aes(x = lab_x, y = lab_y_sd, label = sd_lab),
            size = 5 / .pt, colour = "white") +
  geom_text(data = cohort_tags,
            aes(x = x, y = y, label = label),
            size = (BASE_SIZE + 1) / .pt, fontface = "bold",
            colour = "grey20") +
  scale_fill_manual(values = c("Cases"    = col_case,
                               "Controls" = col_control),
                    breaks = c("Cases", "Controls"),
                    name = NULL) +
  coord_fixed(xlim = c(-2.05, 2.05), ylim = c(-2.05, 2.05)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_void(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    legend.position     = "bottom",
    legend.title        = element_blank(),
    legend.text         = element_text(size = BASE_SIZE),
    legend.key.size     = unit(0.32, "cm"),
    legend.margin       = margin(0, 0, 0, 0),
    legend.box.spacing  = unit(0.15, "cm"),
    plot.margin         = margin(4, 4, 4, 4)
  )
pD
save_panel(pD, "panel_D_age_donut")
cat("done — panel D redesigned (4 equal quarters)\n")
