############################################################
############################################################

rm(list = ls())
gc()

############################
## 0. Packages
############################
pkgs <- c(
  "readxl", "dplyr", "ggplot2", "patchwork",
  "stringr", "tidyr", "forcats", "tibble",
  "ggalluvial"
)
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(tibble)
  library(ggalluvial)
})

############################
## 1. Paths
############################
# Input and output directories
# This script reads the mediation result table from "results/mediation".
# Figure files will be saved in "figures/mediation".

xlsx_file_main <- file.path(
  "results", "mediation",
  "mediation result.xlsx"
)

outdir <- file.path("figures", "mediation")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

############################
## 2. Helper functions
############################
fmt_p_label <- function(x) {
  ifelse(
    is.na(x), "P = NA",
    ifelse(x < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", x)))
  )
}

fmt_beta_p_label <- function(beta, p) {
  paste0("Beta = ", sprintf("%.3f", beta), ", ", fmt_p_label(p))
}

fmt_prop_pct <- function(x, digits = 1) {
  ifelse(
    is.na(x), "NA",
    paste0(sprintf(paste0("%.", digits, "f"), x * 100), "%")
  )
}

pretty_trait <- function(x) {
  dplyr::case_when(
    x == "TailFat_g" ~ "TailFat",
    TRUE ~ x
  )
}

make_source_plotmath <- function(feature, block) {
  feature <- gsub("\\\\", "\\\\\\\\", feature)
  feature <- gsub("'", "\\\\'", feature, fixed = TRUE)
  paste0("italic('", feature, "')~'(", block, ")'")
}

save_pdf_editable <- function(plot_obj, filename, width, height, family = "Helvetica") {
  ggplot2::ggsave(
    filename = filename,
    plot = plot_obj,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    family = family
  )
}

base_family <- "Helvetica"

theme_clean <- theme_bw(base_family = base_family) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    legend.title = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  )

## Soft color palette
col_ile  <- "#D79A93"
col_rum  <- "#7AA6A1"
col_col  <- "#9AA9C9"
col_bar  <- "#B79C84"
col_low  <- "#F3E4D4"
col_high <- "#86AAA5"

## Triangle-plot colors based on the previous style
col_node  <- "#5A8991"
col_edge  <- "#606060"
col_arrow <- "#E32636"
col_guide <- "#D96B76"
col_text  <- "#222222"

############################
## 3. Read data
############################
sheet_names <- excel_sheets(xlsx_file_main)
if (!"combined_sig_only" %in% sheet_names) {
  stop("The worksheet combined_sig_only was not found in the Excel file.")
}

dat_main <- read_excel(xlsx_file_main, sheet = "combined_sig_only")

if (nrow(dat_main) == 0) stop("The combined_sig_only worksheet is empty; figures cannot be generated.")

req_cols <- c(
  "source_block", "source_feature", "mediator_feature", "trait_feature",
  "a_est", "a_p", "b_est", "b_p", "cprime_est", "cprime_p", "prop_med"
)
miss_cols <- setdiff(req_cols, colnames(dat_main))
if (length(miss_cols) > 0) {
  stop("The result table is missing required columns: ", paste(miss_cols, collapse = ", "))
}

if ("prop_med_p" %in% colnames(dat_main)) {
  dat_main <- dat_main %>% mutate(center_p = prop_med_p)
} else if ("acme_p" %in% colnames(dat_main)) {
  dat_main <- dat_main %>% mutate(center_p = acme_p)
} else {
  stop("Neither prop_med_p nor acme_p was found in the result table.")
}

############################
## 4. Figure 1
############################
## Panel A: number of significant chains by compartment and phenotype
fig1_a_dat <- dat_main %>%
  mutate(
    trait_feature = pretty_trait(trait_feature),
    source_block = factor(source_block, levels = c("Rum", "Ile", "Col"))
  ) %>%
  count(source_block, trait_feature, name = "n")

p_fig1_a <- ggplot(fig1_a_dat, aes(x = trait_feature, y = n, fill = source_block)) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.4,
    position = position_dodge(width = 0.78)
  ) +
  scale_fill_manual(values = c("Rum" = col_rum, "Ile" = col_ile, "Col" = col_col)) +
  labs(x = NULL, y = "Significant chains", fill = NULL) +
  theme_clean +
  theme(legend.position = "top")

## Panel B: total number of chains for core mediator genes
fig1_b_dat <- dat_main %>%
  count(mediator_feature, name = "n") %>%
  arrange(desc(n)) %>%
  slice_head(n = 10) %>%
  mutate(mediator_feature = fct_reorder(mediator_feature, n))

p_fig1_b <- ggplot(fig1_b_dat, aes(x = n, y = mediator_feature)) +
  geom_col(width = 0.72, fill = col_bar, color = "black", linewidth = 0.4) +
  labs(x = "Significant chains", y = NULL) +
  theme_clean

## Panel C: heatmap of core mediator genes by phenotype
fig1_c_dat <- dat_main %>%
  mutate(trait_feature = pretty_trait(trait_feature)) %>%
  filter(mediator_feature %in% c("SCARB1", "HSD17B12", "APOB", "AMACR", "HADHB", "PCK2", "FDFT1")) %>%
  count(mediator_feature, trait_feature, name = "n") %>%
  complete(mediator_feature, trait_feature, fill = list(n = 0))

p_fig1_c <- ggplot(fig1_c_dat, aes(x = trait_feature, y = mediator_feature, fill = n)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = n), size = 3.8, family = base_family) +
  scale_fill_gradient(low = col_low, high = col_high) +
  labs(x = NULL, y = NULL, fill = "Count") +
  theme_clean +
  theme(legend.position = "right")

p_fig1 <- p_fig1_a + p_fig1_b + p_fig1_c +
  plot_layout(widths = c(1.15, 1.00, 1.10))

############################
## 5. Selected 9 chains
############################
## Final selected set of 9 chains
targets_main9 <- tibble::tribble(
  ~source_block, ~source_feature,                         ~mediator_feature, ~trait_feature, ~panel_order,
  "Rum",         "Alcaligenes",                           "HSD17B12",        "TBA",          1,
  "Rum",         "Sharpea",                               "SCARB1",          "TBA",          2,
  "Rum",         "Prevotella_7",                          "HSD17B12",        "TBA",          3,
  "Col",         "Fournierella",                          "SCARB1",          "TBA",          4,
  "Col",         "GWE2_31_10",                            "HSD17B12",        "TBA",          5,
  "Col",         "unclassified_Hydrogenoanaerobacterium", "HSD17B12",        "TBA",          6,
  "Rum",         "Klebsiella",                            "AMACR",           "TG",           7,
  "Rum",         "Anaerovibrio",                          "HADHB",           "TailFat_g",    8,
  "Col",         "[Eubacterium]_oxidoreducens_group",     "PCK2",            "GLU",          9
)

dat_plot_main9 <- dat_main %>%
  inner_join(
    targets_main9,
    by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
  ) %>%
  arrange(panel_order) %>%
  mutate(
    panel_id      = seq_len(n()),
    source_lab    = paste0(source_feature, " (", source_block, ")"),
    mediator_lab  = mediator_feature,
    trait_lab     = pretty_trait(trait_feature),
    source_expr   = vapply(seq_len(n()), function(i) {
      make_source_plotmath(source_feature[i], source_block[i])
    }, character(1)),
    lab_a         = fmt_beta_p_label(a_est, a_p),
    lab_b         = fmt_beta_p_label(b_est, b_p),
    lab_c         = fmt_beta_p_label(cprime_est, cprime_p),
    lab_m         = paste0(fmt_prop_pct(prop_med, digits = 1), "\n", fmt_p_label(center_p))
  )

if (nrow(dat_plot_main9) != 9) {
  miss_show <- targets_main9 %>%
    anti_join(
      dat_main,
      by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
    )
  print(miss_show)
  stop("The selected 9 chains were not fully matched to the result table. Please check whether the names are consistent.")
}

############################
## 6. Triangle plotting function
############################
make_triangle_plot <- function(dat_plot, ncol_wrap = 3, source_parse = TRUE, mediator_parse = FALSE) {
  
  ## Geometric parameters for the equilateral triangle based on the previous mediation script
  side_len <- 0.50
  
  x_left  <- 0.25
  y_left  <- 0.19
  
  x_right <- x_left + side_len
  y_right <- y_left
  
  x_top <- (x_left + x_right) / 2
  y_top <- y_left + sqrt(3) / 2 * side_len
  
  centroid_x <- (x_left + x_top + x_right) / 3
  centroid_y <- (y_left + y_top + y_right) / 3
  
  get_inward_normal <- function(x1, y1, x2, y2, cx, cy) {
    vx <- x2 - x1
    vy <- y2 - y1
    vlen <- sqrt(vx^2 + vy^2)
    
    nx1 <- -vy / vlen
    ny1 <-  vx / vlen
    nx2 <-  vy / vlen
    ny2 <- -vx / vlen
    
    mx <- (x1 + x2) / 2
    my <- (y1 + y2) / 2
    
    dot1 <- (cx - mx) * nx1 + (cy - my) * ny1
    if (dot1 > 0) {
      c(nx1, ny1, vlen, vx / vlen, vy / vlen)
    } else {
      c(nx2, ny2, vlen, vx / vlen, vy / vlen)
    }
  }
  
  left_geo  <- get_inward_normal(x_left, y_left, x_top, y_top, centroid_x, centroid_y)
  nx_left   <- left_geo[1]
  ny_left   <- left_geo[2]
  ux_left   <- left_geo[4]
  uy_left   <- left_geo[5]
  
  right_geo <- get_inward_normal(x_top, y_top, x_right, y_right, centroid_x, centroid_y)
  nx_right  <- right_geo[1]
  ny_right  <- right_geo[2]
  ux_right  <- right_geo[4]
  uy_right  <- right_geo[5]
  
  offset_left  <- 0.024
  offset_right <- 0.028
  
  left_start_shrink <- 0.030
  right_end_shrink  <- 0.055
  
  joint_drop <- 0.055
  joint_x <- x_top
  joint_y <- y_top - joint_drop
  
  edge_df <- bind_rows(
    dat_plot %>% transmute(panel_id, x = x_left, y = y_left, xend = x_top,   yend = y_top),
    dat_plot %>% transmute(panel_id, x = x_top,  y = y_top,  xend = x_right, yend = y_right),
    dat_plot %>% transmute(panel_id, x = x_left, y = y_left, xend = x_right, yend = y_right)
  )
  
  guide_df <- dat_plot %>%
    transmute(
      panel_id,
      x = joint_x,
      y = joint_y,
      xend = x_left + nx_left * offset_left + ux_left * left_start_shrink,
      yend = y_left + ny_left * offset_left + uy_left * left_start_shrink
    )
  
  arrow_df <- dat_plot %>%
    transmute(
      panel_id,
      x = joint_x,
      y = joint_y,
      xend = x_right + nx_right * offset_right - ux_right * right_end_shrink,
      yend = y_right + ny_right * offset_right - uy_right * right_end_shrink
    )
  
  node_df <- bind_rows(
    dat_plot %>% transmute(panel_id, x = x_left,  y = y_left,  label = source_expr,  type = "source"),
    dat_plot %>% transmute(panel_id, x = x_top,   y = y_top,   label = mediator_lab, type = "mediator"),
    dat_plot %>% transmute(panel_id, x = x_right, y = y_right, label = trait_lab,    type = "trait")
  )
  
  ann_df <- dat_plot %>% transmute(panel_id, lab_a, lab_b, lab_c, lab_m)
  
  ann_a <- ann_df %>% mutate(
    x = (x_left + x_top) / 2 - 0.020,
    y = (y_left + y_top) / 2 + 0.004,
    angle = 60
  )
  
  ann_b <- ann_df %>% mutate(
    x = (x_top + x_right) / 2 + 0.020,
    y = (y_top + y_right) / 2 + 0.004,
    angle = -60
  )
  
  ann_c <- ann_df %>% mutate(
    x = (x_left + x_right) / 2,
    y = y_left + 0.038,
    angle = 0
  )
  
  ann_m <- ann_df %>% mutate(
    x = x_top,
    y = y_left + (y_top - y_left) * 0.43
  )
  
  theme_tri <- theme_void(base_family = base_family) +
    theme(
      strip.text       = element_blank(),
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.spacing.x  = unit(1.4, "lines"),
      panel.spacing.y  = unit(1.6, "lines"),
      plot.margin      = margin(16, 18, 18, 24)
    )
  
  p <- ggplot() +
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.90,
      color = col_edge,
      lineend = "round"
    ) +
    geom_segment(
      data = guide_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.72,
      color = col_guide,
      lineend = "round"
    ) +
    geom_segment(
      data = arrow_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.92,
      color = col_arrow,
      arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
      lineend = "round"
    ) +
    geom_point(
      data = node_df,
      aes(x = x, y = y),
      size = 9.8,
      shape = 16,
      color = col_node
    ) +
    geom_text(
      data = subset(node_df, type == "source"),
      aes(x = x + 0.010, y = y - 0.118, label = label),
      parse = TRUE,
      family = base_family,
      hjust = 0.5,
      size = 4.8,
      color = col_text
    ) +
    geom_text(
      data = subset(node_df, type == "mediator"),
      aes(x = x, y = y + 0.090, label = label),
      family = base_family,
      size = 4.7,
      color = col_text
    ) +
    geom_text(
      data = subset(node_df, type == "trait"),
      aes(x = x, y = y - 0.118, label = label),
      family = base_family,
      size = 4.7,
      color = col_text
    ) +
    geom_text(
      data = ann_a,
      aes(x = x, y = y, label = lab_a, angle = angle),
      family = base_family,
      size = 3.7,
      color = col_text
    ) +
    geom_text(
      data = ann_b,
      aes(x = x, y = y, label = lab_b, angle = angle),
      family = base_family,
      size = 3.7,
      color = col_text
    ) +
    geom_text(
      data = ann_c,
      aes(x = x, y = y, label = lab_c, angle = angle),
      family = base_family,
      size = 3.7,
      color = col_text
    ) +
    geom_text(
      data = ann_m,
      aes(x = x, y = y, label = lab_m),
      family = base_family,
      fontface = "bold",
      lineheight = 0.95,
      size = 4.1,
      color = col_text
    ) +
    coord_cartesian(
      xlim = c(0.03, 0.97),
      ylim = c(0.03, 0.97),
      clip = "off"
    ) +
    facet_wrap(~ panel_id, ncol = ncol_wrap) +
    theme_tri
  
  return(p)
}

############################
## 7. Sankey plotting function
############################
make_mediation_sankey <- function(dat_sankey) {
  
  sankey_df <- dat_sankey %>%
    transmute(
      Source_show   = paste0(source_feature, " (", source_block, ")"),
      Mediator_show = mediator_feature,
      Trait_show    = pretty_trait(trait_feature),
      Direction_src = factor(source_block, levels = c("Rum", "Ile", "Col")),
      weight        = prop_med
    )
  
  source_order <- sankey_df %>%
    group_by(Direction_src, Source_show) %>%
    summarise(w = sum(weight), .groups = "drop") %>%
    arrange(Direction_src, desc(w), Source_show) %>%
    pull(Source_show) %>%
    unique()
  
  mediator_order <- sankey_df %>%
    group_by(Mediator_show) %>%
    summarise(w = sum(weight), .groups = "drop") %>%
    arrange(desc(w), Mediator_show) %>%
    pull(Mediator_show)
  
  trait_order <- c("TBA", "TG", "GLU", "TailFat")
  trait_order <- trait_order[trait_order %in% unique(sankey_df$Trait_show)]
  
  sankey_df <- sankey_df %>%
    mutate(
      Source_show   = factor(Source_show, levels = source_order),
      Mediator_show = factor(Mediator_show, levels = mediator_order),
      Trait_show    = factor(Trait_show, levels = trait_order)
    )
  
  lodes_df <- ggalluvial::to_lodes_form(
    data = sankey_df,
    axes = c("Source_show", "Mediator_show", "Trait_show"),
    key = "x",
    value = "stratum",
    id = "alluvium"
  )
  
  lodes_df$x <- factor(
    lodes_df$x,
    levels = c("Source_show", "Mediator_show", "Trait_show"),
    labels = c("Source", "Mediator", "Trait")
  )
  
  fill_map <- c(
    "Rum" = col_rum,
    "Ile" = col_ile,
    "Col" = col_col
  )
  
  p_sankey <- ggplot(
    sankey_df,
    aes(
      axis1 = Source_show,
      axis2 = Mediator_show,
      axis3 = Trait_show,
      y = weight
    )
  ) +
    geom_alluvium(
      aes(fill = Direction_src),
      width = 0.13,
      knot.pos = 0.42,
      alpha = 0.84,
      color = "white",
      linewidth = 0.25
    ) +
    geom_stratum(
      data = lodes_df %>% filter(x == "Source"),
      aes(x = x, stratum = stratum, y = weight),
      inherit.aes = FALSE,
      width = 0.56,
      fill = "#EEE8DF",
      color = "grey45",
      linewidth = 0.42
    ) +
    geom_stratum(
      data = lodes_df %>% filter(x == "Mediator"),
      aes(x = x, stratum = stratum, y = weight),
      inherit.aes = FALSE,
      width = 0.20,
      fill = "#EEE8DF",
      color = "grey45",
      linewidth = 0.42
    ) +
    geom_stratum(
      data = lodes_df %>% filter(x == "Trait"),
      aes(x = x, stratum = stratum, y = weight),
      inherit.aes = FALSE,
      width = 0.44,
      fill = "#EEE8DF",
      color = "grey45",
      linewidth = 0.42
    ) +
    geom_text(
      data = lodes_df %>% filter(x == "Source"),
      stat = "stratum",
      aes(x = x, stratum = stratum, y = weight, label = after_stat(stratum)),
      inherit.aes = FALSE,
      family = base_family,
      size = 3.3,
      fontface = "bold",
      lineheight = 0.92
    ) +
    geom_text(
      data = lodes_df %>% filter(x == "Mediator"),
      stat = "stratum",
      aes(x = x, stratum = stratum, y = weight, label = after_stat(stratum)),
      inherit.aes = FALSE,
      family = base_family,
      size = 3.2,
      fontface = "bold",
      lineheight = 0.92
    ) +
    geom_text(
      data = lodes_df %>% filter(x == "Trait"),
      stat = "stratum",
      aes(x = x, stratum = stratum, y = weight, label = after_stat(stratum)),
      inherit.aes = FALSE,
      family = base_family,
      size = 3.4,
      fontface = "bold",
      lineheight = 0.92
    ) +
    scale_fill_manual(values = fill_map, name = NULL) +
    scale_x_discrete(
      limits = c("Source", "Mediator", "Trait"),
      expand = c(0.06, 0.06)
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_family = base_family) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 11, face = "bold", color = "black"),
      legend.position = "top",
      legend.text = element_text(size = 10, color = "black"),
      plot.margin = margin(14, 18, 14, 18)
    )
  
  return(p_sankey)
}

############################
## 8. Build plots
############################
p_fig2 <- make_triangle_plot(
  dat_plot = dat_plot_main9,
  ncol_wrap = 3,
  source_parse = TRUE,
  mediator_parse = FALSE
)

p_fig3 <- make_mediation_sankey(dat_plot_main9)

############################
## 9. Supplementary Sankey S1
############################
## Compact overall version: 12 chains (TBA 6 + TG 2 + GLU 2 + TailFat 1)
targets_sankey_s1 <- tibble::tribble(
  ~source_block, ~source_feature,                         ~mediator_feature, ~trait_feature, ~group_type, ~group_order,
  "Rum",         "Alcaligenes",                           "HSD17B12",        "TBA",          "TBA",       1,
  "Rum",         "Sharpea",                               "SCARB1",          "TBA",          "TBA",       2,
  "Rum",         "Prevotella_7",                          "HSD17B12",        "TBA",          "TBA",       3,
  "Col",         "Fournierella",                          "SCARB1",          "TBA",          "TBA",       4,
  "Col",         "GWE2_31_10",                            "HSD17B12",        "TBA",          "TBA",       5,
  "Col",         "unclassified_Hydrogenoanaerobacterium", "HSD17B12",        "TBA",          "TBA",       6,
  "Rum",         "Klebsiella",                            "AMACR",           "TG",           "TG",        7,
  "Rum",         "Saccharofermentans",                    "AMACR",           "TG",           "TG",        8,
  "Col",         "[Eubacterium]_oxidoreducens_group",     "PCK2",            "GLU",          "GLU",       9,
  "Col",         "Lachnospiraceae_UCG_010",               "FDFT1",           "GLU",          "GLU",       10,
  "Rum",         "Anaerovibrio",                          "HADHB",           "TailFat_g",    "TailFat",   11
)

dat_sankey_s1 <- dat_main %>%
  inner_join(
    targets_sankey_s1,
    by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
  ) %>%
  arrange(group_order)

if (nrow(dat_sankey_s1) != nrow(targets_sankey_s1)) {
  miss_s1 <- targets_sankey_s1 %>%
    anti_join(
      dat_main,
      by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
    )
  print(miss_s1)
  stop("Some chains in the S1 supplementary Sankey plot were not matched to the result table. Please check whether the names are consistent.")
}

p_figS1 <- make_mediation_sankey(dat_sankey_s1)

############################
## 10. Supplementary Sankey S2
############################
## TBA-focused version: 8 chains (Rum 4 + Col 4)
targets_sankey_s2 <- tibble::tribble(
  ~source_block, ~source_feature,                         ~mediator_feature, ~trait_feature, ~group_order,
  "Rum",         "Alcaligenes",                           "HSD17B12",        "TBA",          1,
  "Rum",         "Sharpea",                               "SCARB1",          "TBA",          2,
  "Rum",         "Prevotella_7",                          "HSD17B12",        "TBA",          3,
  "Rum",         "Olsenella",                             "HSD17B12",        "TBA",          4,
  "Col",         "Fournierella",                          "SCARB1",          "TBA",          5,
  "Col",         "GWE2_31_10",                            "HSD17B12",        "TBA",          6,
  "Col",         "unclassified_Hydrogenoanaerobacterium", "HSD17B12",        "TBA",          7,
  "Col",         "Flavonifractor",                        "SCARB1",          "TBA",          8
)

dat_sankey_s2 <- dat_main %>%
  inner_join(
    targets_sankey_s2,
    by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
  ) %>%
  arrange(group_order)

if (nrow(dat_sankey_s2) != nrow(targets_sankey_s2)) {
  miss_s2 <- targets_sankey_s2 %>%
    anti_join(
      dat_main,
      by = c("source_block", "source_feature", "mediator_feature", "trait_feature")
    )
  print(miss_s2)
  stop("Some chains in the S2 supplementary Sankey plot were not matched to the result table. Please check whether the names are consistent.")
}

p_figS2 <- make_mediation_sankey(dat_sankey_s2)

############################
## 11. Save PDFs only
############################
fig1_pdf   <- file.path(outdir, "Fig1_mediation_overall_distribution.pdf")
fig2_pdf   <- file.path(outdir, "Fig2_mediation_triangles_main9_equilateral_3col.pdf")
fig3_pdf   <- file.path(outdir, "Fig3_mediation_sankey_main9.pdf")
figS1_pdf  <- file.path(outdir, "FigS1_mediation_sankey_overall_compressed.pdf")
figS2_pdf  <- file.path(outdir, "FigS2_mediation_sankey_TBA_focused.pdf")

csv_main9  <- file.path(outdir, "Fig2_Fig3_main9_used_table.csv")
csv_s1     <- file.path(outdir, "FigS1_sankey_overall_compressed_used_table.csv")
csv_s2     <- file.path(outdir, "FigS2_sankey_TBA_focused_used_table.csv")

save_pdf_editable(
  p_fig1,
  fig1_pdf,
  width = 14.5,
  height = 5.8,
  family = base_family
)

save_pdf_editable(
  p_fig2,
  fig2_pdf,
  width = 15.5,
  height = 15.0,
  family = base_family
)

save_pdf_editable(
  p_fig3,
  fig3_pdf,
  width = 13.5,
  height = 8.8,
  family = base_family
)

save_pdf_editable(
  p_figS1,
  figS1_pdf,
  width = 14.2,
  height = 9.0,
  family = base_family
)

save_pdf_editable(
  p_figS2,
  figS2_pdf,
  width = 13.8,
  height = 8.8,
  family = base_family
)

write.csv(dat_plot_main9,  csv_main9, row.names = FALSE)
write.csv(dat_sankey_s1,   csv_s1,    row.names = FALSE)
write.csv(dat_sankey_s2,   csv_s2,    row.names = FALSE)

############################
## 12. Export supporting tables
############################
write.csv(fig1_a_dat, file.path(outdir, "Fig1_panelA_trait_by_source_count.csv"), row.names = FALSE)
write.csv(fig1_b_dat, file.path(outdir, "Fig1_panelB_top_mediator_count.csv"), row.names = FALSE)
write.csv(fig1_c_dat, file.path(outdir, "Fig1_panelC_key_mediator_trait_matrix.csv"), row.names = FALSE)

cat("Figure generation completed:\n")
cat("Fig1 PDF：", fig1_pdf, "\n")
cat("Fig2 PDF：", fig2_pdf, "\n")
cat("Fig3 PDF：", fig3_pdf, "\n")
cat("FigS1 PDF：", figS1_pdf, "\n")
cat("FigS2 PDF：", figS2_pdf, "\n")
cat("Output directory: ", outdir, "\n")