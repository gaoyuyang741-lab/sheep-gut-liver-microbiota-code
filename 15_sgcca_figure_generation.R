############################################################
############################################################

rm(list = ls())
gc()

############################
## 0. Packages
############################
pkgs <- c(
  "tidyverse", "readxl", "patchwork", "stringr", "forcats", "scales", "cowplot"
)

need_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(need_install) > 0) {
  install.packages(need_install, dependencies = TRUE)
}

library(tidyverse)
library(readxl)
library(patchwork)
library(stringr)
library(forcats)
library(scales)
library(cowplot)

############################
## 1. Font + PDF
############################
base_family <- "Arial"

if (.Platform$OS.type == "windows") {
  suppressWarnings({
    windowsFonts(Arial = windowsFont("Arial"))
  })
}

save_pdf_editable <- function(plot_obj, filename, width = 10, height = 7) {
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = base_family
  )
  print(plot_obj)
  dev.off()
}

############################
## 2. Paths
############################
# Input and output directories
# This script reads the refined sGCCA result table from "results/sgcca/refined_results".
# Figure files will be saved in "figures/sgcca".

infile <- file.path(
  "results", "sgcca", "refined_results",
  "sGCCA_postprocess_all_traits.xlsx"
)

outdir <- file.path("figures", "sgcca")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(infile))

############################
## 3. Read data
############################
sum_df  <- read_excel(infile, sheet = "summary_all_traits")
core_df <- read_excel(infile, sheet = "core_chains_all")
hub_df  <- read_excel(infile, sheet = "liver_hubs_all")

############################
## 4. Basic clean
############################
sum_df <- sum_df %>%
  mutate(
    Trait = as.character(Trait),
    model_level = as.character(model_level),
    model_pass = as.logical(model_pass),
    trait_liver_pass = as.logical(trait_liver_pass)
  )

core_df <- core_df %>%
  mutate(
    Trait = as.character(Trait),
    source_block = factor(as.character(source_block), levels = c("Rum", "Ile", "Col")),
    source_feature = as.character(source_feature),
    liver_feature = as.character(liver_feature),
    trait_feature = as.character(trait_feature),
    source_to_liver_direction = tolower(as.character(source_to_liver_direction)),
    liver_to_trait_direction  = tolower(as.character(liver_to_trait_direction)),
    source_to_trait_direction = tolower(as.character(source_to_trait_direction))
  )

hub_df <- hub_df %>%
  mutate(
    Trait = as.character(Trait),
    liver_feature = as.character(liver_feature),
    hub_flag = as.logical(hub_flag)
  )

if (!"abs_corr_source_liver" %in% names(core_df)) {
  if ("corr_source_liver" %in% names(core_df)) {
    core_df <- core_df %>% mutate(abs_corr_source_liver = abs(corr_source_liver))
  } else {
    core_df <- core_df %>% mutate(abs_corr_source_liver = NA_real_)
  }
}

if (!"abs_corr_liver_trait" %in% names(core_df)) {
  if ("corr_liver_trait" %in% names(core_df)) {
    core_df <- core_df %>% mutate(abs_corr_liver_trait = abs(corr_liver_trait))
  } else {
    core_df <- core_df %>% mutate(abs_corr_liver_trait = NA_real_)
  }
}

if (!"abs_corr_source_trait" %in% names(core_df)) {
  if ("corr_source_trait" %in% names(core_df)) {
    core_df <- core_df %>% mutate(abs_corr_source_trait = abs(corr_source_trait))
  } else {
    core_df <- core_df %>% mutate(abs_corr_source_trait = NA_real_)
  }
}

if (!"chain_score" %in% names(core_df)) {
  core_df <- core_df %>% mutate(chain_score = abs_corr_source_liver)
}

############################
## 5. Rename trait labels
############################
trait_label_map <- c(
  "GLU" = "GLU",
  "TBA" = "TBA",
  "TG" = "TG",
  "LDL" = "LDL",
  "TC" = "TC",
  "TailFat_g" = "TailFat",
  "TailFat_Carcass_g_per_kg" = "TailFat_Carcass",
  "TailFat_PreLive_g_per_kg" = "TailFat_PreLive"
)

sum_df <- sum_df %>%
  mutate(Trait_show = recode(Trait, !!!trait_label_map))

core_df <- core_df %>%
  mutate(Trait_show = recode(Trait, !!!trait_label_map))

hub_df <- hub_df %>%
  mutate(Trait_show = recode(Trait, !!!trait_label_map))

main_traits <- c("GLU", "TBA", "TailFat_g", "TG")
main_trait_labels <- c(
  "GLU" = "GLU",
  "TBA" = "TBA",
  "TailFat_g" = "TailFat",
  "TG" = "TG"
)

sum_main <- sum_df %>% filter(Trait %in% main_traits)
core_main <- core_df %>% filter(Trait %in% main_traits)
hub_main <- hub_df %>% filter(Trait %in% main_traits)

############################
## 6. Style
############################
col_source_fill <- c(
  "Rum" = "#BFD7A8",
  "Ile" = "#E1C97B",
  "Col" = "#C5D3E3"
)

col_source_text <- c(
  "Rum" = "#6E8558",
  "Ile" = "#9C8340",
  "Col" = "#6E859C"
)

## Close to the previous script
col_dir_line <- c(
  "positive" = "#C98A7D",
  "negative" = "#7FA6C9"
)

col_model <- c(
  "A_model_supported_primary_interpretation" = "#CBB7A7",
  "B_model_partly_supported_secondary_interpretation" = "#DDD5CD"
)

col_liver_node <- "#E8E1D9"
col_trait_node <- "#DDD0C5"

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = base_size + 1, face = "bold", colour = "#333333"),
      axis.text = element_text(size = base_size, colour = "#333333"),
      legend.title = element_text(size = base_size + 1, face = "bold", colour = "#333333"),
      legend.text = element_text(size = base_size, colour = "#333333"),
      strip.text = element_text(size = base_size + 1, face = "bold", colour = "#333333"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(10, 12, 10, 12)
    )
}

wrap_feature <- function(x, width = 22) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\.", " ") %>%
    str_wrap(width = width)
}

safe_max <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

############################
## 7. Figure 1
## clean overview
############################
fig1_df <- sum_df %>%
  mutate(
    Trait_show = factor(
      Trait_show,
      levels = rev(sum_df %>% arrange(trait_liver_corr) %>% pull(Trait_show))
    )
  )

p1 <- ggplot(fig1_df, aes(x = trait_liver_corr, y = Trait_show)) +
  geom_segment(
    aes(x = 0, xend = trait_liver_corr, y = Trait_show, yend = Trait_show),
    linewidth = 0.9, colour = "#CFC8C0"
  ) +
  geom_point(
    aes(fill = model_level),
    shape = 21, size = 5.2, colour = "#666666", stroke = 0.35
  ) +
  geom_text(
    aes(label = sprintf("%.3f", trait_liver_corr)),
    nudge_x = 0.040,   ## further to the right
    size = 3.6,
    family = base_family,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = col_model,
    labels = c(
      "A_model_supported_primary_interpretation" = "Level A",
      "B_model_partly_supported_secondary_interpretation" = "Level B"
    )
  ) +
  scale_x_continuous(
    limits = c(0, max(fig1_df$trait_liver_corr, na.rm = TRUE) + 0.16),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Trait-liver correlation",
    y = NULL,
    fill = "Model level"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold")
  )

save_pdf_editable(
  p1,
  filename = file.path(outdir, "Fig1_main_trait_overview_clean.pdf"),
  width = 8.6,
  height = 5.4
)

############################
## 8. Figure 2
## single-trait network, old-style logic
############################

## Use direct color values here to avoid name mismatches
line_col_positive <- "#C98A7D"
line_col_negative <- "#7FA6C9"

pick_top_sources_by_block <- function(df_trait,
                                      n_rum = 5,
                                      n_ile = 3,
                                      n_col = 5) {
  
  df_rank <- df_trait %>%
    mutate(
      rank_score = 0.70 * dplyr::coalesce(abs_corr_source_liver, 0) +
        0.30 * abs(dplyr::coalesce(chain_score, 0))
    ) %>%
    group_by(source_block, source_feature) %>%
    summarise(rank_score = safe_max(rank_score), .groups = "drop") %>%
    mutate(rank_score = dplyr::coalesce(rank_score, -999)) %>%
    arrange(source_block, desc(rank_score), source_feature)
  
  keep_rum <- df_rank %>%
    filter(source_block == "Rum") %>%
    slice_head(n = n_rum)
  
  keep_ile <- df_rank %>%
    filter(source_block == "Ile") %>%
    slice_head(n = n_ile)
  
  keep_col <- df_rank %>%
    filter(source_block == "Col") %>%
    slice_head(n = n_col)
  
  bind_rows(keep_rum, keep_ile, keep_col) %>%
    select(source_block, source_feature)
}

build_trait_network <- function(dat_trait, trait_name,
                                n_rum = 5, n_ile = 3, n_col = 5) {
  
  d0 <- dat_trait %>%
    filter(Trait == trait_name) %>%
    mutate(
      source_lab = wrap_feature(source_feature, width = 22),
      liver_lab  = wrap_feature(liver_feature, width = 16),
      trait_lab  = recode(trait_name, !!!main_trait_labels)
    )
  
  keep_df <- pick_top_sources_by_block(d0, n_rum = n_rum, n_ile = n_ile, n_col = n_col)
  
  d <- d0 %>%
    semi_join(keep_df, by = c("source_block", "source_feature"))
  
  make_source_nodes <- function(df, block_name, y_start) {
    tmp <- df %>%
      filter(source_block == block_name) %>%
      group_by(source_block, source_feature, source_lab) %>%
      summarise(score = safe_max(abs_corr_source_liver), .groups = "drop") %>%
      mutate(score = dplyr::coalesce(score, 0)) %>%
      arrange(desc(score), source_lab)
    
    if (nrow(tmp) == 0) {
      return(tibble(
        source_block = factor(character(), levels = c("Rum", "Ile", "Col")),
        source_feature = character(),
        source_lab = character(),
        score = numeric(),
        x = numeric(),
        y = numeric()
      ))
    }
    
    tmp %>%
      mutate(
        x = 1,
        y = seq(from = y_start, by = -1.35, length.out = n())
      )
  }
  
  nR <- d %>% filter(source_block == "Rum") %>% distinct(source_feature) %>% nrow()
  nI <- d %>% filter(source_block == "Ile") %>% distinct(source_feature) %>% nrow()
  
  rum_nodes <- make_source_nodes(d, "Rum", y_start = 12.0)
  ile_nodes <- make_source_nodes(d, "Ile", y_start = ifelse(nR > 0, 12.0 - nR * 1.35 - 1.6, 8.0))
  col_nodes <- make_source_nodes(
    d, "Col",
    y_start = ifelse(
      nI > 0,
      (ifelse(nR > 0, 12.0 - nR * 1.35 - 1.6, 8.0)) - nI * 1.35 - 1.6,
      ifelse(nR > 0, 12.0 - nR * 1.35 - 1.6, 8.0)
    )
  )
  
  source_nodes <- bind_rows(rum_nodes, ile_nodes, col_nodes)
  
  liver_nodes <- d %>%
    group_by(liver_feature, liver_lab) %>%
    summarise(score = safe_max(abs_corr_liver_trait), .groups = "drop") %>%
    mutate(score = dplyr::coalesce(score, 0)) %>%
    arrange(desc(score), liver_lab)
  
  if (nrow(liver_nodes) == 0) {
    liver_nodes <- tibble(
      liver_feature = paste0("Hub_", trait_name),
      liver_lab = paste0("Hub_", trait_name),
      score = 0
    )
  }
  
  liver_nodes <- liver_nodes %>%
    mutate(
      x = 2,
      y = seq(
        from = mean(range(source_nodes$y, na.rm = TRUE)) + 1.1,
        by = -1.5,
        length.out = n()
      )
    )
  
  trait_node <- tibble(
    trait_lab = recode(trait_name, !!!main_trait_labels),
    x = 3,
    y = mean(liver_nodes$y)
  )
  
  edge_sl <- d %>%
    left_join(
      source_nodes %>% select(source_block, source_feature, x_s = x, y_s = y),
      by = c("source_block", "source_feature")
    ) %>%
    left_join(
      liver_nodes %>% select(liver_feature, x_l = x, y_l = y),
      by = "liver_feature"
    ) %>%
    filter(!is.na(x_s), !is.na(y_s), !is.na(x_l), !is.na(y_l)) %>%
    mutate(
      line_dir = factor(
        ifelse(source_to_liver_direction == "positive", "positive", "negative"),
        levels = c("positive", "negative")
      ),
      line_weight = dplyr::coalesce(abs_corr_source_liver, 0.12)
    )
  
  edge_lt <- d %>%
    distinct(liver_feature, liver_lab, liver_to_trait_direction, abs_corr_liver_trait) %>%
    left_join(
      liver_nodes %>% select(liver_feature, x_l = x, y_l = y),
      by = "liver_feature"
    ) %>%
    filter(!is.na(x_l), !is.na(y_l)) %>%
    mutate(
      line_dir = factor(
        ifelse(liver_to_trait_direction == "positive", "positive", "negative"),
        levels = c("positive", "negative")
      ),
      line_weight = dplyr::coalesce(abs_corr_liver_trait, 0.12)
    )
  
  list(
    source_nodes = source_nodes,
    liver_nodes = liver_nodes,
    trait_node = trait_node,
    edge_sl = edge_sl,
    edge_lt = edge_lt,
    y_top = max(source_nodes$y, na.rm = TRUE) + 1.8,
    y_bottom = min(source_nodes$y, na.rm = TRUE) - 0.8
  )
}

plot_trait_network <- function(dat_trait, trait_name,
                               n_rum = 5, n_ile = 3, n_col = 5,
                               show_column_titles = TRUE,
                               show_legend = TRUE) {
  
  net <- build_trait_network(
    dat_trait = dat_trait,
    trait_name = trait_name,
    n_rum = n_rum,
    n_ile = n_ile,
    n_col = n_col
  )
  
  source_nodes <- net$source_nodes
  liver_nodes  <- net$liver_nodes
  trait_node   <- net$trait_node
  edge_sl      <- net$edge_sl
  edge_lt      <- net$edge_lt
  y_top        <- net$y_top
  y_bottom     <- net$y_bottom
  
  p <- ggplot() +
    geom_curve(
      data = edge_sl,
      aes(
        x = x_s, y = y_s, xend = x_l, yend = y_l,
        colour = line_dir,
        linetype = line_dir,
        linewidth = line_weight
      ),
      curvature = 0.18,
      alpha = 0.98,
      lineend = "round"
    ) +
    geom_curve(
      data = edge_lt,
      aes(
        x = x_l, y = y_l, xend = trait_node$x[1], yend = trait_node$y[1],
        colour = line_dir,
        linetype = line_dir,
        linewidth = line_weight
      ),
      curvature = 0.12,
      alpha = 0.98,
      lineend = "round"
    ) +
    geom_point(
      data = source_nodes,
      aes(x = x, y = y, fill = source_block),
      shape = 21, size = 4.9, stroke = 0.8, colour = "#5F5F5F"
    ) +
    geom_point(
      data = liver_nodes,
      aes(x = x, y = y),
      shape = 21, size = 6.2, stroke = 0.9,
      fill = col_liver_node, colour = "#5F5F5F"
    ) +
    geom_point(
      data = trait_node,
      aes(x = x, y = y),
      shape = 21, size = 7.0, stroke = 0.95,
      fill = col_trait_node, colour = "#5F5F5F"
    ) +
    geom_text(
      data = source_nodes,
      aes(x = x - 0.12, y = y, label = source_lab),
      colour = unname(col_source_text[as.character(source_nodes$source_block)]),
      family = base_family,
      size = 3.75,
      fontface = "italic",
      hjust = 1,
      vjust = 0.5,
      lineheight = 0.90,
      show.legend = FALSE
    ) +
    geom_text(
      data = liver_nodes,
      aes(x = x, y = y + 0.60, label = liver_lab),
      family = base_family,
      size = 3.95,
      fontface = "bold.italic",
      hjust = 0.5,
      vjust = 0
    ) +
    geom_text(
      data = trait_node,
      aes(x = x + 0.11, y = y, label = trait_lab),
      family = base_family,
      size = 4.45,
      fontface = "bold",
      hjust = 0,
      vjust = 0.5
    ) +
    scale_fill_manual(values = col_source_fill, name = "Source block") +
    scale_colour_manual(
      values = c(
        positive = line_col_positive,
        negative = line_col_negative
      ),
      breaks = c("positive", "negative"),
      labels = c("Positive", "Negative"),
      name = "Direction",
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = c(
        positive = "solid",
        negative = "22"
      ),
      breaks = c("positive", "negative"),
      labels = c("Positive", "Negative"),
      name = "Direction",
      drop = FALSE
    ) +
    scale_linewidth_continuous(range = c(0.45, 1.25), guide = "none") +
    coord_cartesian(
      xlim = c(0.12, 3.45),
      ylim = c(y_bottom, y_top + 0.45),
      clip = "off"
    ) +
    theme_void(base_size = 14) +
    theme(
      text = element_text(family = base_family),
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.box = "horizontal",
      legend.title = element_text(size = 11.8, face = "bold"),
      legend.text = element_text(size = 10.8),
      plot.margin = margin(16, 70, 16, 150)
    )
  
  if (show_column_titles) {
    p <- p +
      annotate(
        "text", x = 1, y = y_top, label = "Microbial source",
        family = base_family, fontface = "bold", size = 4.6
      ) +
      annotate(
        "text", x = 2, y = y_top, label = "Liver hub",
        family = base_family, fontface = "bold", size = 4.6
      ) +
      annotate(
        "text", x = 3, y = y_top, label = "Trait",
        family = base_family, fontface = "bold", size = 4.6
      )
  }
  
  p
}

## Single figure
p2_glu  <- plot_trait_network(core_main, "GLU",        n_rum = 2, n_ile = 0, n_col = 5, show_column_titles = TRUE, show_legend = TRUE)
p2_tba  <- plot_trait_network(core_main, "TBA",        n_rum = 5, n_ile = 0, n_col = 5, show_column_titles = TRUE, show_legend = TRUE)
p2_tail <- plot_trait_network(core_main, "TailFat_g",  n_rum = 4, n_ile = 1, n_col = 5, show_column_titles = TRUE, show_legend = TRUE)
p2_tg   <- plot_trait_network(core_main, "TG",         n_rum = 5, n_ile = 0, n_col = 5, show_column_titles = TRUE, show_legend = TRUE)

save_pdf_editable(p2_glu,  file.path(outdir, "Fig2_network_GLU.pdf"),      width = 8.2, height = 6.8)
save_pdf_editable(p2_tba,  file.path(outdir, "Fig2_network_TBA.pdf"),      width = 8.2, height = 6.8)
save_pdf_editable(p2_tail, file.path(outdir, "Fig2_network_TailFat.pdf"),  width = 8.2, height = 6.8)
save_pdf_editable(p2_tg,   file.path(outdir, "Fig2_network_TG.pdf"),       width = 8.2, height = 6.8)

## Composite layout: subplots without legends
p2_glu_c  <- plot_trait_network(core_main, "GLU",       n_rum = 2, n_ile = 0, n_col = 5, show_column_titles = TRUE,  show_legend = FALSE)
p2_tba_c  <- plot_trait_network(core_main, "TBA",       n_rum = 5, n_ile = 0, n_col = 5, show_column_titles = TRUE,  show_legend = FALSE)
p2_tail_c <- plot_trait_network(core_main, "TailFat_g", n_rum = 4, n_ile = 1, n_col = 5, show_column_titles = FALSE, show_legend = FALSE)
p2_tg_c   <- plot_trait_network(core_main, "TG",        n_rum = 5, n_ile = 0, n_col = 5, show_column_titles = FALSE, show_legend = FALSE)

## Create one complete standalone legend
legend_fill_plot <- ggplot(
  tibble(
    x = 1:3,
    y = 1,
    grp = factor(c("Rum", "Ile", "Col"), levels = c("Rum", "Ile", "Col"))
  ),
  aes(x = x, y = y, fill = grp)
) +
  geom_point(shape = 21, size = 4.8, colour = "#5F5F5F") +
  scale_fill_manual(values = col_source_fill, name = "Source block") +
  theme_void(base_family = base_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11.8, face = "bold"),
    legend.text = element_text(size = 10.8)
  )

legend_line_plot <- ggplot(
  tibble(
    x = c(1, 1),
    xend = c(2, 2),
    y = c(1, 2),
    yend = c(1, 2),
    grp = factor(c("positive", "negative"), levels = c("positive", "negative"))
  ),
  aes(x = x, xend = xend, y = y, yend = yend, colour = grp, linetype = grp)
) +
  geom_segment(linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      positive = line_col_positive,
      negative = line_col_negative
    ),
    breaks = c("positive", "negative"),
    labels = c("Positive", "Negative"),
    name = "Direction",
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = c(
      positive = "solid",
      negative = "22"
    ),
    breaks = c("positive", "negative"),
    labels = c("Positive", "Negative"),
    name = "Direction",
    drop = FALSE
  ) +
  theme_void(base_family = base_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11.8, face = "bold"),
    legend.text = element_text(size = 10.8)
  )

legend_fill  <- cowplot::get_legend(legend_fill_plot)
legend_line  <- cowplot::get_legend(legend_line_plot)
legend_comb  <- patchwork::wrap_elements(full = cowplot::plot_grid(legend_fill, legend_line, nrow = 1, rel_widths = c(1, 1)))

p2_combined <- (p2_glu_c + p2_tba_c) / (p2_tail_c + p2_tg_c) / legend_comb +
  plot_layout(heights = c(1, 1, 0.14))

save_pdf_editable(
  p2_combined,
  file.path(outdir, "Fig2_network_2x2_combined.pdf"),
  width = 13.2,
  height = 10.2
)

############################
## 9. Figure 3
## source mode summary
############################
fig3_df <- core_main %>%
  count(Trait_show, source_block, name = "n_chain") %>%
  complete(
    Trait_show = c("GLU", "TBA", "TailFat", "TG"),
    source_block = levels(core_main$source_block),
    fill = list(n_chain = 0)
  ) %>%
  mutate(
    Trait_show = factor(Trait_show, levels = c("GLU", "TBA", "TailFat", "TG"))
  )

p3 <- ggplot(fig3_df, aes(x = Trait_show, y = n_chain, fill = source_block)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = ifelse(n_chain > 0, n_chain, "")),
    position = position_stack(vjust = 0.5),
    family = base_family,
    size = 3.6,
    colour = "#333333"
  ) +
  scale_fill_manual(values = col_source_fill) +
  labs(
    x = NULL,
    y = "Core chain count",
    fill = "Source block"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(face = "bold")
  )

save_pdf_editable(
  p3,
  file.path(outdir, "Fig3_source_mode_summary.pdf"),
  width = 7.3,
  height = 5.5
)

############################
## 10. Figure 4
## bubble matrix
############################
fig4_df <- hub_main %>%
  filter(hub_flag) %>%
  mutate(
    Trait_show = factor(Trait_show, levels = c("GLU", "TBA", "TailFat", "TG")),
    liver_feature = fct_reorder(liver_feature, mean_chain_score, .fun = max, .desc = TRUE)
  )

p4 <- ggplot(fig4_df, aes(x = Trait_show, y = liver_feature)) +
  geom_point(
    aes(size = n_source_links, fill = mean_chain_score),
    shape = 21, colour = "#6A6A6A", stroke = 0.35, alpha = 0.96
  ) +
  geom_text(
    aes(label = n_source_links),
    family = base_family,
    size = 3.3,
    fontface = "bold",
    colour = "#333333"
  ) +
  scale_fill_gradient(
    low = "#F5F1EB",
    high = "#BA917E",
    name = "Mean chain score"
  ) +
  scale_size_continuous(
    range = c(5.0, 12.5),
    name = "Source link n"
  ) +
  labs(
    x = NULL,
    y = "Liver hub"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold.italic")
  )

save_pdf_editable(
  p4,
  file.path(outdir, "Fig4_liver_hub_bubble_matrix.pdf"),
  width = 7.6,
  height = 5.3
)

############################
## 11. Supplementary Figure S1
############################
figS1_df <- core_main %>%
  count(source_block, source_feature, Trait_show, name = "n_chain") %>%
  group_by(source_block, source_feature) %>%
  mutate(
    total_chain = sum(n_chain, na.rm = TRUE),
    n_trait = sum(n_chain > 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(n_trait >= 2) %>%
  distinct(source_block, source_feature, Trait_show, n_chain, total_chain, n_trait)

top_taxa <- figS1_df %>%
  distinct(source_block, source_feature, total_chain, n_trait) %>%
  arrange(desc(n_trait), desc(total_chain), source_block, source_feature) %>%
  slice_head(n = 15) %>%
  pull(source_feature)

figS1_df2 <- figS1_df %>%
  filter(source_feature %in% top_taxa) %>%
  mutate(
    Trait_show = factor(Trait_show, levels = c("GLU", "TBA", "TailFat", "TG")),
    source_feature = factor(
      source_feature,
      levels = figS1_df %>%
        filter(source_feature %in% top_taxa) %>%
        distinct(source_feature, total_chain, n_trait) %>%
        arrange(n_trait, total_chain) %>%
        pull(source_feature)
    )
  )

pS1 <- ggplot(figS1_df2, aes(x = Trait_show, y = source_feature)) +
  geom_point(
    aes(size = n_chain, fill = source_block),
    shape = 21, colour = "#696969", stroke = 0.35, alpha = 0.96
  ) +
  scale_fill_manual(values = col_source_fill) +
  scale_size_continuous(range = c(3.5, 10.5)) +
  labs(
    x = NULL,
    y = "Recurrent source microbe",
    fill = "Source block",
    size = "Chain count"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic")
  )

save_pdf_editable(
  pS1,
  file.path(outdir, "FigS1_recurrent_source_microbes.pdf"),
  width = 9.0,
  height = 6.8
)

############################
## 12. Export check tables
############################
write.csv(sum_df,    file.path(outdir, "plot_check_summary_all_traits.csv"),   row.names = FALSE, fileEncoding = "UTF-8")
write.csv(core_df,   file.path(outdir, "plot_check_core_chains_all.csv"),      row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hub_df,    file.path(outdir, "plot_check_liver_hubs_all.csv"),       row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fig3_df,   file.path(outdir, "plot_check_fig3_source_mode.csv"),     row.names = FALSE, fileEncoding = "UTF-8")
write.csv(figS1_df2, file.path(outdir, "plot_check_figS1_recurrent.csv"),      row.names = FALSE, fileEncoding = "UTF-8")

message("All figures finished.")
message("Output directory: ", outdir)