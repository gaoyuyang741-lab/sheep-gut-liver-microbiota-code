############################################################
############################################################

rm(list = ls())
gc()

############################
## 0. Packages
############################
pkg_needed <- c(
  "ggplot2", "dplyr", "tidyr", "readxl", "stringr", "forcats",
  "patchwork", "ggrepel", "scales", "tibble", "VennDiagram", "grid"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, dependencies = TRUE)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(forcats)
library(patchwork)
library(ggrepel)
library(scales)
library(tibble)
library(VennDiagram)
library(grid)

############################
## 1. Paths
############################
# Input and output directories
# This script reads pairwise analysis result files from "results/pairwise".
# Figure files will be saved in "figures/pairwise".

indir  <- file.path("results", "pairwise")
outdir <- file.path("figures", "pairwise")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

file_module <- file.path(indir, "Module_summary.xlsx")
file_sig    <- file.path(indir, "Significant_associations.xlsx")
file_kegg   <- file.path(indir, "Liver_KEGG_gene_list.xlsx")

############################
## 2. Global settings
############################
base_family <- "Arial"
base_size   <- 11

theme_pairwise <- function(base_size = 11) {
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      strip.text = element_text(color = "black"),
      legend.key = element_blank()
    )
}

theme_set(theme_pairwise(base_size))

save_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    dpi = 300
  )
}

############################
## 3. Color palettes
############################
col_segment <- c(
  "Rum" = "#8FB7AA",
  "Ile" = "#D7B08A",
  "Col" = "#A8BFA3"
)

col_host <- c(
  "Liver" = "#8FA7C6",
  "Blood" = "#D8A39D",
  "Tail fat" = "#C5AFD7"
)

col_direction <- c(
  "Positive" = "#D5A09A",
  "Negative" = "#8EA8C9"
)

col_set <- c(
  "Rum-specific" = "#8FB7AA",
  "Shared core" = "#C8B07A",
  "Col-specific" = "#A8BFA3"
)

col_category <- c(
  "BA synthesis/regulation" = "#D7B06A",
  "BA transport/modification" = "#E4C995",
  "Lipogenesis" = "#D8A49C",
  "FA oxidation/catabolism" = "#9FBC9A",
  "Gluconeogenesis/energy" = "#9CB2D1",
  "Cholesterol/sterol metabolism" = "#B8A2CC",
  "TG/lipid droplet metabolism" = "#CFB082",
  "Other metabolism" = "#C8C8C8"
)

############################
## 4. Read data
############################
module_df <- readxl::read_excel(file_module, sheet = 1) %>%
  as.data.frame()

sig_df <- readxl::read_excel(file_sig, sheet = 1) %>%
  as.data.frame()

if (file.exists(file_kegg)) {
  kegg_df <- readxl::read_excel(file_kegg, sheet = 1) %>%
    as.data.frame()
} else {
  kegg_df <- NULL
}

############################
## 5. Basic parsing
############################
module_order <- c(
  "Rum_vs_Liver", "Ile_vs_Liver", "Colon_vs_Liver",
  "Rum_vs_Blood", "Ile_vs_Blood", "Colon_vs_Blood",
  "Rum_vs_Tailfat", "Ile_vs_Tailfat", "Colon_vs_Tailfat",
  "Liver_vs_Blood", "Liver_vs_Tailfat"
)

module_label_map <- c(
  "Rum_vs_Liver"     = "Rum–Liver",
  "Ile_vs_Liver"     = "Ile–Liver",
  "Colon_vs_Liver"   = "Col–Liver",
  "Rum_vs_Blood"     = "Rum–Blood",
  "Ile_vs_Blood"     = "Ile–Blood",
  "Colon_vs_Blood"   = "Col–Blood",
  "Rum_vs_Tailfat"   = "Rum–Tail fat",
  "Ile_vs_Tailfat"   = "Ile–Tail fat",
  "Colon_vs_Tailfat" = "Col–Tail fat",
  "Liver_vs_Blood"   = "Liver–Blood",
  "Liver_vs_Tailfat" = "Liver–Tail fat"
)

parse_module <- function(x) {
  sp <- stringr::str_split(x, "_vs_", simplify = TRUE)
  data.frame(
    module = x,
    left_block = sp[, 1],
    right_block = sp[, 2],
    stringsAsFactors = FALSE
  )
}

module_meta <- parse_module(module_df$module)

module_df2 <- module_df %>%
  left_join(module_meta, by = "module") %>%
  mutate(
    module = factor(module, levels = module_order),
    module_label = module_label_map[as.character(module)],
    segment = case_when(
      left_block == "Rum" ~ "Rum",
      left_block == "Ile" ~ "Ile",
      left_block == "Colon" ~ "Col",
      TRUE ~ left_block
    ),
    host = case_when(
      right_block == "Liver" ~ "Liver",
      right_block == "Blood" ~ "Blood",
      right_block == "Tailfat" ~ "Tail fat",
      TRUE ~ right_block
    ),
    significant_rate = significant_pairs / tested_pairs
  )

sig_df2 <- sig_df %>%
  left_join(module_meta, by = "module") %>%
  mutate(
    module = factor(module, levels = module_order),
    module_label = module_label_map[as.character(module)],
    segment = case_when(
      left_block == "Rum" ~ "Rum",
      left_block == "Ile" ~ "Ile",
      left_block == "Colon" ~ "Col",
      TRUE ~ left_block
    ),
    host = case_when(
      right_block == "Liver" ~ "Liver",
      right_block == "Blood" ~ "Blood",
      right_block == "Tailfat" ~ "Tail fat",
      TRUE ~ right_block
    ),
    direction2 = case_when(
      tolower(direction) %in% c("positive", "pos") ~ "Positive",
      tolower(direction) %in% c("negative", "neg") ~ "Negative",
      TRUE ~ NA_character_
    )
  )

module_main <- module_df2 %>%
  filter(
    left_block %in% c("Rum", "Ile", "Colon"),
    right_block %in% c("Liver", "Blood", "Tailfat")
  )

sig_main <- sig_df2 %>%
  filter(
    left_block %in% c("Rum", "Ile", "Colon"),
    right_block %in% c("Liver", "Blood", "Tailfat")
  )

############################
## 6. FIGURE 1
## Diverging bar chart for module-level counts
############################
fig1_df <- module_main %>%
  select(module, module_label, significant_pairs, positive_pairs, negative_pairs) %>%
  pivot_longer(
    cols = c(positive_pairs, negative_pairs),
    names_to = "direction_type",
    values_to = "n"
  ) %>%
  mutate(
    direction = ifelse(direction_type == "positive_pairs", "Positive", "Negative"),
    value = ifelse(direction == "Negative", -n, n)
  )

fig1_levels <- module_main %>%
  arrange(module) %>%
  pull(module_label)

fig1_df <- fig1_df %>%
  mutate(module_label = factor(module_label, levels = fig1_levels))

fig1_total <- module_main %>%
  mutate(module_label = factor(module_label, levels = fig1_levels)) %>%
  arrange(module_label)

p_fig1 <- ggplot(fig1_df, aes(x = module_label, y = value, fill = direction)) +
  geom_col(width = 0.64, color = "black", linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "black") +
  geom_text(
    data = fig1_total,
    aes(
      x = module_label,
      y = pmax(positive_pairs, 3) + 8,
      label = significant_pairs
    ),
    inherit.aes = FALSE,
    family = base_family,
    size = 3.4
  ) +
  scale_fill_manual(values = col_direction) +
  labs(
    x = NULL,
    y = "Significant associations",
    fill = NULL
  ) +
  coord_flip() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

save_pdf(
  p_fig1,
  file.path(outdir, "Fig1_pairwise_module_counts.pdf"),
  width = 8.0,
  height = 5.8
)

############################
## 7. FIGURE 2
## Connection spectrum by gut segment
############################
fig2_df <- module_main %>%
  group_by(segment, host) %>%
  summarise(significant_pairs = sum(significant_pairs), .groups = "drop") %>%
  mutate(
    segment = factor(segment, levels = c("Rum", "Ile", "Col")),
    host = factor(host, levels = c("Liver", "Blood", "Tail fat"))
  )

fig2_total <- fig2_df %>%
  group_by(segment) %>%
  summarise(total = sum(significant_pairs), .groups = "drop")

p_fig2 <- ggplot(fig2_df, aes(x = segment, y = significant_pairs, fill = host)) +
  geom_col(width = 0.60, color = "black", linewidth = 0.3) +
  geom_text(
    data = fig2_total,
    aes(x = segment, y = total + 6, label = total),
    inherit.aes = FALSE,
    family = base_family,
    size = 3.5
  ) +
  scale_fill_manual(values = col_host) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = scales::comma
  ) +
  labs(
    x = NULL,
    y = "Significant associations",
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

save_pdf(
  p_fig2,
  file.path(outdir, "Fig2_pairwise_connection_spectrum.pdf"),
  width = 5.8,
  height = 5.4
)

############################
## 8. FIGURE 3
## Bubble plot of key liver genes across gut segments
############################
liver_sig <- sig_main %>%
  filter(host == "Liver")

gene_seg_stat <- liver_sig %>%
  group_by(segment, feature_y) %>%
  summarise(
    n_sig = n(),
    max_abs_rho = max(abs_rho, na.rm = TRUE),
    mean_rho = mean(rho, na.rm = TRUE),
    .groups = "drop"
  )

gene_rank <- gene_seg_stat %>%
  group_by(feature_y) %>%
  summarise(
    total_n = sum(n_sig),
    present_seg = n(),
    max_abs_rho = max(max_abs_rho),
    .groups = "drop"
  ) %>%
  arrange(desc(present_seg), desc(total_n), desc(max_abs_rho))

rum_liver_genes <- unique(liver_sig$feature_y[liver_sig$segment == "Rum"])
col_liver_genes <- unique(liver_sig$feature_y[liver_sig$segment == "Col"])

common_genes <- intersect(rum_liver_genes, col_liver_genes)
rum_only_genes <- setdiff(rum_liver_genes, col_liver_genes)
col_only_genes <- setdiff(col_liver_genes, rum_liver_genes)

rum_only_top <- gene_seg_stat %>%
  filter(segment == "Rum", feature_y %in% rum_only_genes) %>%
  arrange(desc(n_sig), desc(max_abs_rho)) %>%
  slice_head(n = 5) %>%
  pull(feature_y)

col_only_top <- gene_seg_stat %>%
  filter(segment == "Col", feature_y %in% col_only_genes) %>%
  arrange(desc(n_sig), desc(max_abs_rho)) %>%
  slice_head(n = 5) %>%
  pull(feature_y)

gene_show <- unique(c(common_genes, rum_only_top, col_only_top))

gene_show <- gene_rank %>%
  filter(feature_y %in% gene_show) %>%
  slice_head(n = 20) %>%
  pull(feature_y)

fig3_df <- gene_seg_stat %>%
  filter(feature_y %in% gene_show) %>%
  mutate(
    segment = factor(segment, levels = c("Rum", "Ile", "Col"))
  )

gene_order_fig3 <- gene_rank %>%
  filter(feature_y %in% gene_show) %>%
  arrange(desc(present_seg), desc(total_n), desc(max_abs_rho)) %>%
  pull(feature_y)

fig3_df$feature_y <- factor(fig3_df$feature_y, levels = rev(gene_order_fig3))

p_fig3 <- ggplot(fig3_df, aes(x = segment, y = feature_y)) +
  geom_point(
    aes(size = n_sig, fill = mean_rho),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  scale_size_continuous(range = c(2.4, 8.0)) +
  scale_fill_gradient2(
    low = "#8FA7C6",
    mid = "white",
    high = "#D7A19A",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "No. of associated genera",
    fill = "Mean rho"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9, face = "italic"),
    plot.margin = margin(8, 8, 8, 8)
  )

save_pdf(
  p_fig3,
  file.path(outdir, "Fig3_pairwise_liver_gene_bubble.pdf"),
  width = 7.0,
  height = 6.8
)

############################
## 9. FIGURE 4
## Standard circular Venn + functional map
############################
gene_category_map <- c(
  "SCARB1"   = "BA transport/modification",
  "SLC10A1"  = "BA transport/modification",
  "SULT2A1"  = "BA transport/modification",
  "UGT2B10"  = "BA transport/modification",
  "UGT2C1"   = "BA transport/modification",
  "CYP8B1"   = "BA synthesis/regulation",
  "NR0B2"    = "BA synthesis/regulation",
  "HADHA"    = "FA oxidation/catabolism",
  "HADHB"    = "FA oxidation/catabolism",
  "PCK1"     = "Gluconeogenesis/energy",
  "PCK2"     = "Gluconeogenesis/energy",
  "LPL"      = "TG/lipid droplet metabolism",
  "HSD17B12" = "Lipogenesis",
  "PNPLA2"   = "TG/lipid droplet metabolism",
  "CYP7A1"   = "BA synthesis/regulation",
  "CYP7B1"   = "BA synthesis/regulation",
  "NR1H3"    = "BA synthesis/regulation",
  "NR1H4"    = "BA synthesis/regulation",
  "FASN"     = "Lipogenesis",
  "ACACA"    = "Lipogenesis",
  "HMGCS1"   = "Cholesterol/sterol metabolism",
  "INSIG2"   = "Cholesterol/sterol metabolism",
  "PKLR"     = "Gluconeogenesis/energy",
  "SQLE"     = "Cholesterol/sterol metabolism",
  "CPT1A"    = "FA oxidation/catabolism",
  "CPT2"     = "FA oxidation/catabolism",
  "ACOX2"    = "FA oxidation/catabolism",
  "AKR1D1"   = "BA transport/modification",
  "AMACR"    = "FA oxidation/catabolism",
  "CYP27A1"  = "BA synthesis/regulation",
  "PC"       = "Gluconeogenesis/energy",
  "SORT1"    = "TG/lipid droplet metabolism",
  "UGT1A1"   = "BA transport/modification",
  "UGT2A3"   = "BA transport/modification",
  "UGT2B31"  = "BA transport/modification",
  "ACAT1"    = "Cholesterol/sterol metabolism",
  "AGPAT2"   = "Lipogenesis"
)

func_order <- c(
  "BA synthesis/regulation",
  "BA transport/modification",
  "Cholesterol/sterol metabolism",
  "FA oxidation/catabolism",
  "Gluconeogenesis/energy",
  "Lipogenesis",
  "TG/lipid droplet metabolism",
  "Other metabolism"
)

rum_df <- data.frame(
  set_group = "Rum-specific",
  gene = rum_only_genes,
  stringsAsFactors = FALSE
)

shared_df <- data.frame(
  set_group = "Shared core",
  gene = common_genes,
  stringsAsFactors = FALSE
)

col_df <- data.frame(
  set_group = "Col-specific",
  gene = col_only_genes,
  stringsAsFactors = FALSE
)

fig4_df <- bind_rows(rum_df, shared_df, col_df) %>%
  mutate(
    category = gene_category_map[gene],
    category = ifelse(is.na(category) & !is.na(gene), "Other metabolism", category),
    category = factor(category, levels = func_order),
    set_group = factor(set_group, levels = c("Rum-specific", "Shared core", "Col-specific"))
  ) %>%
  arrange(set_group, category, gene)

pad_and_rank <- function(df, grp) {
  df_sub <- df %>%
    filter(set_group == grp) %>%
    arrange(category, gene)
  n_sub <- nrow(df_sub)
  df_sub$row_id <- seq_len(n_sub)
  df_sub
}

fig4_r <- pad_and_rank(fig4_df, "Rum-specific")
fig4_s <- pad_and_rank(fig4_df, "Shared core")
fig4_c <- pad_and_rank(fig4_df, "Col-specific")

max_n <- max(nrow(fig4_r), nrow(fig4_s), nrow(fig4_c))

pad_to_max <- function(df_sub, grp, max_n) {
  if (nrow(df_sub) < max_n) {
    df_sub <- bind_rows(
      df_sub,
      data.frame(
        set_group = grp,
        gene = NA_character_,
        category = factor(NA_character_, levels = func_order),
        row_id = seq(nrow(df_sub) + 1, max_n),
        stringsAsFactors = FALSE
      )
    )
  }
  df_sub
}

fig4_r <- pad_to_max(fig4_r, "Rum-specific", max_n)
fig4_s <- pad_to_max(fig4_s, "Shared core", max_n)
fig4_c <- pad_to_max(fig4_c, "Col-specific", max_n)

fig4_plot_df <- bind_rows(fig4_r, fig4_s, fig4_c) %>%
  mutate(
    set_group = factor(set_group, levels = c("Rum-specific", "Shared core", "Col-specific"))
  )

## ---------- top venn ----------
## Use a square viewport so circles stay true circles.
make_top_venn_grob <- function() {
  grid::grid.grabExpr({
    grid::grid.newpage()
    
    ## Outer canvas
    grid::pushViewport(grid::viewport(
      x = 0.5, y = 0.50,
      width = unit(1, "npc"),
      height = unit(1, "npc"),
      just = c("center", "center")
    ))
    
    ## Inner square region: this is the key to avoid ellipses
    grid::pushViewport(grid::viewport(
      x = 0.38, y = 0.50,
      width = unit(1.3, "snpc"),
      height = unit(1.3, "snpc"),
      just = c("center", "center")
    ))
    
    venn_g <- VennDiagram::draw.pairwise.venn(
      area1 = length(col_liver_genes),   # left circle = Col
      area2 = length(rum_liver_genes),   # right circle = Rum
      cross.area = length(common_genes),
      category = c("Col", "Rum"),
      fill = c("#A8BFA3", "#8FB7AA"),
      alpha = c(0.45, 0.45),
      col = c("black", "black"),
      lwd = 1.35,
      lty = "solid",
      scaled = FALSE,
      rotation.degree = 0,
      cex = 2.0,
      fontface = "bold",
      fontfamily = base_family,
      cat.cex = 1.9,
      cat.fontface = "bold",
      cat.fontfamily = base_family,
      cat.col = c("black", "black"),
      cat.pos = c(105, 75),
      cat.dist = c(0.06, 0.06),
      margin = 0.02,
      ind = FALSE
    )
    
    grid::grid.draw(venn_g)
    
    grid::popViewport()
    grid::popViewport()
  })
}

venn_top_grob <- make_top_venn_grob()

p4_top <- patchwork::wrap_elements(full = venn_top_grob) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

## ---------- bottom functional map ----------
p4_bottom <- ggplot() +
  geom_tile(
    data = subset(fig4_plot_df, !is.na(gene)),
    aes(x = set_group, y = row_id),
    width = 0.78,
    height = 0.78,
    fill = "white",
    color = "black",
    linewidth = 0.25
  ) +
  geom_tile(
    data = subset(fig4_plot_df, !is.na(gene)),
    aes(x = as.numeric(set_group) - 0.31, y = row_id, fill = category),
    width = 0.08,
    height = 0.78,
    color = NA
  ) +
  geom_point(
    data = subset(fig4_plot_df, !is.na(gene)),
    aes(x = as.numeric(set_group) - 0.23, y = row_id, fill = category),
    shape = 21,
    size = 1.9,
    color = "black",
    stroke = 0.2
  ) +
  geom_text(
    data = subset(fig4_plot_df, !is.na(gene)),
    aes(x = set_group, y = row_id, label = gene),
    family = base_family,
    fontface = "italic",
    size = 3.0
  ) +
  scale_fill_manual(values = col_category, drop = FALSE) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, fill = "Functional class") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10.8),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 8, 8, 8)
  )

p_fig4 <- p4_top / p4_bottom + plot_layout(heights = c(1.95, 4.05))

save_pdf(
  p_fig4,
  file.path(outdir, "Fig4_pairwise_overlap_functional_map.pdf"),
  width = 9.4,
  height = 8.3
)

############################
## 10. FIGURE 5
## Cleaner bipartite network plots
############################
build_bipartite_network_plot <- function(
    data_sig,
    module_name,
    header_label,
    top_gene_n = 8,
    top_microbe_n = 9,
    max_edge_n = 22,
    node_size = 3.8
) {
  df0 <- data_sig %>%
    filter(module == module_name)
  
  if (nrow(df0) == 0) {
    return(ggplot() + theme_void())
  }
  
  gene_rank0 <- df0 %>%
    group_by(feature_y) %>%
    summarise(
      edge_n = n(),
      max_abs_rho = max(abs_rho, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(edge_n), desc(max_abs_rho)) %>%
    slice_head(n = top_gene_n)
  
  microbe_rank0 <- df0 %>%
    group_by(feature_x) %>%
    summarise(
      edge_n = n(),
      max_abs_rho = max(abs_rho, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(edge_n), desc(max_abs_rho)) %>%
    slice_head(n = top_microbe_n)
  
  df1 <- df0 %>%
    filter(
      feature_y %in% gene_rank0$feature_y,
      feature_x %in% microbe_rank0$feature_x
    ) %>%
    arrange(desc(abs_rho)) %>%
    slice_head(n = max_edge_n)
  
  if (nrow(df1) < 8) {
    df1 <- df0 %>%
      filter(feature_y %in% gene_rank0$feature_y) %>%
      arrange(desc(abs_rho)) %>%
      slice_head(n = max_edge_n)
    
    microbe_rank0 <- df1 %>%
      group_by(feature_x) %>%
      summarise(
        edge_n = n(),
        max_abs_rho = max(abs_rho, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(edge_n), desc(max_abs_rho)) %>%
      slice_head(n = top_microbe_n)
    
    df1 <- df1 %>%
      filter(feature_x %in% microbe_rank0$feature_x)
  }
  
  microbes <- df1 %>%
    group_by(feature_x) %>%
    summarise(weight = sum(abs_rho, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(weight), feature_x)
  
  genes <- df1 %>%
    group_by(feature_y) %>%
    summarise(weight = sum(abs_rho, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(weight), feature_y)
  
  microbes <- microbes %>%
    mutate(
      x = 0,
      y = rev(seq(1, by = 1.15, length.out = n()))
    )
  
  genes <- genes %>%
    mutate(
      x = 1,
      y = rev(seq(1, by = 1.15, length.out = n()))
    )
  
  edges <- df1 %>%
    left_join(microbes, by = c("feature_x")) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(genes, by = c("feature_y")) %>%
    rename(x2 = x, y2 = y) %>%
    mutate(
      direction2 = factor(direction2, levels = c("Positive", "Negative"))
    )
  
  node_m <- microbes %>%
    transmute(name = feature_x, type = "Genus", x, y)
  
  node_g <- genes %>%
    transmute(name = feature_y, type = "Liver gene", x, y)
  
  nodes <- bind_rows(node_m, node_g)
  
  y_max <- max(c(node_m$y, node_g$y)) + 0.8
  y_min <- min(c(node_m$y, node_g$y)) - 0.8
  
  p <- ggplot() +
    geom_segment(
      data = edges,
      aes(
        x = x1, y = y1,
        xend = x2, yend = y2,
        linewidth = abs_rho,
        color = direction2
      ),
      alpha = 0.72,
      lineend = "round"
    ) +
    geom_point(
      data = subset(nodes, type == "Genus"),
      aes(x = x, y = y),
      shape = 21,
      size = node_size,
      fill = "#E6D6C7",
      color = "black",
      stroke = 0.35,
      show.legend = FALSE
    ) +
    geom_point(
      data = subset(nodes, type == "Liver gene"),
      aes(x = x, y = y),
      shape = 21,
      size = node_size,
      fill = "#CFD9E7",
      color = "black",
      stroke = 0.35,
      show.legend = FALSE
    ) +
    geom_text(
      data = node_m,
      aes(x = -0.18, y = y, label = name),
      family = base_family,
      fontface = "italic",
      hjust = 1,
      size = 2.9
    ) +
    geom_text(
      data = node_g,
      aes(x = 1.18, y = y, label = name),
      family = base_family,
      fontface = "italic",
      hjust = 0,
      size = 3.0
    ) +
    scale_color_manual(values = col_direction, drop = FALSE) +
    scale_linewidth(range = c(0.45, 1.25), guide = "none") +
    coord_cartesian(
      xlim = c(-0.52, 1.52),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = NULL,
      color = "Direction"
    ) +
    theme_void(base_family = base_family) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      plot.margin = margin(6, 34, 6, 34)
    )
  
  tag <- ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = header_label,
      family = base_family,
      size = 4.7
    ) +
    theme_void()
  
  tag / p + plot_layout(heights = c(0.12, 1))
}

p5_left <- build_bipartite_network_plot(
  data_sig = sig_main,
  module_name = "Rum_vs_Liver",
  header_label = "Rum–Liver",
  top_gene_n = 8,
  top_microbe_n = 9,
  max_edge_n = 22,
  node_size = 3.8
)

p5_right <- build_bipartite_network_plot(
  data_sig = sig_main,
  module_name = "Colon_vs_Liver",
  header_label = "Col–Liver",
  top_gene_n = 8,
  top_microbe_n = 9,
  max_edge_n = 22,
  node_size = 3.8
)

p_fig5 <- p5_left | p5_right

save_pdf(
  p_fig5,
  file.path(outdir, "Fig5_pairwise_dual_network.pdf"),
  width = 13.6,
  height = 6.9
)

############################
## 11. Message
############################
message("Final PDF figures have been saved to: ", outdir)
message("Files generated:")
message("1) Fig1_pairwise_module_counts.pdf")
message("2) Fig2_pairwise_connection_spectrum.pdf")
message("3) Fig3_pairwise_liver_gene_bubble.pdf")
message("4) Fig4_pairwise_overlap_functional_map.pdf")
message("5) Fig5_pairwise_dual_network.pdf")