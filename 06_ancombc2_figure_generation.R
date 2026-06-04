############################################################
# Output:
# 01_pattern_phylum_stacked_no_other_with_IT_0.001.pdf
# 02_selected_genus_paired_trajectory_0.001.pdf
# 03_overall_sig_genus_detection_venn_0.001.pdf
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

############################
# 0. packages
############################
pkg_needed <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "tibble", "stringr",
  "forcats", "ggplot2", "scales", "purrr", "grid",
  "ggVennDiagram","VennDiagram", "eulerr"
)

pkg_to_install <- pkg_needed[!pkg_needed %in% installed.packages()[, "Package"]]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, dependencies = TRUE, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

############################
# 1. paths
############################
# Input and output directories
# Please place the required input files in "data/feature_taxa/input".
# Intermediate result tables should be placed in "results/feature_taxa/q0.001".
# Figure files will be saved in "figures/feature_taxa".

result_dir <- file.path("results", "feature_taxa", "q0.001")
prep_dir   <- file.path("data", "feature_taxa", "input")

overall_sig_fp <- file.path(result_dir, "06_Friedman_overall_significant_genus.xlsx")
pairwise_fp    <- file.path(result_dir, "07_pairwise_paired_wilcoxon_for_overall_sig_genus.xlsx")
main_tbl_fp    <- file.path(result_dir, "08_spatial_pattern_classification_main_table.xlsx")

abund_fp <- file.path(prep_dir, "genus_abundance_3group_merged.xlsx")
meta_fp  <- file.path(prep_dir, "metadata_3group.xlsx")

out_dir <- file.path("figures", "feature_taxa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_files <- c(overall_sig_fp, pairwise_fp, main_tbl_fp, abund_fp, meta_fp)
miss_files <- need_files[!file.exists(need_files)]
if (length(miss_files) > 0) {
  stop("The following files do not exist; please check the paths:\n", paste(miss_files, collapse = "\n"))
}

############################
# 2. global settings
############################
q_cutoff_plot <- 0.001
base_family <- "sans"

# group colors
col_rum <- "#7BA8A3"
col_ile <- "#E6C39A"
col_col <- "#C98D7A"

# line colors
col_line_ind <- "#C9C9C9"
col_median   <- "#4F4F4F"

# pattern CN / EN
pattern_order_cn <- c(
  "Foregut-enriched", "Ileal-transition", "Hindgut-enriched",
  "Foregut-hindgut coordinated", "Gradient", "Complex/undetermined"
)

pattern_label_en <- c(
  "Foregut-enriched"   = "Foregut-enriched",
  "Ileal-transition"   = "Ileal-transition",
  "Hindgut-enriched"   = "Hindgut-enriched",
  "Foregut-hindgut coordinated" = "Foregut-hindgut coordinated",
  "Gradient"       = "Gradient",
  "Complex/undetermined" = "Complex/other"
)

pattern_order_en <- c(
  "Foregut-enriched",
  "Ileal-transition",
  "Hindgut-enriched",
  "Foregut-hindgut coordinated",
  "Gradient",
  "Complex/other"
)

# Fig2 fixed 8 genera
traj_genus_tbl <- tibble::tribble(
  ~Genus,             ~Pattern,
  "Prevotella_7",     "Foregut-enriched",
  "Butyrivibrio",     "Foregut-enriched",
  "Fournierella",     "Hindgut-enriched",
  "Treponema",        "Hindgut-enriched",
  "Anaerovibrio",     "Foregut-hindgut coordinated",
  "Ruminococcus",     "Foregut-hindgut coordinated",
  "Succiniclasticum", "Gradient",
  "Romboutsia",       "Gradient"
)

pattern_abbr_map <- c(
  "Foregut-enriched"   = "FGH",
  "Ileal-transition"   = "IT",
  "Hindgut-enriched"   = "HGH",
  "Foregut-hindgut coordinated" = "FHC",
  "Gradient"       = "Grad",
  "Complex/undetermined" = "Other"
)

# Fig. 1 phylum colors (excluding Other)
phylum_cols <- c(
  "Firmicutes" = "#8FB9A8",
  "Bacteroidota" = "#D8B77E",
  "Proteobacteria" = "#C98D7A",
  "Actinobacteriota" = "#9C8FB8",
  "Spirochaetota" = "#8FA7C6"
)

############################
# 3. helper functions
############################
theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "plain", color = "black"),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

save_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_obj,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

facet_label_expr <- function(x) {
  vapply(x, function(s) {
    genus <- sub(" \\(.*\\)$", "", s)
    abbr  <- sub("^.*\\(", "", sub("\\)$", "", s))
    genus <- gsub("'", "", genus)
    paste0("italic('", genus, "')~' (", abbr, ")'")
  }, character(1))
}

############################
# 4. read data
############################
overall_sig_tbl <- readxl::read_excel(overall_sig_fp, sheet = 1)
pairwise_tbl    <- readxl::read_excel(pairwise_fp, sheet = 1)
main_tbl        <- readxl::read_excel(main_tbl_fp, sheet = 1)
abund_df        <- readxl::read_excel(abund_fp, sheet = 1)
meta_df         <- readxl::read_excel(meta_fp, sheet = 1)

colnames(abund_df)[1] <- "Genus"
colnames(meta_df) <- c("SampleID", "Group", "SheepID")

meta_df <- meta_df %>%
  mutate(
    SampleID = as.character(SampleID),
    SheepID  = as.character(SheepID),
    Group    = factor(Group, levels = c("Rum", "Ile", "Col"))
  )

############################
# 5. main table
############################
main_tbl_001 <- main_tbl %>%
  filter(!is.na(q_friedman), q_friedman <= q_cutoff_plot) %>%
  left_join(
    pairwise_tbl %>% dplyr::select(any_of(c("Genus", "p_RI", "p_RC", "p_IC", "q_RI", "q_RC", "q_IC", "sig_RI", "sig_RC", "sig_IC"))),
    by = "Genus",
    suffix = c("", "_pair")
  ) %>%
  mutate(
    pattern_class_cn = dplyr::case_when(
      is.na(pattern_class) ~ "Complex/undetermined",
      TRUE ~ as.character(pattern_class)
    ),
    pattern_class_en = unname(pattern_label_en[pattern_class_cn]),
    pattern_class_en = factor(pattern_class_en, levels = pattern_order_en),
    med_Rum = safe_num(med_Rum),
    med_Ile = safe_num(med_Ile),
    med_Col = safe_num(med_Col),
    Phylum = as.character(Phylum),
    Phylum = ifelse(is.na(Phylum) | Phylum == "", NA_character_, Phylum)
  )

sig_genus_001 <- unique(main_tbl_001$Genus)

traj_missing <- setdiff(traj_genus_tbl$Genus, sig_genus_001)
if (length(traj_missing) > 0) {
  warning("The following Fig. 2 candidate genera were not included in the q < 0.001 main analysis set:\n", paste(traj_missing, collapse = ", "))
}

traj_genus_use <- traj_genus_tbl %>%
  filter(Genus %in% sig_genus_001) %>%
  mutate(
    Pattern = factor(Pattern, levels = pattern_order_cn),
    Pattern_abbr = unname(pattern_abbr_map[as.character(Pattern)])
  ) %>%
  arrange(Pattern, Genus)

############################
# 6. abundance for Fig2 / Fig3
############################
sample_cols <- intersect(colnames(abund_df), meta_df$SampleID)
if (length(sample_cols) == 0) {
  stop("No sample columns matched between the abundance table and metadata; please check SampleID.")
}

abund_df2 <- abund_df %>%
  dplyr::select(Genus, all_of(sample_cols))

otu_mat <- abund_df2 %>%
  tibble::column_to_rownames("Genus") %>%
  as.matrix()
storage.mode(otu_mat) <- "numeric"

rel_mat <- sweep(otu_mat, 2, colSums(otu_mat), "/")
rel_mat[is.na(rel_mat)] <- 0

rel_long <- as.data.frame(rel_mat) %>%
  tibble::rownames_to_column("Genus") %>%
  tidyr::pivot_longer(-Genus, names_to = "SampleID", values_to = "RelAbund") %>%
  left_join(meta_df, by = "SampleID") %>%
  mutate(
    Group = factor(Group, levels = c("Rum", "Ile", "Col")),
    plot_value = log10(RelAbund * 1e6 + 1)
  )

############################
# 7. Fig. 1: pattern-by-phylum stacked bar (including IT and excluding Other)
############################
fig1_use_patterns <- c(
  "Foregut-enriched",
  "Ileal-transition",
  "Hindgut-enriched",
  "Foregut-hindgut coordinated",
  "Gradient"
)

fig1_df <- main_tbl_001 %>%
  filter(pattern_class_en %in% fig1_use_patterns) %>%
  mutate(
    Phylum_plot = dplyr::case_when(
      Phylum %in% names(phylum_cols) ~ Phylum,
      TRUE ~ NA_character_
    ),
    pattern_class_en = factor(pattern_class_en, levels = fig1_use_patterns)
  ) %>%
  filter(!is.na(Phylum_plot)) %>%
  count(pattern_class_en, Phylum_plot, name = "n")

p_fig1 <- ggplot(fig1_df, aes(x = pattern_class_en, y = n, fill = Phylum_plot)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = phylum_cols, drop = FALSE) +
  labs(x = NULL, y = "Number of significant genera") +
  theme_paper(12.5) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11.5, color = "black"),
    legend.key.size = grid::unit(0.95, "lines"),
    axis.text.x = element_text(size = 10.8, angle = 15, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12.5)
  )

save_pdf(p_fig1, "01_pattern_phylum_stacked_no_other_with_IT_0.001.pdf", width = 9.8, height = 6.8)

############################
# 8. Fig2: 8 genera paired trajectory
############################
traj_df <- rel_long %>%
  filter(Genus %in% traj_genus_use$Genus) %>%
  left_join(traj_genus_use, by = "Genus") %>%
  mutate(
    Pattern = factor(Pattern, levels = pattern_order_cn),
    Genus = factor(Genus, levels = traj_genus_use$Genus),
    facet_lab_chr = paste0(as.character(Genus), " (", Pattern_abbr, ")")
  )

traj_med <- traj_df %>%
  group_by(Genus, facet_lab_chr, Pattern, Group) %>%
  summarise(
    med = median(plot_value, na.rm = TRUE),
    .groups = "drop"
  )

p_fig2 <- ggplot(traj_df, aes(x = Group, y = plot_value, group = SheepID)) +
  geom_line(color = col_line_ind, linewidth = 0.35, alpha = 0.82) +
  geom_point(
    aes(fill = Group),
    shape = 21, color = "black", stroke = 0.25,
    size = 1.9, alpha = 0.96
  ) +
  geom_line(
    data = traj_med,
    aes(x = Group, y = med, group = 1),
    inherit.aes = FALSE,
    color = col_median,
    linewidth = 0.95
  ) +
  geom_point(
    data = traj_med,
    aes(x = Group, y = med),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = col_median,
    stroke = 0.60,
    size = 2.3
  ) +
  scale_fill_manual(values = c("Rum" = col_rum, "Ile" = col_ile, "Col" = col_col)) +
  labs(x = NULL, y = expression(log[10]("relative abundance" %*% 10^6 + 1))) +
  facet_wrap(
    ~ facet_lab_chr,
    ncol = 4,
    scales = "free_y",
    labeller = as_labeller(facet_label_expr, label_parsed)
  ) +
  theme_paper(11.8) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 15, color = "black"),
    legend.key.size = grid::unit(1.25, "lines"),
    strip.text = element_text(size = 10.8, lineheight = 0.92),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_text(size = 10.2, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3.8, shape = 21),
      nrow = 1,
      byrow = TRUE
    )
  )

save_pdf(p_fig2, "02_selected_genus_paired_trajectory_0.001.pdf", width = 13.4, height = 7.6)

############################
# 9. Fig3: Venn of overall significant genera by detection
# Definition:
# R set = genera among overall significant taxa with abundance > 0 in any Rum sample
# I set = genera among overall significant taxa with abundance > 0 in any Ile sample
# C set = genera among overall significant taxa with abundance > 0 in any Col sample
# Therefore:
# R-only = detected only in Rum and not detected in Ile or Col
############################
overall_abund_long <- abund_df2 %>%
  filter(Genus %in% sig_genus_001) %>%
  tidyr::pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(meta_df, by = "SampleID") %>%
  mutate(
    Abundance = safe_num(Abundance)
  )

detect_tbl <- overall_abund_long %>%
  group_by(Genus, Group) %>%
  summarise(
    detected = any(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = detected,
    values_fill = FALSE
  )

set_R <- detect_tbl %>% filter(Rum) %>% pull(Genus) %>% unique()
set_I <- detect_tbl %>% filter(Ile) %>% pull(Genus) %>% unique()
set_C <- detect_tbl %>% filter(Col) %>% pull(Genus) %>% unique()

# Colors: soft palette
venn_fill_cols <- c("#7BA8A3", "#E6C39A", "#C98D7A")

# Aspect ratio: close to square for a more balanced layout
venn_width  <- 7.2
venn_height <- 6.6

venn_grob <- VennDiagram::venn.diagram(
  x = list(
    Rum = set_R,
    Ile = set_I,
    Col = set_C
  ),
  filename = NULL,
  imagetype = "pdf",
  fill = venn_fill_cols,
  alpha = c(0.45, 0.45, 0.45),
  col = c("#5F8F89", "#D2AE78", "#B77766"),
  lwd = 1.2,
  cex = 1.8,                # intersection-number size
  fontfamily = base_family,
  fontface = "plain",
  cat.cex = 1.6,            # set-name size
  cat.fontfamily = base_family,
  cat.fontface = "plain",
  cat.dist = c(0.055, 0.055, 0.055),
  cat.pos = c(-20, 20, 180),
  margin = 0.08
)

pdf(
  file = file.path(out_dir, "03_overall_sig_genus_detection_venn_0.001.pdf"),
  width = venn_width,
  height = venn_height,
  family = base_family
)
grid::grid.newpage()
grid::grid.draw(venn_grob)
dev.off()

############################
# 11. console message
############################
cat("\n================ Completed =================\n")
cat("Output directory: ", out_dir, "\n")
cat("Output files:\n")
cat("01_pattern_phylum_stacked_no_other_with_IT_0.001.pdf\n")
cat("02_selected_genus_paired_trajectory_0.001.pdf\n")
cat("03_overall_sig_genus_detection_venn_0.001.pdf\n")
cat("Fig1_Fig2_Fig3_helper_tables.xlsx\n")
cat("========================================\n")