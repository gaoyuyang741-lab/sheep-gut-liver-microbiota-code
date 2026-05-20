############################################################
# 输出：
#   Fig1_pattern_ring_KO_Pathway_L2.pdf
#   Fig2_filtered_L2_function_division.pdf
#   Fig4_representative_pathway_pattern_trajectories.pdf
#   Fig5_summary_sankey_RCsource_fixed.pdf
#   PICRUSt2_rebuilt_plot_tables.xlsx
############################################################

rm(list = ls())
gc()
graphics.off()

############################
# 0. 加载程序包
############################
need_pkgs <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "stringr", "forcats",
  "ggplot2", "scales", "showtext", "sysfonts", "ggalluvial"
)

to_install <- need_pkgs[!need_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(ggalluvial)
})

options(stringsAsFactors = FALSE)

############################
# 1. 路径设置
############################
find_file <- function(primary, fallback = NULL) {
  if (!is.null(primary) && file.exists(primary)) return(primary)
  if (!is.null(fallback) && file.exists(fallback)) return(fallback)
  stop("File not found:\n  ", primary,
       if (!is.null(fallback)) paste0("\n  ", fallback) else "")
}

fp_picrust <- find_file(
  file.path("results", "picrust2", "PICRUSt2_核心结果_合并.xlsx"),
  "PICRUSt2_核心结果_合并.xlsx"
)

fp_genus <- find_file(
  file.path("results", "feature_taxa", "q0.001", "08_spatial_pattern_classification_main_table.xlsx"),
  "08_spatial_pattern_classification_main_table.xlsx"
)

outdir <- file.path("figures", "picrust2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################
# 2. 字体与保存函数
############################
font_family <- "sans"
font_candidates <- c(
  "C:/Windows/Fonts/arial.ttf",
  "C:/Windows/Fonts/Arial.ttf",
  "C:/Windows/Fonts/msyh.ttc",
  "C:/Windows/Fonts/msyh.ttf"
)
font_file <- font_candidates[file.exists(font_candidates)][1]

if (!is.na(font_file) && length(font_file) > 0) {
  try({
    sysfonts::font_add(family = "figfont", regular = font_file)
    showtext::showtext_auto()
    font_family <- "figfont"
  }, silent = TRUE)
}

pdf_family <- "Arial"
if (.Platform$OS.type == "windows") {
  suppressWarnings({
    windowsFonts(Arial = windowsFont("Arial"))
  })
}

save_pdf_editable <- function(plot_obj, filename, width = 10, height = 8, family = pdf_family) {
  showtext::showtext_auto(FALSE)
  grDevices::cairo_pdf(filename = filename, width = width, height = height, family = family)
  print(plot_obj)
  dev.off()
  showtext::showtext_auto(TRUE)
}

theme_sci <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = "black"),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", color = "black", size = 12),
      axis.text = element_text(color = "black", size = 11),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", color = "black", size = 11),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 11),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(10, 14, 10, 14)
    )
}

############################
# 3. 配色
############################
col_fgh   <- "#C58A8A"
col_fhc   <- "#B8A18C"
col_grad  <- "#A9B6A2"
col_hgh   <- "#8FAFC3"
col_it    <- "#C7B2D6"
col_other <- "#D8D8D8"

pattern_palette <- c(
  "FGH"   = col_fgh,
  "FHC"   = col_fhc,
  "Grad"  = col_grad,
  "HGH"   = col_hgh,
  "IT"    = col_it,
  "Other" = col_other
)

site_palette <- c(
  "Rum" = "#8DAEAA",
  "Ile" = "#C7B2D6",
  "Col" = "#A9B6D3"
)

source_cols <- c(
  "R source" = "#C98E8E",
  "C source" = "#B9A38E"
)

############################
# 4. 辅助函数
############################
to_pattern_en <- function(x) {
  dplyr::case_when(
    x %in% c("前肠偏高型", "FGH") ~ "FGH",
    x %in% c("前后肠协同型", "FHC") ~ "FHC",
    x %in% c("梯度型", "Grad") ~ "Grad",
    x %in% c("后肠偏高型", "HGH") ~ "HGH",
    x %in% c("回肠过渡型", "Ile transition", "IT") ~ "IT",
    TRUE ~ as.character(x)
  )
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 1e-4, format(x, scientific = TRUE, digits = 2), sprintf("%.4f", x))
  )
}

rescale01 <- function(x, to = c(0.3, 1)) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(mean(to), length(x)))
  (x - rng[1]) / diff(rng) * diff(to) + to[1]
}

guess_first <- function(cands, cn) {
  hit <- cands[cands %in% cn]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

score_from_table <- function(df, median_cols, q_col, prev_col = NULL) {
  med_max <- apply(df[, median_cols, drop = FALSE], 1, max, na.rm = TRUE)
  qv <- pmax(df[[q_col]], 1e-300)
  if (is.null(prev_col) || !prev_col %in% colnames(df)) {
    prev <- rep(1, nrow(df))
  } else {
    prev <- pmax(df[[prev_col]], 1e-6)
  }
  raw <- (-log10(qv)) * (med_max + 1e-12) * prev
  rescale01(raw, to = c(0.35, 1))
}

clean_label <- function(x, width = 24) {
  x %>%
    str_replace_all("_", " ") %>%
    str_wrap(width = width)
}

scale_within_pathway <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

shorten_pathway_fixed <- function(x) {
  map <- c(
    "Biosynthesis of nucleotide sugars" = "Biosynthesis of\nnucleotide sugars",
    "Nicotinate and nicotinamide metabolism" = "Nicotinate and\nnicotinamide metab.",
    "Protein digestion and absorption" = "Protein digestion\nand absorption",
    "N-Glycan biosynthesis" = "N-Glycan\nbiosynth.",
    "Alanine, aspartate and glutamate metabolism" = "Alanine, aspartate\nand glutamate metab."
  )
  y <- unname(map[x])
  y[is.na(y)] <- x[is.na(y)]
  y
}

wrap_two_lines <- function(x, width = 24) {
  x <- str_replace_all(x, "_", " ")
  x <- str_wrap(x, width = width)
  vapply(
    strsplit(x, "\n", fixed = TRUE),
    function(z) {
      if (length(z) <= 2) {
        paste(z, collapse = "\n")
      } else {
        paste(c(z[1], paste(z[2:length(z)], collapse = " ")), collapse = "\n")
      }
    },
    character(1)
  )
}

############################
# 5. 读入数据
############################
sheet_names <- excel_sheets(fp_picrust)
req_sheets <- c("summary_pattern", "pathway_all_results", "pathway_sig_results", "KO_all_results")
missing_sheets <- setdiff(req_sheets, sheet_names)
if (length(missing_sheets) > 0) {
  stop("PICRUSt2 结果文件缺少 sheet: ", paste(missing_sheets, collapse = ", "))
}

sum_pattern <- read_excel(fp_picrust, sheet = "summary_pattern")
path_all    <- read_excel(fp_picrust, sheet = "pathway_all_results")
path_sig    <- read_excel(fp_picrust, sheet = "pathway_sig_results")
ko_all      <- read_excel(fp_picrust, sheet = "KO_all_results")
genus_main  <- read_excel(fp_genus)

sum_pattern$pattern_class <- to_pattern_en(sum_pattern$pattern_class)
path_all$pattern_class    <- to_pattern_en(path_all$pattern_class)
path_sig$pattern_class    <- to_pattern_en(path_sig$pattern_class)
ko_all$pattern_class      <- to_pattern_en(ko_all$pattern_class)
if ("pattern_class" %in% colnames(genus_main)) {
  genus_main$pattern_en <- to_pattern_en(genus_main$pattern_class)
}

############################
# 6. FIGURE 1
############################
fig1_levels_keep <- c("KO", "Pathway", "L2")

fig1_df <- sum_pattern %>%
  filter(
    level %in% fig1_levels_keep,
    !is.na(pattern_class),
    !pattern_class %in% c("Other", "NA", "")
  ) %>%
  mutate(
    pattern_class = factor(pattern_class, levels = c("FGH", "FHC", "Grad", "HGH", "IT")),
    level = factor(level, levels = fig1_levels_keep)
  ) %>%
  group_by(level) %>%
  mutate(
    pct = n / sum(n),
    ymax = cumsum(pct),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2
  ) %>%
  ungroup()

fig1_center <- fig1_df %>%
  group_by(level) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  mutate(center_lab = paste0(as.character(level), "\n(n = ", total_n, ")"))

p1 <- ggplot(fig1_df, aes(x = 2, y = pct, fill = pattern_class)) +
  geom_col(color = "white", linewidth = 0.8, width = 1) +
  coord_polar(theta = "y") +
  xlim(0.7, 2.5) +
  facet_wrap(~ level, nrow = 1) +
  geom_text(
    aes(y = ymid, label = ifelse(pct >= 0.08, paste0(pattern_class, "\n", percent(pct, accuracy = 0.1)), "")),
    size = 3.5, family = font_family, lineheight = 0.95
  ) +
  geom_text(
    data = fig1_center,
    aes(x = 0.95, y = 0, label = center_lab),
    inherit.aes = FALSE,
    family = font_family, fontface = "bold", size = 4.1
  ) +
  scale_fill_manual(values = pattern_palette, drop = FALSE) +
  theme_void(base_family = font_family) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11),
    strip.text = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

save_pdf_editable(
  p1,
  file.path(outdir, "Fig1_pattern_ring_KO_Pathway_L2.pdf"),
  width = 12.0, height = 4.8
)

############################
# 7. FIGURE 2
############################
drop_l2 <- c(
  "Cancer: overview",
  "Cancer: specific types",
  "Infectious disease: bacterial",
  "Infectious disease: viral",
  "Neurodegenerative disease",
  "Nervous system",
  "Substance dependence",
  "Drug resistance: antimicrobial",
  "Drug resistance: antineoplastic",
  "Cardiovascular diseases",
  "Aging",
  "Sensory system"
)

keep_priority <- c(
  "Global and overview maps",
  "Carbohydrate metabolism",
  "Amino acid metabolism",
  "Lipid metabolism",
  "Glycan biosynthesis and metabolism",
  "Metabolism of cofactors and vitamins",
  "Metabolism of terpenoids and polyketides",
  "Biosynthesis of other secondary metabolites",
  "Xenobiotics biodegradation and metabolism",
  "Energy metabolism",
  "Digestive system",
  "Endocrine system",
  "Immune system",
  "Signal transduction",
  "Cell growth and death",
  "Endocrine and metabolic disease"
)

fig2_df <- path_sig %>%
  filter(
    !is.na(pattern_class),
    !pattern_class %in% c("NA", "", "Other"),
    !is.na(L2),
    !L2 %in% drop_l2,
    L2 %in% keep_priority
  ) %>%
  count(L2, pattern_class, name = "n_path") %>%
  complete(
    L2 = keep_priority,
    pattern_class = c("FGH", "FHC", "Grad", "HGH", "IT"),
    fill = list(n_path = 0)
  ) %>%
  group_by(L2) %>%
  mutate(total = sum(n_path)) %>%
  ungroup() %>%
  filter(total > 0) %>%
  mutate(
    L2 = factor(L2, levels = rev(keep_priority[keep_priority %in% unique(L2)])),
    pattern_class = factor(pattern_class, levels = c("FGH", "FHC", "Grad", "HGH", "IT"))
  )

tot_lab <- fig2_df %>%
  distinct(L2, total)

p2 <- ggplot(fig2_df, aes(x = n_path, y = L2, fill = pattern_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.5) +
  geom_text(
    data = tot_lab,
    aes(x = total + 0.35, y = L2, label = total),
    inherit.aes = FALSE,
    hjust = 0, size = 4.0, family = font_family, fontface = "bold"
  ) +
  scale_fill_manual(values = pattern_palette[c("FGH", "FHC", "Grad", "HGH", "IT")], drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of significant pathways", y = NULL, fill = "Pattern") +
  theme_sci(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10.4)
  )

save_pdf_editable(
  p2,
  file.path(outdir, "Fig2_filtered_L2_function_division.pdf"),
  width = 9.6, height = 7.2
)

############################
# 8. FIGURE 4
############################
selected_pathways <- c(
  "Inositol phosphate metabolism",
  "Biosynthesis of nucleotide sugars",
  "Protein digestion and absorption",
  "N-Glycan biosynthesis",
  "Pentose phosphate pathway",
  "Alanine, aspartate and glutamate metabolism"
)

fig4_base <- path_all %>%
  filter(pathway_name %in% selected_pathways) %>%
  mutate(pathway_name = factor(pathway_name, levels = selected_pathways))

if (nrow(fig4_base) > 0) {
  fig4_df <- fig4_base %>%
    transmute(
      pathway_name,
      pattern_class = factor(pattern_class, levels = c("FGH", "FHC", "Grad", "HGH", "IT")),
      overall_q,
      Rum = median_Rum,
      Ile = median_Ile,
      Col = median_Col
    ) %>%
    pivot_longer(cols = c(Rum, Ile, Col), names_to = "Site", values_to = "MedianAbundance") %>%
    group_by(pathway_name) %>%
    mutate(z = scale_within_pathway(MedianAbundance)) %>%
    ungroup() %>%
    mutate(
      Site = factor(Site, levels = c("Rum", "Ile", "Col")),
      facet_lab = paste0("[", pattern_class, "] ", pathway_name, "\nq = ", fmt_p(overall_q))
    )
  
  facet_levels <- fig4_df %>%
    distinct(pathway_name, facet_lab) %>%
    arrange(pathway_name) %>%
    pull(facet_lab)
  
  fig4_df$facet_lab <- factor(fig4_df$facet_lab, levels = facet_levels)
  
  p4 <- ggplot(fig4_df, aes(x = Site, y = z, group = pathway_name)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5, linetype = "dashed") +
    geom_line(linewidth = 0.9, color = "grey45") +
    geom_point(aes(fill = Site), shape = 21, size = 3.2, color = "black", stroke = 0.35) +
    facet_wrap(~ facet_lab, ncol = 3, scales = "fixed") +
    scale_fill_manual(values = site_palette) +
    labs(x = NULL, y = "Within-pathway standardized abundance (z-score)") +
    theme_sci(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10.1),
      axis.text.x = element_text(face = "bold", size = 12)
    )
  
  save_pdf_editable(
    p4,
    file.path(outdir, "Fig4_representative_pathway_pattern_trajectories.pdf"),
    width = 10.8, height = 7.4
  )
} else {
  warning("Fig4 未匹配到任何 pathway，请检查 pathway_name。")
}

############################
# 9. FIGURE 5
# 强制固定 5 个 pathway，不再让绘图自己推断
############################
genus_nodes <- tibble::tribble(
  ~Genus,                   ~Genus_source,
  "Butyrivibrio",           "R source",
  "Paraprevotella",         "R source",
  "Ruminococcus",           "R source",
  "Anaerovibrio",           "R source",
  "Treponema",              "C source",
  "Oscillospira",           "C source",
  "Colidextribacter",       "C source",
  "Prevotellaceae_UCG_003", "C source"
)

edge_genus_ko <- tibble::tribble(
  ~Genus,                   ~KO,
  "Butyrivibrio",           "K18677",
  "Butyrivibrio",           "K14974",
  "Paraprevotella",         "K18030",
  "Ruminococcus",           "K01278",
  "Ruminococcus",           "K00721",
  "Anaerovibrio",           "K01278",
  "Treponema",              "K00850",
  "Treponema",              "K01915",
  "Oscillospira",           "K01915",
  "Colidextribacter",       "K00850",
  "Prevotellaceae_UCG_003", "K00721"
)

edge_ko_path <- tibble::tribble(
  ~KO,       ~Pathway,
  "K18677",  "Biosynthesis of nucleotide sugars",
  "K14974",  "Nicotinate and nicotinamide metabolism",
  "K18030",  "Nicotinate and nicotinamide metabolism",
  "K01278",  "Protein digestion and absorption",
  "K00721",  "N-Glycan biosynthesis",
  "K00850",  "Alanine, aspartate and glutamate metabolism",
  "K01915",  "Alanine, aspartate and glutamate metabolism"
)

pathway_keep_fixed <- c(
  "Biosynthesis of nucleotide sugars",
  "Nicotinate and nicotinamide metabolism",
  "Protein digestion and absorption",
  "N-Glycan biosynthesis",
  "Alanine, aspartate and glutamate metabolism"
)

genus_cn <- colnames(genus_main)
genus_col  <- guess_first(c("Genus", "genus", "feature", "taxon"), genus_cn)
g_q_col    <- guess_first(c("q_friedman", "overall_q", "q", "p_adj"), genus_cn)
g_prev_col <- guess_first(c("prevalence_all", "prevalence", "prev"), genus_cn)
g_med_r    <- guess_first(c("med_Rum", "median_Rum", "Rum_median"), genus_cn)
g_med_i    <- guess_first(c("med_Ile", "median_Ile", "Ile_median"), genus_cn)
g_med_c    <- guess_first(c("med_Col", "median_Col", "Col_median"), genus_cn)

if (is.na(genus_col) || is.na(g_q_col) || is.na(g_med_r) || is.na(g_med_i) || is.na(g_med_c)) {
  stop("08_spatial_pattern_classification_main_table.xlsx 缺少必要列，请检查 genus / q / median_Rum/Ile/Col 列名。")
}

genus_sel <- genus_main %>%
  rename(Genus = !!genus_col) %>%
  filter(Genus %in% genus_nodes$Genus) %>%
  mutate(
    score = score_from_table(
      df = .,
      median_cols = c(g_med_r, g_med_i, g_med_c),
      q_col = g_q_col,
      prev_col = if (!is.na(g_prev_col)) g_prev_col else NULL
    )
  ) %>%
  select(Genus, score) %>%
  left_join(genus_nodes, by = "Genus")

ko_sel <- ko_all %>%
  filter(feature %in% unique(c(edge_genus_ko$KO, edge_ko_path$KO))) %>%
  mutate(
    score = score_from_table(
      df = .,
      median_cols = c("median_Rum", "median_Ile", "median_Col"),
      q_col = "overall_q",
      prev_col = if ("prevalence" %in% colnames(.)) "prevalence" else NULL
    )
  ) %>%
  transmute(KO = feature, score)

path_sel <- path_all %>%
  filter(pathway_name %in% pathway_keep_fixed) %>%
  mutate(
    score = score_from_table(
      df = .,
      median_cols = c("median_Rum", "median_Ile", "median_Col"),
      q_col = "overall_q",
      prev_col = if ("prevalence" %in% colnames(.)) "prevalence" else NULL
    )
  ) %>%
  transmute(Pathway = pathway_name, score)

fig5_df <- edge_genus_ko %>%
  left_join(genus_sel %>% rename(genus_score = score), by = "Genus") %>%
  left_join(ko_sel %>% rename(ko_score = score), by = "KO") %>%
  inner_join(edge_ko_path, by = "KO") %>%
  filter(Pathway %in% pathway_keep_fixed) %>%
  left_join(path_sel %>% rename(path_score = score), by = "Pathway") %>%
  mutate(
    Genus_source = genus_nodes$Genus_source[match(Genus, genus_nodes$Genus)],
    genus_score = ifelse(is.na(genus_score), 0.5, genus_score),
    ko_score    = ifelse(is.na(ko_score), 0.5, ko_score),
    path_score  = ifelse(is.na(path_score), 0.5, path_score),
    weight_raw  = sqrt(sqrt(genus_score * ko_score) * sqrt(ko_score * path_score)),
    weight      = rescale01(weight_raw, to = c(0.8, 2.8)),
    Genus_source = factor(Genus_source, levels = c("R source", "C source")),
    Genus_show = wrap_two_lines(Genus, width = 22),
    KO_show = factor(KO, levels = unique(edge_ko_path$KO)),
    Pathway_show_raw = Pathway,
    Pathway_show = shorten_pathway_fixed(Pathway)
  ) %>%
  select(Genus_source, Genus_show, KO_show, Pathway_show_raw, Pathway_show, weight)

# 强制固定 5 个右侧 pathway 顺序
pathway_show_levels <- shorten_pathway_fixed(pathway_keep_fixed)

# 排序
genus_order <- fig5_df %>%
  group_by(Genus_source, Genus_show) %>%
  summarise(w = sum(weight), .groups = "drop") %>%
  arrange(Genus_source, desc(w), Genus_show) %>%
  pull(Genus_show) %>%
  unique()

ko_order <- fig5_df %>%
  group_by(KO_show) %>%
  summarise(w = sum(weight), .groups = "drop") %>%
  arrange(desc(w), KO_show) %>%
  pull(KO_show) %>%
  as.character()

fig5_df <- fig5_df %>%
  mutate(
    Genus_show = factor(Genus_show, levels = genus_order),
    KO_show = factor(as.character(KO_show), levels = ko_order),
    Pathway_show = factor(Pathway_show, levels = pathway_show_levels)
  )

p5 <- ggplot(
  fig5_df,
  aes(
    axis1 = Genus_show,
    axis2 = KO_show,
    axis3 = Pathway_show,
    y = weight
  )
) +
  geom_alluvium(
    aes(fill = Genus_source),
    width = 0.13,
    knot.pos = 0.42,
    alpha = 0.84,
    color = "white",
    linewidth = 0.25
  ) +
  geom_stratum(
    width = 0.32,
    fill = "#EEE8DF",
    color = "grey45",
    linewidth = 0.40
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    family = font_family,
    size = 12 / 2.845,
    fontface = "bold",
    lineheight = 0.92
  ) +
  scale_fill_manual(values = source_cols, drop = FALSE) +
  scale_x_discrete(
    limits = c("Genus", "KO", "Pathway"),
    expand = c(0.06, 0.02)
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 11, color = "black"),
    plot.margin = margin(14, 18, 14, 18)
  )

save_pdf_editable(
  p5,
  file.path(outdir, "Fig5_summary_sankey_RCsource_fixed.pdf"),
  width = 12.4,
  height = 7.2
)

############################
# 10. 导出核查表
############################
wb <- createWorkbook()

addWorksheet(wb, "Fig1_pattern_ring")
writeData(wb, "Fig1_pattern_ring", fig1_df %>% select(level, pattern_class, n, pct))

addWorksheet(wb, "Fig2_filtered_L2")
writeData(wb, "Fig2_filtered_L2", fig2_df)

if (exists("fig4_base") && nrow(fig4_base) > 0) {
  addWorksheet(wb, "Fig4_selected_pathways")
  writeData(
    wb, "Fig4_selected_pathways",
    path_all %>%
      filter(pathway_name %in% selected_pathways) %>%
      select(pathway_name, pattern_class, overall_q, median_Rum, median_Ile, median_Col, L1, L2)
  )
}

addWorksheet(wb, "Fig5_genus_nodes")
writeData(wb, "Fig5_genus_nodes", genus_nodes)

addWorksheet(wb, "Fig5_edge_genus_ko")
writeData(wb, "Fig5_edge_genus_ko", edge_genus_ko)

addWorksheet(wb, "Fig5_edge_ko_path")
writeData(wb, "Fig5_edge_ko_path", edge_ko_path)

addWorksheet(wb, "Fig5_plot_df")
writeData(wb, "Fig5_plot_df", fig5_df)

saveWorkbook(
  wb,
  file.path(outdir, "PICRUSt2_rebuilt_plot_tables.xlsx"),
  overwrite = TRUE
)

############################
# 11. 完成提示
############################
message("全部完成。输出目录：", normalizePath(outdir, winslash = "/"))
message("已输出：")
message("  Fig1_pattern_ring_KO_Pathway_L2.pdf")
message("  Fig2_filtered_L2_function_division.pdf")
if (exists("fig4_base") && nrow(fig4_base) > 0) {
  message("  Fig4_representative_pathway_pattern_trajectories.pdf")
}
message("  Fig5_summary_sankey_RCsource_fixed.pdf")
message("  PICRUSt2_rebuilt_plot_tables.xlsx")