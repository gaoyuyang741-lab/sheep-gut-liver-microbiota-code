############################################################
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

## =========================
## 0. 环境准备
## =========================
my_tmp <- "D:/R_tmp"
dir.create(my_tmp, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = my_tmp)
Sys.setenv(TMP = my_tmp)
Sys.setenv(TEMP = my_tmp)

pkgs <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "ggplot2",
  "stringr", "forcats", "scales", "patchwork",
  "zCompositions", "compositions", "tibble", "grid"
)

need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(need) > 0) install.packages(need, dependencies = TRUE)

library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)
library(patchwork)
library(zCompositions)
library(compositions)
library(tibble)
library(grid)

## =========================
## 1. 路径设置
## =========================
# Input and output directories
# This script reads coordinated-genera result files from "results/coordinated_genera".
# Required abundance and metadata files should be placed in "data/coordinated_genera/input".
# Figure files will be saved in "figures/coordinated_genera".

res_fp <- file.path(
  "results", "coordinated_genera",
  "02_连续性结构分析_结果文件.xlsx"
)

genus_fp <- file.path(
  "data", "coordinated_genera", "input",
  "genus_abundance_3group_merged.xlsx"
)

meta_fp <- file.path(
  "data", "coordinated_genera", "input",
  "metadata_3group.xlsx"
)

out_dir <- file.path("figures", "coordinated_genera")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(res_fp))   stop("结果文件不存在：", res_fp)
if (!file.exists(genus_fp)) stop("genus 丰度表不存在：", genus_fp)
if (!file.exists(meta_fp))  stop("metadata 文件不存在：", meta_fp)

## =========================
## 2. 配色（温和）
## =========================
col_tcg_fill <- "#D9CCE8"
col_tcg_line <- "#A68DBE"

col_fcg_fill <- "#C9D8EE"
col_fcg_line <- "#7F9CCB"

col_elg_fill <- "#F3D8C2"
col_elg_line <- "#D39A6A"

col_hcg_fill <- "#CFE4E2"
col_hcg_line <- "#6FA6A0"

pair_line_map <- c(
  "Rumen–Ileum" = "#7F9CCB",
  "Rumen–Colon" = "#D39A6A",
  "Ileum–Colon" = "#6FA6A0"
)

pair_fill_map <- c(
  "Rumen–Ileum" = "#C9D8EE",
  "Rumen–Colon" = "#F3D8C2",
  "Ileum–Colon" = "#CFE4E2"
)

class_fill_map <- c(
  "TCG" = col_tcg_fill,
  "FCG" = col_fcg_fill,
  "ELG" = col_elg_fill,
  "HCG" = col_hcg_fill
)

col_text   <- "#333333"
col_border <- "#9A9A9A"
col_grid   <- "#ECECEC"
col_na     <- "#F5F5F5"

## =========================
## 3. 通用函数
## =========================
theme_sci_soft <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(colour = col_grid, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = col_border, linewidth = 0.6),
      axis.line = element_line(colour = col_border, linewidth = 0.4),
      axis.ticks = element_line(colour = col_border, linewidth = 0.4),
      axis.text = element_text(colour = col_text),
      axis.title = element_text(colour = col_text, face = "plain"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(colour = col_text),
      strip.background = element_rect(fill = "#F5F5F5", colour = "#D8D8D8"),
      strip.text = element_text(colour = col_text, face = "plain"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
           formatC(x, format = "f", digits = 3))
  )
}

sig_lab <- function(p_used) {
  dplyr::case_when(
    is.na(p_used)     ~ "ns",
    p_used < 0.001    ~ "***",
    p_used < 0.01     ~ "**",
    p_used < 0.05     ~ "*",
    TRUE              ~ "ns"
  )
}

## =========================
## 4. 读取结果表
## =========================
pattern_df <- openxlsx::read.xlsx(res_fp, sheet = "Pattern_table")

if (!"Genus" %in% colnames(pattern_df)) {
  stop("Pattern_table 中缺少 Genus 列")
}
if (!"structure_grade" %in% colnames(pattern_df)) {
  stop("Pattern_table 中缺少 structure_grade 列")
}

need_cols <- c("r_RI", "q_RI", "r_RC", "q_RC", "r_IC", "q_IC")
miss_need_cols <- setdiff(need_cols, colnames(pattern_df))
if (length(miss_need_cols) > 0) {
  stop("Pattern_table 缺少必要列：", paste(miss_need_cols, collapse = ", "))
}

pattern_df <- pattern_df %>%
  dplyr::mutate(
    Class = dplyr::case_when(
      structure_grade %in% c("三部位贯通强证据", "三部位贯通支持证据") ~ "TCG",
      structure_grade == "前段连续强证据" ~ "FCG",
      structure_grade == "前后呼应强证据" ~ "ELG",
      structure_grade == "后段连续强证据" ~ "HCG",
      TRUE ~ NA_character_
    )
  )

## =========================
## 5. 读取 genus 丰度表与 metadata，计算 CLR
## =========================
abund <- readxl::read_excel(genus_fp, sheet = 1)
meta  <- readxl::read_excel(meta_fp, sheet = 1)

req_meta_cols <- c("SampleID", "Group", "SheepID")
miss_meta_cols <- setdiff(req_meta_cols, colnames(meta))
if (length(miss_meta_cols) > 0) {
  stop("metadata 缺少必要列：", paste(miss_meta_cols, collapse = ", "))
}
if (!"Genus" %in% colnames(abund)) {
  stop("丰度表缺少 Genus 列")
}

meta_use <- meta %>%
  dplyr::filter(Group %in% c("Rum", "Ile", "Col")) %>%
  dplyr::filter(!is.na(SheepID)) %>%
  dplyr::distinct(SampleID, .keep_all = TRUE)

sample_cols <- intersect(colnames(abund), meta_use$SampleID)
if (length(sample_cols) == 0) {
  stop("丰度表与 metadata 没有重叠的样本列")
}

abund_use <- abund %>%
  dplyr::select(Genus, dplyr::all_of(sample_cols)) %>%
  dplyr::mutate(dplyr::across(-Genus, as.numeric)) %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")

mat_df <- abund_use %>%
  tibble::column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame(check.names = FALSE)

mat_df$SampleID <- rownames(mat_df)

mat_df <- mat_df %>%
  dplyr::left_join(meta_use, by = "SampleID")

paired_sheep <- mat_df %>%
  dplyr::distinct(SheepID, Group) %>%
  dplyr::count(SheepID) %>%
  dplyr::filter(n == 3) %>%
  dplyr::pull(SheepID)

mat_df <- mat_df %>%
  dplyr::filter(SheepID %in% paired_sheep) %>%
  dplyr::mutate(
    Group_show = factor(
      Group,
      levels = c("Rum", "Ile", "Col"),
      labels = c("Rumen", "Ileum", "Colon")
    )
  ) %>%
  dplyr::arrange(SheepID, Group_show)

genus_cols <- setdiff(colnames(mat_df), c("SampleID", "Group", "SheepID", "Group_show"))
comp_mat <- as.matrix(mat_df[, genus_cols, drop = FALSE])
mode(comp_mat) <- "numeric"
comp_mat[is.na(comp_mat)] <- 0

comp_repl <- zCompositions::cmultRepl(comp_mat, method = "CZM", output = "prop")
clr_mat <- compositions::clr(comp_repl)
clr_df <- as.data.frame(clr_mat, check.names = FALSE)

clr_df <- cbind(
  mat_df[, c("SampleID", "Group", "SheepID", "Group_show")],
  clr_df
)

## =========================
## 6. 图1：4个经典代表菌，每个菌一行，三边都画
## =========================
core_genus_map <- data.frame(
  Genus = c(
    "Olsenella",
    "[Ruminococcus]_gauvreauii_group",
    "Christensenellaceae_R_7_group",
    "Mycoplasma"
  ),
  Class = c("TCG", "TCG", "FCG", "HCG"),
  stringsAsFactors = FALSE
)

miss_core <- setdiff(core_genus_map$Genus, colnames(clr_df))
if (length(miss_core) > 0) {
  stop("以下核心代表菌不在 CLR 数据中：", paste(miss_core, collapse = ", "))
}

make_pair_df <- function(genus_name) {
  tmp <- clr_df %>%
    dplyr::select(SheepID, Group_show, dplyr::all_of(genus_name)) %>%
    dplyr::rename(CLR = dplyr::all_of(genus_name)) %>%
    tidyr::pivot_wider(names_from = Group_show, values_from = CLR)
  tmp
}

pair_spec <- data.frame(
  Pair = c("Rumen–Ileum", "Rumen–Colon", "Ileum–Colon"),
  xcol = c("Rumen", "Rumen", "Ileum"),
  ycol = c("Ileum", "Colon", "Colon"),
  xlab = c("Rumen CLR", "Rumen CLR", "Ileum CLR"),
  ylab = c("Ileum CLR", "Colon CLR", "Colon CLR"),
  rcol = c("r_RI", "r_RC", "r_IC"),
  qcol = c("q_RI", "q_RC", "q_IC"),
  stringsAsFactors = FALSE
)

scatter_df <- dplyr::bind_rows(lapply(seq_len(nrow(core_genus_map)), function(i) {
  genus_name <- core_genus_map$Genus[i]
  class_name <- core_genus_map$Class[i]
  wide_df <- make_pair_df(genus_name)
  
  dplyr::bind_rows(lapply(seq_len(nrow(pair_spec)), function(j) {
    dd <- wide_df[, c("SheepID", pair_spec$xcol[j], pair_spec$ycol[j])]
    colnames(dd) <- c("SheepID", "x", "y")
    dd$Genus <- genus_name
    dd$Class <- class_name
    dd$Pair  <- pair_spec$Pair[j]
    dd$xlab  <- pair_spec$xlab[j]
    dd$ylab  <- pair_spec$ylab[j]
    dd
  }))
}))

scatter_df <- scatter_df %>%
  dplyr::mutate(
    Genus = factor(
      Genus,
      levels = c(
        "Olsenella",
        "[Ruminococcus]_gauvreauii_group",
        "Christensenellaceae_R_7_group",
        "Mycoplasma"
      )
    ),
    Pair = factor(Pair, levels = c("Rumen–Ileum", "Rumen–Colon", "Ileum–Colon"))
  )

anno_scatter <- dplyr::bind_rows(lapply(seq_len(nrow(core_genus_map)), function(i) {
  genus_name <- core_genus_map$Genus[i]
  class_name <- core_genus_map$Class[i]
  rowi <- pattern_df %>% dplyr::filter(Genus == genus_name)
  
  data.frame(
    Genus = genus_name,
    Class = class_name,
    Pair  = pair_spec$Pair,
    lab   = c(
      paste0("r = ", fmt_num(rowi$r_RI[1]), "\nq = ", fmt_p(rowi$q_RI[1]), ifelse(sig_lab(rowi$q_RI[1]) == "ns", "", sig_lab(rowi$q_RI[1]))),
      paste0("r = ", fmt_num(rowi$r_RC[1]), "\nq = ", fmt_p(rowi$q_RC[1]), ifelse(sig_lab(rowi$q_RC[1]) == "ns", "", sig_lab(rowi$q_RC[1]))),
      paste0("r = ", fmt_num(rowi$r_IC[1]), "\nq = ", fmt_p(rowi$q_IC[1]), ifelse(sig_lab(rowi$q_IC[1]) == "ns", "", sig_lab(rowi$q_IC[1])))
    ),
    stringsAsFactors = FALSE
  )
}))

anno_scatter <- anno_scatter %>%
  dplyr::mutate(
    Genus = factor(
      Genus,
      levels = c(
        "Olsenella",
        "[Ruminococcus]_gauvreauii_group",
        "Christensenellaceae_R_7_group",
        "Mycoplasma"
      )
    ),
    Pair = factor(Pair, levels = c("Rumen–Ileum", "Rumen–Colon", "Ileum–Colon"))
  )

p1 <- ggplot(scatter_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = Pair),
    size = 2.0,
    alpha = 0.82
  ) +
  geom_smooth(
    aes(colour = Pair),
    method = "lm",
    se = FALSE,
    linewidth = 0.7,
    alpha = 0.9
  ) +
  facet_grid(
    rows = vars(Genus),
    cols = vars(Pair),
    scales = "free"
  ) +
  scale_colour_manual(values = pair_line_map, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_sci_soft(base_size = 11.2) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 10.2),
    strip.text.y = element_text(size = 10.2, angle = 0),
    axis.text = element_text(size = 9.3)
  ) +
  geom_text(
    data = anno_scatter,
    aes(x = -Inf, y = Inf, label = lab),
    hjust = -0.05, vjust = 1.1,
    inherit.aes = FALSE,
    size = 3.2,
    colour = col_text,
    family = "sans"
  )

pdf(
  file = file.path(out_dir, "Fig1_core_genera_all_three_pairs_scatter.pdf"),
  width = 10.8,
  height = 9.8,
  family = "sans",
  useDingbats = FALSE
)
print(p1)
dev.off()

## =========================
## 7. 图2：9个代表菌 × 3条边 相关证据矩阵图
## =========================
rep9 <- data.frame(
  Genus = c(
    "Olsenella",
    "[Ruminococcus]_gauvreauii_group",
    "Candidatus_Saccharimonas",
    "Bifidobacterium",
    "Bacteroides",
    "Christensenellaceae_R_7_group",
    "Parabacteroides",
    "Prevotellaceae_NK3B31_group",
    "Mycoplasma"
  ),
  Class = c(
    "TCG", "TCG", "TCG",
    "FCG", "FCG", "FCG",
    "ELG", "ELG", "HCG"
  ),
  stringsAsFactors = FALSE
)

miss_rep9 <- setdiff(rep9$Genus, pattern_df$Genus)
if (length(miss_rep9) > 0) {
  stop("以下代表菌不在 Pattern_table 中：", paste(miss_rep9, collapse = ", "))
}

mat_df_plot <- pattern_df %>%
  dplyr::filter(Genus %in% rep9$Genus) %>%
  dplyr::select(Genus, Class, r_RI, q_RI, r_RC, q_RC, r_IC, q_IC) %>%
  dplyr::left_join(rep9, by = "Genus", suffix = c("", "_keep")) %>%
  dplyr::mutate(Class = Class_keep) %>%
  dplyr::select(-Class_keep) %>%
  tidyr::pivot_longer(
    cols = c(r_RI, q_RI, r_RC, q_RC, r_IC, q_IC),
    names_to = c(".value", "PairCode"),
    names_pattern = "(r|q)_(RI|RC|IC)"
  ) %>%
  dplyr::mutate(
    Pair = dplyr::recode(
      PairCode,
      "RI" = "Rumen–Ileum",
      "RC" = "Rumen–Colon",
      "IC" = "Ileum–Colon"
    ),
    sig_txt = sig_lab(q),
    Pair = factor(Pair, levels = c("Rumen–Ileum", "Rumen–Colon", "Ileum–Colon"))
  )

row_order <- c(
  "Olsenella",
  "[Ruminococcus]_gauvreauii_group",
  "Candidatus_Saccharimonas",
  "Bifidobacterium",
  "Bacteroides",
  "Christensenellaceae_R_7_group",
  "Parabacteroides",
  "Prevotellaceae_NK3B31_group",
  "Mycoplasma"
)

mat_df_plot <- mat_df_plot %>%
  dplyr::mutate(
    Genus = factor(Genus, levels = rev(row_order))
  )

class_anno <- rep9 %>%
  dplyr::mutate(
    Pair = factor("Class", levels = c("Class", "Rumen–Ileum", "Rumen–Colon", "Ileum–Colon")),
    Genus = factor(Genus, levels = rev(row_order)),
    fill_col = dplyr::recode(
      Class,
      "TCG" = col_tcg_fill,
      "FCG" = col_fcg_fill,
      "ELG" = col_elg_fill,
      "HCG" = col_hcg_fill
    )
  )

plot_matrix <- mat_df_plot %>%
  dplyr::mutate(
    Pair = factor(Pair, levels = c("Class", "Rumen–Ileum", "Rumen–Colon", "Ileum–Colon"))
  )

p2 <- ggplot() +
  geom_tile(
    data = class_anno,
    aes(x = Pair, y = Genus),
    fill = class_anno$fill_col,
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    data = class_anno,
    aes(x = Pair, y = Genus, label = Class),
    colour = col_text,
    size = 3.15,
    family = "sans"
  ) +
  geom_tile(
    data = plot_matrix,
    aes(x = Pair, y = Genus, fill = r),
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    data = plot_matrix,
    aes(x = Pair, y = Genus, label = sig_txt),
    colour = col_text,
    size = 3.2,
    family = "sans"
  ) +
  scale_fill_gradient2(
    low = "#8FB3D9",
    mid = "white",
    high = "#D8A47F",
    midpoint = 0,
    na.value = col_na,
    limits = c(-1, 1)
  ) +
  labs(x = NULL, y = NULL, fill = "Spearman r") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(colour = col_text, angle = 0, hjust = 0.5),
    axis.text.y = element_text(colour = col_text),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_text(colour = col_text),
    legend.text = element_text(colour = col_text),
    plot.margin = margin(10, 14, 10, 10)
  )

pdf(
  file = file.path(out_dir, "Fig2_representative_9genera_correlation_matrix_signif_only.pdf"),
  width = 8.2,
  height = 6.6,
  family = "sans",
  useDingbats = FALSE
)
print(p2)
dev.off()

message("完成：已输出两个 PDF 至 ", out_dir)