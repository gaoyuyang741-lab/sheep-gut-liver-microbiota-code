# =========================================================

# =========================
# 0. 安装/加载包
# =========================
need_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "stringr",
  "readxl", "grid", "data.table"
)
new_pkgs <- need_pkgs[!(need_pkgs %in% installed.packages()[, "Package"])]
if(length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(grid)
library(data.table)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate

# =========================
# 1. 文件路径
# =========================
# Input and output directories
# Please place the raw alpha-diversity tables in "data/diversity/alpha/input".
# Pairwise test results should be placed in "results/diversity/alpha".
# Figure files will be saved in "figures/diversity/alpha".

alpha_fp_ri <- file.path(
  "data", "diversity", "alpha", "input",
  "RI_alpha_diversity.xls"
)

alpha_fp_c <- file.path(
  "data", "diversity", "alpha", "input",
  "C_alpha_diversity.xls"
)

pairwise_fp <- file.path(
  "results", "diversity", "alpha",
  "alpha_pairwise_wilcox.xlsx"
)

out_dir <- file.path("figures", "diversity", "alpha")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# =========================
# 2. 全局作图参数
# =========================
base_family <- "Helvetica"

col_fill <- c(
  "Rum"   = "#BFD6D1",
  "Ile"   = "#E8C9C4",
  "Colon" = "#D8CFE6"
)

col_line <- c(
  "Rum"   = "#6F9F98",
  "Ile"   = "#C99690",
  "Colon" = "#8E7EAF"
)

pair_line_color <- "grey85"
pair_line_alpha <- 0.90
pair_line_width <- 0.28

point_size         <- 2.3
point_jitter_width <- 0.015
point_stroke       <- 0.60

box_width      <- 0.42
box_line_width <- 0.90

bracket_line_width <- 0.90
star_text_size     <- 5.2

x_title_size <- 14
y_text_size  <- 15

legend_text_size <- 11
legend_key_size  <- 0.65

plot_width  <- 4.8
plot_height <- 4.4
plot_dpi    <- 600

# =========================
# 3. 读取 alpha 文件函数
# =========================
read_alpha_file <- function(in_file) {
  df <- tryCatch(
    read.delim(
      in_file,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(df)) {
    df <- tryCatch(
      data.table::fread(
        in_file,
        data.table = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  }
  
  if (is.null(df)) {
    stop(paste0("无法读取文件：", in_file))
  }
  
  names(df) <- str_replace_all(names(df), "\\s+", "_")
  names(df) <- str_replace_all(names(df), "-", "_")
  
  if (!"Sample_ID" %in% names(df)) {
    if ("SampleID" %in% names(df)) {
      df <- df %>% rename(Sample_ID = SampleID)
    } else if ("Sample" %in% names(df)) {
      df <- df %>% rename(Sample_ID = Sample)
    } else {
      names(df)[1] <- "Sample_ID"
    }
  }
  
  pd_candidates <- c("PD_whole_tree", "PD_whole_", "PD_whole_tree_index", "PD_whole")
  pd_col <- intersect(pd_candidates, names(df))
  if (length(pd_col) == 0) {
    pd_col2 <- grep("^PD", names(df), value = TRUE)
    if (length(pd_col2) == 0) {
      stop(paste0("找不到 PD 列。当前列名：", paste(names(df), collapse = ", ")))
    } else {
      pd_col <- pd_col2[1]
    }
  } else {
    pd_col <- pd_col[1]
  }
  
  df <- df %>% rename(PD_whole_tree = all_of(pd_col))
  
  num_cols <- c("Feature", "ACE", "Chao1", "Simpson", "Shannon", "PD_whole_tree", "Coverage")
  df <- df %>%
    mutate(across(any_of(num_cols), ~ suppressWarnings(as.numeric(.))))
  
  return(df)
}

# =========================
# 4. 读取原始数据
# =========================
alpha_ri <- read_alpha_file(alpha_fp_ri)
alpha_c  <- read_alpha_file(alpha_fp_c)
alpha    <- bind_rows(alpha_ri, alpha_c)

# =========================
# 5. 整理 alpha 原始数据
# =========================
alpha2 <- alpha %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Site = case_when(
      str_detect(Sample_ID, "^R-") ~ "Rum",
      str_detect(Sample_ID, "^I-") ~ "Ile",
      str_detect(Sample_ID, "^C-") ~ "Colon",
      TRUE ~ NA_character_
    ),
    pair_id = str_replace(Sample_ID, "^[RIC]-", "")
  )

triplet_check <- alpha2 %>%
  filter(!is.na(Site)) %>%
  count(pair_id, Site) %>%
  pivot_wider(names_from = Site, values_from = n, values_fill = 0)

for (s in c("Rum", "Ile", "Colon")) {
  if (!s %in% names(triplet_check)) triplet_check[[s]] <- 0
}

valid_triplets <- triplet_check %>%
  filter(Rum == 1, Ile == 1, Colon == 1) %>%
  pull(pair_id)

alpha2 <- alpha2 %>%
  filter(pair_id %in% valid_triplets) %>%
  mutate(
    Site = factor(Site, levels = c("Rum", "Ile", "Colon"))
  )

# =========================
# 6. 读取并整理两两比较结果
# =========================
pairwise_stat <- read_excel(pairwise_fp)

if (!"index" %in% colnames(pairwise_stat)) {
  stop("alpha_pairwise_wilcox.xlsx 中未找到列名：'index'")
}
if (!"contrast" %in% colnames(pairwise_stat)) {
  stop("alpha_pairwise_wilcox.xlsx 中未找到列名：'contrast'")
}

pairwise_stat2 <- pairwise_stat %>%
  rename(metric = index) %>%
  mutate(
    p_used = case_when(
      "p_fdr_within_index" %in% colnames(.) ~ p_fdr_within_index,
      "p_fdr_global" %in% colnames(.) ~ p_fdr_global,
      "p_raw" %in% colnames(.) ~ p_raw,
      TRUE ~ NA_real_
    ),
    sig_text = case_when(
      is.na(p_used)   ~ "ns",
      p_used < 0.001  ~ "***",
      p_used < 0.01   ~ "**",
      p_used < 0.05   ~ "*",
      TRUE            ~ "ns"
    ),
    x1 = case_when(
      contrast == "Rum vs Ile"   ~ 1,
      contrast == "Rum vs Colon" ~ 1,
      contrast == "Ile vs Colon" ~ 2,
      TRUE ~ NA_real_
    ),
    x2 = case_when(
      contrast == "Rum vs Ile"   ~ 2,
      contrast == "Rum vs Colon" ~ 3,
      contrast == "Ile vs Colon" ~ 3,
      TRUE ~ NA_real_
    )
  )

# =========================
# 7. 指标名映射
# =========================
xlab_map <- c(
  "ACE"           = "ACE index",
  "Chao1"         = "Chao1 index",
  "Shannon"       = "Shannon index",
  "Simpson"       = "Simpson index",
  "PD_whole_tree" = "PD whole tree"
)

# =========================
# 8. 单指标作图函数
# =========================
plot_alpha_triplet_pairwise <- function(dat, stat_tab, metric_name,
                                        width = plot_width,
                                        height = plot_height,
                                        dpi = plot_dpi) {
  
  if(!metric_name %in% colnames(dat)){
    stop(paste0("原始 alpha 数据中不存在指标列：", metric_name))
  }
  
  plot_df <- dat %>%
    dplyr::select(pair_id, Site, dplyr::all_of(metric_name)) %>%
    dplyr::rename(value = dplyr::all_of(metric_name)) %>%
    dplyr::filter(!is.na(value), !is.na(Site), !is.na(pair_id)) %>%
    dplyr::mutate(
      Site = factor(Site, levels = c("Rum", "Ile", "Colon"))
    )
  
  stat_sub <- stat_tab %>%
    filter(metric == metric_name) %>%
    arrange(x2 - x1, x1, x2)
  
  xlab_name <- if(metric_name %in% names(xlab_map)) xlab_map[[metric_name]] else metric_name
  
  y_min <- min(plot_df$value, na.rm = TRUE)
  y_max <- max(plot_df$value, na.rm = TRUE)
  y_range <- y_max - y_min
  if(y_range == 0) y_range <- abs(y_max) * 0.1 + 1
  
  # 三层括号高度
  bracket_base <- y_max + y_range * 0.08
  bracket_step <- y_range * 0.10
  upper_lim    <- y_max + y_range * 0.42
  
  p <- ggplot(plot_df, aes(x = Site, y = value)) +
    geom_line(
      aes(group = pair_id),
      color = pair_line_color,
      alpha = pair_line_alpha,
      linewidth = pair_line_width
    ) +
    geom_boxplot(
      aes(fill = Site, color = Site),
      width = box_width,
      outlier.shape = NA,
      linewidth = box_line_width
    ) +
    geom_point(
      aes(fill = Site, color = Site),
      size = point_size,
      shape = 21,
      stroke = point_stroke,
      position = position_jitter(width = point_jitter_width, height = 0)
    ) +
    scale_fill_manual(
      values = col_fill,
      breaks = c("Rum", "Ile", "Colon"),
      labels = c("Rum", "Ile", "Colon"),
      name = NULL
    ) +
    scale_color_manual(
      values = col_line,
      breaks = c("Rum", "Ile", "Colon"),
      labels = c("Rum", "Ile", "Colon"),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(y_min, upper_lim),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = xlab_name, y = NULL, title = NULL) +
    theme_classic(base_size = 14, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      
      axis.title.x = element_text(
        size = x_title_size,
        color = "black",
        family = base_family,
        margin = margin(t = 8)
      ),
      axis.title.y = element_blank(),
      
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(
        size = y_text_size,
        color = "black",
        family = base_family
      ),
      
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 1.0, color = "black"),
      axis.line    = element_line(linewidth = 1.0, color = "black"),
      axis.ticks.length = unit(0.20, "cm"),
      
      legend.position = c(1.06, 0.88),
      legend.justification = c(1, 1),
      legend.text = element_text(
        size = legend_text_size,
        color = "black",
        family = base_family
      ),
      legend.key = element_blank(),
      legend.key.size = unit(legend_key_size, "cm"),
      legend.background = element_blank(),
      
      panel.grid = element_blank(),
      plot.margin = margin(t = 16, r = 14, b = 10, l = 12)
    )
  
  # 叠加三条两两比较括号
  if (nrow(stat_sub) > 0) {
    for (i in seq_len(nrow(stat_sub))) {
      x1 <- stat_sub$x1[i]
      x2 <- stat_sub$x2[i]
      lab <- stat_sub$sig_text[i]
      
      y_now <- bracket_base + (i - 1) * bracket_step
      
      p <- p +
        annotate(
          "segment",
          x = x1, xend = x2,
          y = y_now, yend = y_now,
          linewidth = bracket_line_width, color = "black"
        ) +
        annotate(
          "segment",
          x = x1, xend = x1,
          y = y_now - y_range * 0.018,
          yend = y_now + y_range * 0.018,
          linewidth = bracket_line_width, color = "black"
        ) +
        annotate(
          "segment",
          x = x2, xend = x2,
          y = y_now - y_range * 0.018,
          yend = y_now + y_range * 0.018,
          linewidth = bracket_line_width, color = "black"
        ) +
        annotate(
          "text",
          x = (x1 + x2) / 2,
          y = y_now + y_range * 0.03,
          label = lab,
          size = star_text_size,
          family = base_family
        )
    }
  }
  
  png_fp <- file.path(out_dir, paste0(metric_name, "_triplet_pairwise.png"))
  ggsave(
    filename = png_fp,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  pdf_fp <- file.path(out_dir, paste0(metric_name, "_triplet_pairwise.pdf"))
  grDevices::pdf(
    file = pdf_fp,
    width = width,
    height = height,
    family = base_family,
    useDingbats = FALSE,
    version = "1.4",
    paper = "special",
    compress = FALSE
  )
  print(p)
  dev.off()
  
  return(p)
}

# =========================
# 9. 批量输出
# =========================
metrics_to_plot <- c("ACE", "Chao1", "Shannon", "Simpson", "PD_whole_tree")

for(met in metrics_to_plot){
  plot_alpha_triplet_pairwise(
    dat = alpha2,
    stat_tab = pairwise_stat2,
    metric_name = met
  )
}