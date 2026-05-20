# =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

# -----------------------------
# 0. 加载程序包
# -----------------------------
pkgs <- c("readxl", "dplyr", "ggplot2", "grid")
need <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(need) > 0) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------------
# 1. 文件路径
# -----------------------------
# Input and output directories
# This script reads beta-diversity PCoA coordinate files from "results/diversity/beta/03_PCoA_coordinates".
# Figure files will be saved in "figures/diversity/beta".

base_dir  <- file.path("results", "diversity", "beta")
coord_dir <- file.path(base_dir, "03_PCoA_coordinates")
outdir    <- file.path("figures", "diversity", "beta")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bray_fp <- file.path(coord_dir, "bray_curtis_PCoA坐标.xlsx")
jacc_fp <- file.path(coord_dir, "binary_jaccard_PCoA坐标.xlsx")

# -----------------------------
# 2. 读取坐标文件
# -----------------------------
bray_plot_df <- readxl::read_excel(bray_fp)
jacc_plot_df <- readxl::read_excel(jacc_fp)

# -----------------------------
# 3. 数据整理
# -----------------------------
clean_pcoa_df <- function(df){
  req_cols <- c("SampleID", "PCoA1", "PCoA2", "Site")
  miss_cols <- setdiff(req_cols, colnames(df))
  if(length(miss_cols) > 0){
    stop("坐标文件缺少以下列：", paste(miss_cols, collapse = ", "))
  }
  
  df <- df %>%
    dplyr::mutate(
      SampleID = as.character(SampleID),
      Site = as.character(Site)
    )
  
  df$Site <- dplyr::case_when(
    df$Site %in% c("Rum", "Rumen", "rum", "rumen") ~ "Rum",
    df$Site %in% c("Ile", "Ileum", "ile", "ileum") ~ "Ile",
    df$Site %in% c("Col", "Colon", "col", "colon") ~ "Col",
    TRUE ~ df$Site
  )
  
  # 图例顺序
  df$Site <- factor(df$Site, levels = c("Rum", "Ile", "Col"))
  
  return(df)
}

bray_plot_df <- clean_pcoa_df(bray_plot_df)
jacc_plot_df <- clean_pcoa_df(jacc_plot_df)

# -----------------------------
# 4. 解释度
# 这里请按你的真实结果改数字
# -----------------------------
bray_pc1_var <- 11.17
bray_pc2_var <- 9.64

jacc_pc1_var <- 3.11
jacc_pc2_var <- 2.78

# -----------------------------
# 5. 重新配色
# 新方案：低饱和、柔和、期刊风
# -----------------------------
fill_cols <- c(
  "Ile" = "#C7D4E2",   # 灰蓝填充
  "Rum" = "#D7CADF",   # 灰紫填充
  "Col" = "#D7E1D6"    # 灰绿填充
)

point_cols <- c(
  "Ile" = "#7C9BB8",   # 灰蓝点/线
  "Rum" = "#8A5A9E",   # 灰紫点/线
  "Col" = "#7F9A7A"    # 灰绿点/线
)

# -----------------------------
# 6. 主题函数
# -----------------------------
pcoa_theme_ref <- function(){
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_blank(),
      
      axis.title = element_text(size = 14, colour = "black"),
      axis.text  = element_text(size = 11.5, colour = "black"),
      
      axis.line  = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = grid::unit(0.14, "cm"),
      
      panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 0.6),
      
      legend.title = element_blank(),
      legend.text  = element_text(size = 11),
      legend.position = c(1.02, 1.03),
      legend.justification = c(1, 1),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

# -----------------------------
# 7. 作图函数
# -----------------------------
make_pcoa_plot_ref_3site <- function(df, xvar, yvar, xlab_txt, ylab_txt){
  
  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    
    stat_ellipse(
      aes(fill = Site, colour = Site),
      geom = "polygon",
      type = "norm",
      level = 0.95,
      alpha = 0.38,
      linewidth = 0.85,
      show.legend = FALSE
    ) +
    
    geom_point(
      aes(colour = Site),
      size = 2.6,
      alpha = 0.95
    ) +
    
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_colour_manual(values = point_cols, drop = FALSE) +
    
    labs(
      x = xlab_txt,
      y = ylab_txt
    ) +
    
    coord_fixed() +
    pcoa_theme_ref()
  
  return(p)
}

# -----------------------------
# 8. 生成图形
# -----------------------------
p_bray <- make_pcoa_plot_ref_3site(
  df = bray_plot_df,
  xvar = "PCoA1",
  yvar = "PCoA2",
  xlab_txt = paste0("PCoA1 (", bray_pc1_var, "%)"),
  ylab_txt = paste0("PCoA2 (", bray_pc2_var, "%)")
)

p_jacc <- make_pcoa_plot_ref_3site(
  df = jacc_plot_df,
  xvar = "PCoA1",
  yvar = "PCoA2",
  xlab_txt = paste0("PCoA1 (", jacc_pc1_var, "%)"),
  ylab_txt = paste0("PCoA2 (", jacc_pc2_var, "%)")
)

# -----------------------------
# 9. 显示
# -----------------------------
print(p_bray)
print(p_jacc)

# -----------------------------
# 10. 保存
# -----------------------------
ggsave(
  filename = file.path(outdir, "PCoA_BrayCurtis_PC1_PC2_refstyle_3site_v2.tiff"),
  plot = p_bray,
  width = 5.2, height = 4.4, dpi = 600, compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "PCoA_BrayCurtis_PC1_PC2_refstyle_3site_v2.pdf"),
  plot = p_bray,
  width = 5.2, height = 4.4
)

ggsave(
  filename = file.path(outdir, "PCoA_BrayCurtis_PC1_PC2_refstyle_3site_v2.png"),
  plot = p_bray,
  width = 5.2, height = 4.4, dpi = 600
)

ggsave(
  filename = file.path(outdir, "PCoA_BinaryJaccard_PC1_PC2_refstyle_3site_v2.tiff"),
  plot = p_jacc,
  width = 5.2, height = 4.4, dpi = 600, compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "PCoA_BinaryJaccard_PC1_PC2_refstyle_3site_v2.pdf"),
  plot = p_jacc,
  width = 5.2, height = 4.4
)

ggsave(
  filename = file.path(outdir, "PCoA_BinaryJaccard_PC1_PC2_refstyle_3site_v2.png"),
  plot = p_jacc,
  width = 5.2, height = 4.4, dpi = 600
)

# -----------------------------
# 11. 保存作图数据
# -----------------------------
write.csv(
  bray_plot_df,
  file.path(outdir, "BrayCurtis_PCoA_plot_data_3site_v2.csv"),
  row.names = FALSE
)

write.csv(
  jacc_plot_df,
  file.path(outdir, "BinaryJaccard_PCoA_plot_data_3site_v2.csv"),
  row.names = FALSE
)

cat("Done! Files saved in:\n", outdir, "\n")