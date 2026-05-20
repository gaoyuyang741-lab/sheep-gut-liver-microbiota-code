# =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

# -------------------------------
# 0. 加载程序包
# -------------------------------
pkg_needed <- c(
  "readxl", "openxlsx", "dplyr", "stringr", "vegan",
  "ggplot2", "ape", "tibble", "purrr"
)
pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, dependencies = TRUE)
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

# -------------------------------
# 1. 路径设置
# -------------------------------
# 1. Path settings
# Please place the required beta-diversity input files in "data/diversity/beta/input".
# Output files will be saved in "results/diversity/beta".

in_dir  <- file.path("data", "diversity", "beta", "input")
out_dir <- file.path("results", "diversity", "beta")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

abund_file <- file.path(in_dir, "all_abundance.xlsx")
meta_file  <- file.path(in_dir, "metadata.xlsx")

# -------------------------------
# 2. 输出文件夹设置
# -------------------------------
dir.create(file.path(out_dir, "01_整理后的输入表"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "02_距离矩阵"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "03_PCoA坐标"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "04_统计结果"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "05_图形"), showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 3. 读取总丰度表
# -------------------------------
abund_raw <- read.xlsx(abund_file)

tax_cols <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
miss_tax <- setdiff(tax_cols, colnames(abund_raw))
if (length(miss_tax) > 0) {
  stop("总丰度表缺少以下分类列：", paste(miss_tax, collapse = ", "))
}

sample_cols <- setdiff(colnames(abund_raw), tax_cols)
if (length(sample_cols) == 0) {
  stop("总丰度表中没有检测到样本列。")
}

# 分类列清理
for (cc in tax_cols) {
  abund_raw[[cc]] <- trimws(as.character(abund_raw[[cc]]))
  abund_raw[[cc]][is.na(abund_raw[[cc]])] <- ""
}

# 样本列转数值
for (cc in sample_cols) {
  abund_raw[[cc]] <- suppressWarnings(as.numeric(abund_raw[[cc]]))
  abund_raw[[cc]][is.na(abund_raw[[cc]])] <- 0
}

# 构建特征键
abund_raw$FeatureID <- apply(abund_raw[, tax_cols, drop = FALSE], 1, function(x) {
  paste(x, collapse = "|||")
})

# 检查 FeatureID 是否重复
dup_feature <- abund_raw %>%
  dplyr::count(FeatureID, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(dup_feature) > 0) {
  message("检测到重复 FeatureID，已自动按 FeatureID 合并并对样本列求和。")
  
  abund_raw <- abund_raw %>%
    dplyr::group_by(FeatureID, across(all_of(tax_cols))) %>%
    dplyr::summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# -------------------------------
# 4. 读取 metadata
# -------------------------------
meta <- read.xlsx(meta_file)

need_meta_cols <- c("SampleID", "Group", "SheepID")
miss_meta <- setdiff(need_meta_cols, colnames(meta))
if (length(miss_meta) > 0) {
  stop("metadata 缺少以下列：", paste(miss_meta, collapse = ", "))
}

meta <- meta %>%
  dplyr::select(SampleID, Group, SheepID) %>%
  dplyr::mutate(
    SampleID = trimws(as.character(SampleID)),
    Group    = trimws(as.character(Group)),
    SheepID  = trimws(as.character(SheepID))
  )

# 把 Group 重命名为 Site
meta <- meta %>%
  dplyr::rename(Site = Group)

# Site 顺序固定
meta$Site <- factor(meta$Site, levels = c("Rum", "Ile", "Col"))

# 从 SampleID 中提取真正处理组
# 例如 R-BA1-1 / I-CM-3 / C-CON-6 -> BA1 / CM / CON
meta$Treat <- stringr::str_split_fixed(meta$SampleID, "-", 3)[, 2]
meta$Treat <- factor(meta$Treat, levels = c("BA1", "BA2", "CM", "CMPS", "CON"))

# -------------------------------
# 5. 核对样本一致性并重排
# -------------------------------
sample_in_abund <- sample_cols
sample_in_meta  <- meta$SampleID

miss_in_meta  <- setdiff(sample_in_abund, sample_in_meta)
miss_in_abund <- setdiff(sample_in_meta, sample_in_abund)

if (length(miss_in_meta) > 0) {
  stop("总丰度表中以下样本在 metadata 中不存在：\n", paste(miss_in_meta, collapse = ", "))
}
if (length(miss_in_abund) > 0) {
  stop("metadata 中以下样本在总丰度表中不存在：\n", paste(miss_in_abund, collapse = ", "))
}

# 按 metadata 顺序重排丰度表样本列
abund <- abund_raw[, c(tax_cols, "FeatureID", meta$SampleID), drop = FALSE]

# 转为 样本 × 特征 矩阵
otu_mat <- t(as.matrix(abund[, meta$SampleID, drop = FALSE]))
colnames(otu_mat) <- abund$FeatureID
rownames(otu_mat) <- meta$SampleID

# 检查是否为非负数
if (any(otu_mat < 0, na.rm = TRUE)) {
  stop("丰度矩阵中存在负值，请检查输入表。")
}

# 去掉全 0 特征
feat_sum <- colSums(otu_mat, na.rm = TRUE)
otu_mat <- otu_mat[, feat_sum > 0, drop = FALSE]

# 去掉全 0 样本
sample_sum <- rowSums(otu_mat, na.rm = TRUE)
if (any(sample_sum == 0)) {
  zero_samples <- names(sample_sum)[sample_sum == 0]
  stop("以下样本总丰度为 0，无法进行 β 多样性分析：\n", paste(zero_samples, collapse = ", "))
}

# 导出整理后的输入表
input_abund_export <- data.frame(
  SampleID = rownames(otu_mat),
  otu_mat,
  check.names = FALSE
)
write.xlsx(input_abund_export,
           file.path(out_dir, "01_整理后的输入表", "01_样本x特征丰度矩阵.xlsx"),
           rowNames = FALSE)

write.xlsx(meta,
           file.path(out_dir, "01_整理后的输入表", "02_整理后的metadata.xlsx"),
           rowNames = FALSE)

# -------------------------------
# 6. 定义函数
# -------------------------------

# 6.1 计算 PCoA 坐标
get_pcoa_df <- function(dist_obj, meta_df) {
  pcoa_res <- ape::pcoa(dist_obj)
  
  coord <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
  colnames(coord) <- c("PCoA1", "PCoA2")
  coord$SampleID <- rownames(coord)
  
  var_exp <- round(pcoa_res$values$Relative_eig[1:2] * 100, 2)
  
  coord <- coord %>%
    dplyr::left_join(meta_df, by = "SampleID")
  
  return(list(coord = coord, var_exp = var_exp, pcoa_res = pcoa_res))
}

# 6.2 画 PCoA 图
plot_pcoa <- function(coord_df, var_exp, title_text, out_pdf, out_png) {
  p <- ggplot(coord_df, aes(x = PCoA1, y = PCoA2, color = Site, fill = Site)) +
    stat_ellipse(geom = "polygon", alpha = 0.18, level = 0.95, linewidth = 0.4) +
    geom_point(size = 3, alpha = 0.9) +
    theme_bw(base_size = 14) +
    labs(
      x = paste0("PCoA1 (", var_exp[1], "%)"),
      y = paste0("PCoA2 (", var_exp[2], "%)")
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  ggsave(out_pdf, p, width = 7.2, height = 5.8, units = "in")
  ggsave(out_png, p, width = 7.2, height = 5.8, units = "in", dpi = 600)
  
  return(p)
}

# 6.3 PERMANOVA
run_permanova <- function(dist_obj, meta_df) {
  # 主模型：部位效应
  adonis_site <- vegan::adonis2(dist_obj ~ Site, data = meta_df, permutations = 9999)
  
  # 控制处理组后的部位效应
  adonis_site_treat <- vegan::adonis2(dist_obj ~ Treat + Site, data = meta_df, permutations = 9999)
  
  # 受限置换：按 SheepID 分层
  adonis_site_strata <- vegan::adonis2(
    dist_obj ~ Site,
    data = meta_df,
    permutations = 9999,
    strata = meta_df$SheepID
  )
  
  return(list(
    adonis_site = as.data.frame(adonis_site),
    adonis_site_treat = as.data.frame(adonis_site_treat),
    adonis_site_strata = as.data.frame(adonis_site_strata)
  ))
}

# 6.4 PERMDISP
run_permdisp <- function(dist_obj, meta_df) {
  bd <- vegan::betadisper(dist_obj, group = meta_df$Site)
  bd_perm <- vegan::permutest(bd, permutations = 9999)
  bd_anova <- anova(bd)
  
  # 各样本到组中心距离
  dist_to_centroid <- data.frame(
    SampleID = names(bd$distances),
    DistanceToCentroid = as.numeric(bd$distances),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(meta_df, by = "SampleID")
  
  # 组均值
  centroid_group_mean <- dist_to_centroid %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      N = dplyr::n(),
      Mean = mean(DistanceToCentroid, na.rm = TRUE),
      SD = sd(DistanceToCentroid, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    bd = bd,
    bd_perm_tab = as.data.frame(bd_perm$tab),
    bd_anova_tab = as.data.frame(bd_anova),
    dist_to_centroid = dist_to_centroid,
    centroid_group_mean = centroid_group_mean
  ))
}

# 6.5 两两比较
run_pairwise_adonis <- function(dist_obj, meta_df) {
  pairs <- list(
    c("Rum", "Ile"),
    c("Rum", "Col"),
    c("Ile", "Col")
  )
  
  res_list <- list()
  
  for (i in seq_along(pairs)) {
    pair_now <- pairs[[i]]
    sub_meta <- meta_df %>% dplyr::filter(Site %in% pair_now)
    
    sub_dist <- as.matrix(dist_obj)[sub_meta$SampleID, sub_meta$SampleID, drop = FALSE]
    sub_dist <- as.dist(sub_dist)
    
    adonis_pair <- vegan::adonis2(
      sub_dist ~ Site,
      data = sub_meta,
      permutations = 9999,
      strata = sub_meta$SheepID
    )
    
    bd_pair <- vegan::betadisper(sub_dist, group = sub_meta$Site)
    bd_pair_perm <- vegan::permutest(bd_pair, permutations = 9999)
    
    tmp1 <- as.data.frame(adonis_pair)
    tmp1$Comparison <- paste(pair_now, collapse = "_vs_")
    tmp1$Method <- "PERMANOVA"
    
    tmp2 <- as.data.frame(bd_pair_perm$tab)
    tmp2$Comparison <- paste(pair_now, collapse = "_vs_")
    tmp2$Method <- "PERMDISP"
    
    res_list[[paste0("adonis_", i)]] <- tmp1
    res_list[[paste0("disp_", i)]]   <- tmp2
  }
  
  pairwise_adonis <- dplyr::bind_rows(res_list[grepl("^adonis_", names(res_list))])
  pairwise_disp   <- dplyr::bind_rows(res_list[grepl("^disp_", names(res_list))])
  
  return(list(
    pairwise_adonis = pairwise_adonis,
    pairwise_disp = pairwise_disp
  ))
}

# 6.6 保存距离矩阵
save_dist_matrix <- function(dist_obj, fp) {
  dm <- as.matrix(dist_obj)
  write.table(
    dm,
    file = fp,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

# -------------------------------
# 7. 计算两种距离
# -------------------------------

# bray_curtis
dist_bray <- vegan::vegdist(otu_mat, method = "bray")

# binary_jaccard
otu_bin <- otu_mat
otu_bin[otu_bin > 0] <- 1
dist_jaccard <- vegan::vegdist(otu_bin, method = "jaccard", binary = TRUE)

# 保存距离矩阵
save_dist_matrix(dist_bray,
                 file.path(out_dir, "02_距离矩阵", "allsample.bray_curtis_dm.txt"))
save_dist_matrix(dist_jaccard,
                 file.path(out_dir, "02_距离矩阵", "allsample.binary_jaccard_dm.txt"))

# -------------------------------
# 8. Bray-Curtis 分析
# -------------------------------
bray_pcoa <- get_pcoa_df(dist_bray, meta)
write.xlsx(bray_pcoa$coord,
           file.path(out_dir, "03_PCoA坐标", "bray_curtis_PCoA坐标.xlsx"),
           rowNames = FALSE)

plot_pcoa(
  coord_df = bray_pcoa$coord,
  var_exp = bray_pcoa$var_exp,
  title_text = "Bray-Curtis",
  out_pdf = file.path(out_dir, "05_图形", "PCoA_BrayCurtis_3site.pdf"),
  out_png = file.path(out_dir, "05_图形", "PCoA_BrayCurtis_3site.png")
)

bray_permanova <- run_permanova(dist_bray, meta)
bray_permdisp  <- run_permdisp(dist_bray, meta)
bray_pairwise  <- run_pairwise_adonis(dist_bray, meta)

wb_bray <- createWorkbook()

addWorksheet(wb_bray, "1_PERMANOVA_Site")
writeData(wb_bray, "1_PERMANOVA_Site", bray_permanova$adonis_site, rowNames = TRUE)

addWorksheet(wb_bray, "2_PERMANOVA_Treat_Site")
writeData(wb_bray, "2_PERMANOVA_Treat_Site", bray_permanova$adonis_site_treat, rowNames = TRUE)

addWorksheet(wb_bray, "3_PERMANOVA_Site_strata")
writeData(wb_bray, "3_PERMANOVA_Site_strata", bray_permanova$adonis_site_strata, rowNames = TRUE)

addWorksheet(wb_bray, "4_PERMDISP_ANOVA")
writeData(wb_bray, "4_PERMDISP_ANOVA", bray_permdisp$bd_anova_tab, rowNames = TRUE)

addWorksheet(wb_bray, "5_PERMDISP_permutest")
writeData(wb_bray, "5_PERMDISP_permutest", bray_permdisp$bd_perm_tab, rowNames = TRUE)

addWorksheet(wb_bray, "6_到组中心距离")
writeData(wb_bray, "6_到组中心距离", bray_permdisp$dist_to_centroid, rowNames = FALSE)

addWorksheet(wb_bray, "7_组中心距离均值")
writeData(wb_bray, "7_组中心距离均值", bray_permdisp$centroid_group_mean, rowNames = FALSE)

addWorksheet(wb_bray, "8_pairwise_PERMANOVA")
writeData(wb_bray, "8_pairwise_PERMANOVA", bray_pairwise$pairwise_adonis, rowNames = TRUE)

addWorksheet(wb_bray, "9_pairwise_PERMDISP")
writeData(wb_bray, "9_pairwise_PERMDISP", bray_pairwise$pairwise_disp, rowNames = TRUE)

saveWorkbook(wb_bray,
             file.path(out_dir, "04_统计结果", "bray_curtis_统计结果.xlsx"),
             overwrite = TRUE)

# -------------------------------
# 9. Binary Jaccard 分析
# -------------------------------
jaccard_pcoa <- get_pcoa_df(dist_jaccard, meta)
write.xlsx(jaccard_pcoa$coord,
           file.path(out_dir, "03_PCoA坐标", "binary_jaccard_PCoA坐标.xlsx"),
           rowNames = FALSE)

plot_pcoa(
  coord_df = jaccard_pcoa$coord,
  var_exp = jaccard_pcoa$var_exp,
  title_text = "Binary Jaccard",
  out_pdf = file.path(out_dir, "05_图形", "PCoA_BinaryJaccard_3site.pdf"),
  out_png = file.path(out_dir, "05_图形", "PCoA_BinaryJaccard_3site.png")
)

jaccard_permanova <- run_permanova(dist_jaccard, meta)
jaccard_permdisp  <- run_permdisp(dist_jaccard, meta)
jaccard_pairwise  <- run_pairwise_adonis(dist_jaccard, meta)

wb_jaccard <- createWorkbook()

addWorksheet(wb_jaccard, "1_PERMANOVA_Site")
writeData(wb_jaccard, "1_PERMANOVA_Site", jaccard_permanova$adonis_site, rowNames = TRUE)

addWorksheet(wb_jaccard, "2_PERMANOVA_Treat_Site")
writeData(wb_jaccard, "2_PERMANOVA_Treat_Site", jaccard_permanova$adonis_site_treat, rowNames = TRUE)

addWorksheet(wb_jaccard, "3_PERMANOVA_Site_strata")
writeData(wb_jaccard, "3_PERMANOVA_Site_strata", jaccard_permanova$adonis_site_strata, rowNames = TRUE)

addWorksheet(wb_jaccard, "4_PERMDISP_ANOVA")
writeData(wb_jaccard, "4_PERMDISP_ANOVA", jaccard_permdisp$bd_anova_tab, rowNames = TRUE)

addWorksheet(wb_jaccard, "5_PERMDISP_permutest")
writeData(wb_jaccard, "5_PERMDISP_permutest", jaccard_permdisp$bd_perm_tab, rowNames = TRUE)

addWorksheet(wb_jaccard, "6_到组中心距离")
writeData(wb_jaccard, "6_到组中心距离", jaccard_permdisp$dist_to_centroid, rowNames = FALSE)

addWorksheet(wb_jaccard, "7_组中心距离均值")
writeData(wb_jaccard, "7_组中心距离均值", jaccard_permdisp$centroid_group_mean, rowNames = FALSE)

addWorksheet(wb_jaccard, "8_pairwise_PERMANOVA")
writeData(wb_jaccard, "8_pairwise_PERMANOVA", jaccard_pairwise$pairwise_adonis, rowNames = TRUE)

addWorksheet(wb_jaccard, "9_pairwise_PERMDISP")
writeData(wb_jaccard, "9_pairwise_PERMDISP", jaccard_pairwise$pairwise_disp, rowNames = TRUE)

saveWorkbook(wb_jaccard,
             file.path(out_dir, "04_统计结果", "binary_jaccard_统计结果.xlsx"),
             overwrite = TRUE)

# -------------------------------
# 10. 总体质控输出
# -------------------------------
qc_summary <- data.frame(
  指标 = c(
    "总样本数",
    "Rum样本数",
    "Ile样本数",
    "Col样本数",
    "总特征数（去除全0特征后）",
    "是否存在 metadata 与丰度表样本不一致",
    "是否使用系统发育树",
    "是否计算 UniFrac"
  ),
  数值 = c(
    nrow(meta),
    sum(meta$Site == "Rum"),
    sum(meta$Site == "Ile"),
    sum(meta$Site == "Col"),
    ncol(otu_mat),
    "否",
    "否",
    "否"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.xlsx(qc_summary,
           file.path(out_dir, "04_统计结果", "0_分析质控汇总.xlsx"),
           rowNames = FALSE)

# -------------------------------
# 11. 控制台提示
# -------------------------------
cat("三部位 β 多样性分析已完成。\n")
cat("输出目录：", out_dir, "\n")
cat("主要结果包括：\n")
cat("1. 02_距离矩阵/allsample.bray_curtis_dm.txt\n")
cat("2. 02_距离矩阵/allsample.binary_jaccard_dm.txt\n")
cat("3. 03_PCoA坐标/bray_curtis_PCoA坐标.xlsx\n")
cat("4. 03_PCoA坐标/binary_jaccard_PCoA坐标.xlsx\n")
cat("5. 04_统计结果/bray_curtis_统计结果.xlsx\n")
cat("6. 04_统计结果/binary_jaccard_统计结果.xlsx\n")
cat("7. 05_图形/PCoA_BrayCurtis_3site.pdf\n")
cat("8. 05_图形/PCoA_BinaryJaccard_3site.pdf\n")