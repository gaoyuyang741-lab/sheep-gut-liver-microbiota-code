############################################################

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

############################
## 0. 安装 / 加载包
############################
pkg_needed <- c(
  "readxl", "openxlsx", "dplyr", "tibble", "stringr",
  "purrr", "zCompositions", "compositions"
)

pkg_new <- pkg_needed[!pkg_needed %in% installed.packages()[, "Package"]]
if (length(pkg_new) > 0) {
  install.packages(pkg_new, dependencies = TRUE)
}

invisible(lapply(pkg_needed, library, character.only = TRUE))

############################
## 1. 文件路径
############################
# Input and output directories
# Please place the required input files in "data/pairwise/input".
# Output files will be saved in "results/pairwise".

base_dir <- file.path("data", "pairwise", "input")
out_dir  <- file.path("results", "pairwise")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rum_fp    <- file.path(base_dir, "Rum_genus_abundance.xlsx")
ile_fp    <- file.path(base_dir, "Ile_genus_abundance.xlsx")
col_fp    <- file.path(base_dir, "Colon_genus_abundance.xlsx")
liver_fp  <- file.path(base_dir, "肝脏筛选基因.xlsx")
blood_fp  <- file.path(base_dir, "blood_phenotype_sGCCA_input.xlsx")
tail_fp   <- file.path(base_dir, "尾脂数据_16s编号转换.xlsx")

############################
## 2. 可调参数
############################
prevalence_cut <- 0.50
min_nonzero_n  <- NULL   # 如不需要额外限制，保持 NULL

# 显著性筛选阈值
q_cut   <- 0.05
rho_cut <- 0.30

############################
## 3. 基础函数
############################

# 3.1 样本名统一
clean_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x
}

# 3.2 数值 z-score（零方差列设为0）
zscore_df <- function(df) {
  out <- as.data.frame(df, check.names = FALSE)
  for (j in seq_len(ncol(out))) {
    x <- suppressWarnings(as.numeric(out[[j]]))
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      out[[j]] <- 0
    } else {
      out[[j]] <- (x - m) / s
    }
  }
  out
}

# 3.3 读 genus 丰度表：第一列必须是 Genus，后面是样本列
read_genus_table <- function(fp, sheet = 1) {
  df <- readxl::read_excel(fp, sheet = sheet)
  df <- as.data.frame(df)
  
  if (!"Genus" %in% colnames(df)) {
    stop("文件中未找到 'Genus' 列：", fp)
  }
  
  df <- df[!is.na(df$Genus) & df$Genus != "", , drop = FALSE]
  
  sample_cols <- setdiff(colnames(df), "Genus")
  sample_cols <- clean_sample_id(sample_cols)
  colnames(df) <- c("Genus", sample_cols)
  
  for (cc in sample_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }
  df[is.na(df)] <- 0
  
  # 同名 genus 合并
  df <- df |>
    dplyr::group_by(Genus) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat) <- df$Genus
  mode(mat) <- "numeric"
  
  return(mat)  # 行=Genus，列=Sample
}

# 3.4 prevalence 过滤（按原始丰度 >0 的样本比例）
prevalence_filter <- function(mat, prevalence_cut = 0.30, min_nonzero_n = NULL) {
  prev <- rowMeans(mat > 0, na.rm = TRUE)
  nonzero_n <- rowSums(mat > 0, na.rm = TRUE)
  
  keep <- prev >= prevalence_cut
  if (!is.null(min_nonzero_n)) {
    keep <- keep & (nonzero_n >= min_nonzero_n)
  }
  
  info <- data.frame(
    Feature = rownames(mat),
    prevalence = prev,
    nonzero_n = nonzero_n,
    kept = keep,
    row.names = NULL,
    check.names = FALSE
  )
  
  list(
    mat = mat[keep, , drop = FALSE],
    info = info
  )
}

# 3.5 CZM补零 + CLR
# genus x sample -> t() -> cmultRepl(CZM) -> clr -> 保持 sample x genus
clr_transform_czm <- function(mat_genus_by_sample) {
  if (nrow(mat_genus_by_sample) < 2) {
    stop("过滤后特征数 < 2，无法做 CLR。请放宽 prevalence_cut。")
  }
  
  x <- t(mat_genus_by_sample)  # 行=样本，列=genus
  
  x_nozero <- zCompositions::cmultRepl(
    X = x,
    method = "CZM",
    output = "p-counts",
    label = 0
  )
  
  x_clr <- compositions::clr(x_nozero, base = exp(1))
  x_clr <- as.matrix(x_clr)
  
  return(x_clr)  # 行=样本，列=genus
}

# 3.6 读取 blood
read_blood_traits <- function(fp) {
  df <- readxl::read_excel(fp, sheet = "selected_raw_table")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  need_cols <- c("Sample16S", "TC", "TG", "LDL", "TBA", "GLU")
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("blood 文件缺少列：", paste(miss_cols, collapse = ", "))
  }
  
  out <- df[, need_cols, drop = FALSE]
  colnames(out)[1] <- "SampleID"
  out$SampleID <- clean_sample_id(out$SampleID)
  
  trait_cols <- setdiff(colnames(out), "SampleID")
  for (cc in trait_cols) {
    out[[cc]] <- suppressWarnings(as.numeric(out[[cc]]))
  }
  
  out <- out |>
    dplyr::group_by(SampleID) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out
}

# 3.7 读取 liver formal log2TPM
read_liver_strict <- function(fp) {
  df <- readxl::read_excel(fp, sheet = "formal_log2TPM", skip = 1)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  if (!"Gene" %in% colnames(df)) {
    stop("liver formal_log2TPM 文件中未找到 'Gene' 列。")
  }
  
  df <- df[!is.na(df$Gene) & df$Gene != "", , drop = FALSE]
  
  sample_cols <- setdiff(colnames(df), "Gene")
  sample_cols <- clean_sample_id(sample_cols)
  colnames(df) <- c("Gene", sample_cols)
  
  for (cc in sample_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }
  
  # 同名基因如有重复，取均值
  df <- df |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat) <- df$Gene
  mode(mat) <- "numeric"
  
  # 转成 sample x gene
  mat <- t(mat)
  
  as.data.frame(mat, check.names = FALSE) |>
    tibble::rownames_to_column("SampleID")
}

# 3.8 读取 tailfat（表头从第5行开始）
read_tailfat_traits <- function(fp) {
  df <- readxl::read_excel(fp, sheet = "尾脂_16s编号", skip = 4)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  need_cols <- c("瘤胃和肠道16s编号", "尾脂g", "尾脂/胴体重g/kg", "尾脂/宰前活重g/kg")
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("tailfat 文件缺少列：", paste(miss_cols, collapse = ", "))
  }
  
  out <- df[, need_cols, drop = FALSE]
  colnames(out) <- c("SampleID", "Tailfat_g", "Tailfat_per_Carcass_gkg", "Tailfat_per_LiveWeight_gkg")
  
  out$SampleID <- clean_sample_id(out$SampleID)
  
  trait_cols <- setdiff(colnames(out), "SampleID")
  for (cc in trait_cols) {
    out[[cc]] <- suppressWarnings(as.numeric(out[[cc]]))
  }
  
  out <- out |>
    dplyr::group_by(SampleID) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out
}

# 3.9 对齐两个矩阵（行=样本）
align_two_blocks <- function(df1, df2) {
  s1 <- clean_sample_id(df1$SampleID)
  s2 <- clean_sample_id(df2$SampleID)
  
  common_samples <- intersect(s1, s2)
  common_samples <- sort(common_samples)
  
  if (length(common_samples) == 0) {
    stop("两个 block 没有共同样本。")
  }
  
  df1$SampleID <- s1
  df2$SampleID <- s2
  
  a1 <- df1[match(common_samples, df1$SampleID), , drop = FALSE]
  a2 <- df2[match(common_samples, df2$SampleID), , drop = FALSE]
  
  rownames(a1) <- a1$SampleID
  rownames(a2) <- a2$SampleID
  
  a1 <- a1[, setdiff(colnames(a1), "SampleID"), drop = FALSE]
  a2 <- a2[, setdiff(colnames(a2), "SampleID"), drop = FALSE]
  
  list(
    x = a1,
    y = a2,
    common_samples = common_samples
  )
}

# 3.10 Spearman 两两相关
run_pairwise_spearman <- function(df_x, df_y, module_name) {
  x_names <- colnames(df_x)
  y_names <- colnames(df_y)
  
  res_list <- vector("list", length(x_names) * length(y_names))
  idx <- 1L
  
  for (fx in x_names) {
    x <- suppressWarnings(as.numeric(df_x[[fx]]))
    
    for (fy in y_names) {
      y <- suppressWarnings(as.numeric(df_y[[fy]]))
      
      ok <- complete.cases(x, y)
      n_ok <- sum(ok)
      
      if (n_ok < 3) {
        rho <- NA_real_
        pval <- NA_real_
      } else {
        ct <- suppressWarnings(
          stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
        )
        rho  <- unname(ct$estimate)
        pval <- ct$p.value
      }
      
      res_list[[idx]] <- data.frame(
        module = module_name,
        feature_x = fx,
        feature_y = fy,
        n = n_ok,
        rho = rho,
        p = pval,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  res <- dplyr::bind_rows(res_list)
  res$q <- p.adjust(res$p, method = "BH")
  res$abs_rho <- abs(res$rho)
  res$direction <- ifelse(
    is.na(res$rho), NA,
    ifelse(res$rho > 0, "positive",
           ifelse(res$rho < 0, "negative", "zero"))
  )
  res
}

# 3.11 模块摘要
make_module_summary <- function(res_df, q_cut = 0.05, rho_cut = 0.30) {
  sig <- res_df |>
    dplyr::filter(!is.na(q), !is.na(rho), q < q_cut, abs(rho) >= rho_cut)
  
  data.frame(
    module = unique(res_df$module),
    total_pairs = nrow(res_df),
    tested_pairs = sum(!is.na(res_df$p)),
    significant_pairs = nrow(sig),
    positive_pairs = sum(sig$rho > 0, na.rm = TRUE),
    negative_pairs = sum(sig$rho < 0, na.rm = TRUE),
    unique_feature_x = length(unique(sig$feature_x)),
    unique_feature_y = length(unique(sig$feature_y)),
    stringsAsFactors = FALSE
  )
}

############################
## 4. 读取数据
############################
cat("======================================\n")
cat("读取数据...\n")
cat("======================================\n")

rum_raw_mat <- read_genus_table(rum_fp, sheet = 1)
ile_raw_mat <- read_genus_table(ile_fp, sheet = 1)
col_raw_mat <- read_genus_table(col_fp, sheet = 1)

blood_df_raw   <- read_blood_traits(blood_fp)
liver_df_raw   <- read_liver_strict(liver_fp)
tailfat_df_raw <- read_tailfat_traits(tail_fp)

############################
## 5. 微生物 prevalence过滤 + CZM + CLR + z-score
############################
cat("======================================\n")
cat("微生物预处理：prevalence -> CZM -> CLR -> z-score\n")
cat("======================================\n")

rum_prev <- prevalence_filter(
  rum_raw_mat,
  prevalence_cut = prevalence_cut,
  min_nonzero_n = min_nonzero_n
)

ile_prev <- prevalence_filter(
  ile_raw_mat,
  prevalence_cut = prevalence_cut,
  min_nonzero_n = min_nonzero_n
)

col_prev <- prevalence_filter(
  col_raw_mat,
  prevalence_cut = prevalence_cut,
  min_nonzero_n = min_nonzero_n
)

rum_clr <- clr_transform_czm(rum_prev$mat)
ile_clr <- clr_transform_czm(ile_prev$mat)
col_clr <- clr_transform_czm(col_prev$mat)

rum_clr_z <- zscore_df(as.data.frame(rum_clr, check.names = FALSE)) |>
  tibble::rownames_to_column("SampleID")

ile_clr_z <- zscore_df(as.data.frame(ile_clr, check.names = FALSE)) |>
  tibble::rownames_to_column("SampleID")

col_clr_z <- zscore_df(as.data.frame(col_clr, check.names = FALSE)) |>
  tibble::rownames_to_column("SampleID")

rum_clr_z$SampleID <- clean_sample_id(rum_clr_z$SampleID)
ile_clr_z$SampleID <- clean_sample_id(ile_clr_z$SampleID)
col_clr_z$SampleID <- clean_sample_id(col_clr_z$SampleID)

############################
## 6. liver / blood / tailfat 做 z-score
############################
cat("======================================\n")
cat("宿主层数据 z-score...\n")
cat("======================================\n")

liver_df <- liver_df_raw
blood_df <- blood_df_raw
tailfat_df <- tailfat_df_raw

liver_df$SampleID   <- clean_sample_id(liver_df$SampleID)
blood_df$SampleID   <- clean_sample_id(blood_df$SampleID)
tailfat_df$SampleID <- clean_sample_id(tailfat_df$SampleID)

liver_df_z <- liver_df
liver_df_z[, setdiff(colnames(liver_df_z), "SampleID")] <-
  zscore_df(liver_df_z[, setdiff(colnames(liver_df_z), "SampleID"), drop = FALSE])

blood_df_z <- blood_df
blood_df_z[, setdiff(colnames(blood_df_z), "SampleID")] <-
  zscore_df(blood_df_z[, setdiff(colnames(blood_df_z), "SampleID"), drop = FALSE])

tailfat_df_z <- tailfat_df
tailfat_df_z[, setdiff(colnames(tailfat_df_z), "SampleID")] <-
  zscore_df(tailfat_df_z[, setdiff(colnames(tailfat_df_z), "SampleID"), drop = FALSE])

############################
## 7. QC summary
############################
cat("======================================\n")
cat("生成 QC summary...\n")
cat("======================================\n")

all_sample_sets <- list(
  Rum = rum_clr_z$SampleID,
  Ile = ile_clr_z$SampleID,
  Colon = col_clr_z$SampleID,
  Liver = liver_df_z$SampleID,
  Blood = blood_df_z$SampleID,
  Tailfat = tailfat_df_z$SampleID
)

common_all_6 <- Reduce(intersect, all_sample_sets)

qc_overview <- data.frame(
  Block = c("Rum", "Ile", "Colon", "Liver", "Blood", "Tailfat"),
  n_samples = c(
    length(unique(rum_clr_z$SampleID)),
    length(unique(ile_clr_z$SampleID)),
    length(unique(col_clr_z$SampleID)),
    length(unique(liver_df_z$SampleID)),
    length(unique(blood_df_z$SampleID)),
    length(unique(tailfat_df_z$SampleID))
  ),
  n_features = c(
    ncol(rum_clr_z) - 1,
    ncol(ile_clr_z) - 1,
    ncol(col_clr_z) - 1,
    ncol(liver_df_z) - 1,
    ncol(blood_df_z) - 1,
    ncol(tailfat_df_z) - 1
  ),
  stringsAsFactors = FALSE
)

microbe_filter_summary <- data.frame(
  Block = c("Rum", "Ile", "Colon"),
  raw_features = c(nrow(rum_raw_mat), nrow(ile_raw_mat), nrow(col_raw_mat)),
  kept_features = c(nrow(rum_prev$mat), nrow(ile_prev$mat), nrow(col_prev$mat)),
  prevalence_cut = prevalence_cut,
  min_nonzero_n = ifelse(is.null(min_nonzero_n), NA, min_nonzero_n),
  stringsAsFactors = FALSE
)

all_block_common_summary <- data.frame(
  Metric = c("Common samples across Rum/Ile/Colon/Liver/Blood/Tailfat"),
  Value = c(length(common_all_6)),
  stringsAsFactors = FALSE
)

############################
## 8. 11个模块对齐并分析
############################
cat("======================================\n")
cat("开始 11 个模块 Spearman 相关分析...\n")
cat("======================================\n")

module_def <- list(
  list(name = "Rum_vs_Liver",      x = rum_clr_z,   y = liver_df_z),
  list(name = "Ile_vs_Liver",      x = ile_clr_z,   y = liver_df_z),
  list(name = "Colon_vs_Liver",    x = col_clr_z,   y = liver_df_z),
  
  list(name = "Rum_vs_Blood",      x = rum_clr_z,   y = blood_df_z),
  list(name = "Ile_vs_Blood",      x = ile_clr_z,   y = blood_df_z),
  list(name = "Colon_vs_Blood",    x = col_clr_z,   y = blood_df_z),
  
  list(name = "Rum_vs_Tailfat",    x = rum_clr_z,   y = tailfat_df_z),
  list(name = "Ile_vs_Tailfat",    x = ile_clr_z,   y = tailfat_df_z),
  list(name = "Colon_vs_Tailfat",  x = col_clr_z,   y = tailfat_df_z),
  
  list(name = "Liver_vs_Blood",    x = liver_df_z,  y = blood_df_z),
  list(name = "Liver_vs_Tailfat",  x = liver_df_z,  y = tailfat_df_z)
)

all_results <- list()
sig_results <- list()
module_summary_list <- list()
module_common_n <- list()

for (i in seq_along(module_def)) {
  m <- module_def[[i]]
  
  cat("\n--------------------------------------\n")
  cat("Module:", m$name, "\n")
  cat("--------------------------------------\n")
  
  aligned <- align_two_blocks(m$x, m$y)
  
  df_x <- aligned$x
  df_y <- aligned$y
  common_samples <- aligned$common_samples
  
  cat("Common samples:", length(common_samples), "\n")
  cat("n_feature_x:", ncol(df_x), "\n")
  cat("n_feature_y:", ncol(df_y), "\n")
  
  res <- run_pairwise_spearman(df_x, df_y, module_name = m$name)
  
  sig <- res |>
    dplyr::filter(!is.na(q), !is.na(rho), q < q_cut, abs(rho) >= rho_cut) |>
    dplyr::arrange(q, dplyr::desc(abs_rho))
  
  all_results[[m$name]] <- res
  sig_results[[m$name]] <- sig
  module_summary_list[[m$name]] <- make_module_summary(res, q_cut = q_cut, rho_cut = rho_cut)
  
  module_common_n[[m$name]] <- data.frame(
    module = m$name,
    n_common_samples = length(common_samples),
    feature_x_n = ncol(df_x),
    feature_y_n = ncol(df_y),
    stringsAsFactors = FALSE
  )
}

all_results_df   <- dplyr::bind_rows(all_results)
sig_results_df   <- dplyr::bind_rows(sig_results)
module_summary_df <- dplyr::bind_rows(module_summary_list)
module_common_n_df <- dplyr::bind_rows(module_common_n)

############################
## 9. 导出 feature 名单 / 预处理矩阵
############################
rum_feature_list <- data.frame(
  Block = "Rum",
  Feature = colnames(rum_clr_z)[colnames(rum_clr_z) != "SampleID"],
  stringsAsFactors = FALSE
)

ile_feature_list <- data.frame(
  Block = "Ile",
  Feature = colnames(ile_clr_z)[colnames(ile_clr_z) != "SampleID"],
  stringsAsFactors = FALSE
)

col_feature_list <- data.frame(
  Block = "Colon",
  Feature = colnames(col_clr_z)[colnames(col_clr_z) != "SampleID"],
  stringsAsFactors = FALSE
)

liver_feature_list <- data.frame(
  Block = "Liver",
  Feature = colnames(liver_df_z)[colnames(liver_df_z) != "SampleID"],
  stringsAsFactors = FALSE
)

blood_feature_list <- data.frame(
  Block = "Blood",
  Feature = colnames(blood_df_z)[colnames(blood_df_z) != "SampleID"],
  stringsAsFactors = FALSE
)

tailfat_feature_list <- data.frame(
  Block = "Tailfat",
  Feature = colnames(tailfat_df_z)[colnames(tailfat_df_z) != "SampleID"],
  stringsAsFactors = FALSE
)

feature_list_df <- dplyr::bind_rows(
  rum_feature_list, ile_feature_list, col_feature_list,
  liver_feature_list, blood_feature_list, tailfat_feature_list
)

############################
## 10. 保存 Excel 结果
############################
cat("======================================\n")
cat("写出 Excel 结果...\n")
cat("======================================\n")

wb <- openxlsx::createWorkbook()

## 10.1 QC
openxlsx::addWorksheet(wb, "QC_overview")
openxlsx::writeData(wb, "QC_overview", qc_overview)

openxlsx::addWorksheet(wb, "Microbe_filter_summary")
openxlsx::writeData(wb, "Microbe_filter_summary", microbe_filter_summary)

openxlsx::addWorksheet(wb, "All_block_common")
openxlsx::writeData(wb, "All_block_common", all_block_common_summary)

openxlsx::addWorksheet(wb, "Module_common_samples")
openxlsx::writeData(wb, "Module_common_samples", module_common_n_df)

## 10.2 feature 名单
openxlsx::addWorksheet(wb, "Feature_list")
openxlsx::writeData(wb, "Feature_list", feature_list_df)

openxlsx::addWorksheet(wb, "Liver_KEGG_gene_list")
openxlsx::writeData(wb, "Liver_KEGG_gene_list", liver_feature_list)

## 10.3 prevalence 信息
openxlsx::addWorksheet(wb, "Rum_prevalence_info")
openxlsx::writeData(wb, "Rum_prevalence_info", rum_prev$info)

openxlsx::addWorksheet(wb, "Ile_prevalence_info")
openxlsx::writeData(wb, "Ile_prevalence_info", ile_prev$info)

openxlsx::addWorksheet(wb, "Colon_prevalence_info")
openxlsx::writeData(wb, "Colon_prevalence_info", col_prev$info)

## 10.4 预处理矩阵
openxlsx::addWorksheet(wb, "Rum_CLR_z")
openxlsx::writeData(wb, "Rum_CLR_z", rum_clr_z)

openxlsx::addWorksheet(wb, "Ile_CLR_z")
openxlsx::writeData(wb, "Ile_CLR_z", ile_clr_z)

openxlsx::addWorksheet(wb, "Colon_CLR_z")
openxlsx::writeData(wb, "Colon_CLR_z", col_clr_z)

openxlsx::addWorksheet(wb, "Liver_z")
openxlsx::writeData(wb, "Liver_z", liver_df_z)

openxlsx::addWorksheet(wb, "Blood_z")
openxlsx::writeData(wb, "Blood_z", blood_df_z)

openxlsx::addWorksheet(wb, "Tailfat_z")
openxlsx::writeData(wb, "Tailfat_z", tailfat_df_z)

## 10.5 模块摘要
openxlsx::addWorksheet(wb, "Module_summary")
openxlsx::writeData(wb, "Module_summary", module_summary_df)

## 10.6 全部相关结果
openxlsx::addWorksheet(wb, "All_associations")
openxlsx::writeData(wb, "All_associations", all_results_df)

## 10.7 显著结果
openxlsx::addWorksheet(wb, "Significant_associations")
openxlsx::writeData(wb, "Significant_associations", sig_results_df)

## 10.8 各模块单独表
for (nm in names(all_results)) {
  sheet_nm <- paste0(substr(nm, 1, 28), "_all")
  sheet_nm <- gsub("[\\\\/:?*\\[\\]]", "_", sheet_nm)
  openxlsx::addWorksheet(wb, sheet_nm)
  openxlsx::writeData(wb, sheet_nm, all_results[[nm]])
}

for (nm in names(sig_results)) {
  sheet_nm <- paste0(substr(nm, 1, 28), "_sig")
  sheet_nm <- gsub("[\\\\/:?*\\[\\]]", "_", sheet_nm)
  openxlsx::addWorksheet(wb, sheet_nm)
  openxlsx::writeData(wb, sheet_nm, sig_results[[nm]])
}

## 基础样式
header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sh in openxlsx::sheets(wb)) {
  nc <- tryCatch(ncol(openxlsx::readWorkbook(wb, sheet = sh)), error = function(e) 0)
  if (nc > 0) {
    openxlsx::addStyle(
      wb, sh, style = header_style,
      rows = 1, cols = 1:nc, gridExpand = TRUE
    )
    openxlsx::setColWidths(wb, sh, cols = 1:nc, widths = "auto")
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
  }
}

out_xlsx <- file.path(out_dir, "sGCCA_prelude_Spearman_11modules_results.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

############################
## 11. 另存几个常用单文件
############################
openxlsx::write.xlsx(
  x = liver_feature_list,
  file = file.path(out_dir, "Liver_KEGG_gene_list.xlsx"),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  x = module_summary_df,
  file = file.path(out_dir, "Module_summary.xlsx"),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  x = sig_results_df,
  file = file.path(out_dir, "Significant_associations.xlsx"),
  overwrite = TRUE
)

############################
## 12. 控制台输出
############################
cat("\n======================================\n")
cat("运行完成。\n")
cat("结果目录：", out_dir, "\n")
cat("总结果文件：", out_xlsx, "\n")
cat("显著阈值：q <", q_cut, "且 |rho| >=", rho_cut, "\n")
cat("Rum 保留 genus 数：", nrow(rum_prev$mat), "\n")
cat("Ile 保留 genus 数：", nrow(ile_prev$mat), "\n")
cat("Colon 保留 genus 数：", nrow(col_prev$mat), "\n")
cat("Liver KEGG gene 数：", ncol(liver_df_z) - 1, "\n")
cat("Blood trait 数：", ncol(blood_df_z) - 1, "\n")
cat("Tailfat trait 数：", ncol(tailfat_df_z) - 1, "\n")
cat("六个 block 共同样本数：", length(common_all_6), "\n")
cat("======================================\n")