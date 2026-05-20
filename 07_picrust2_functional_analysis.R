# =========================================================
# PICRUSt2 三部位“功能梯度与功能分工”基础数据处理脚本（精简输出版）
# 说明：
# 1）无用的过程表全部合并到一个 Excel
# 2）真正有用的结果表全部合并到一个 Excel
# 3）不作图，只输出基础处理与统计整理结果
# =========================================================

rm(list = ls())
gc()

# =========================================================
# 0. 加载程序包
# =========================================================
pkg_needed <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr",
  "readxl", "openxlsx", "tibble"
)

pkg_to_install <- pkg_needed[!pkg_needed %in% rownames(installed.packages())]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, dependencies = TRUE)
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(datatable.fread.datatable = FALSE)

# =========================================================
# 1. 参数设置
# =========================================================
# Input and output directories
# Please place the required input files in "data/picrust2/input".
# Output files will be saved in "results/picrust2".

input_dir <- file.path("data", "picrust2", "input")
out_dir   <- file.path("results", "picrust2")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 输入文件
ri_ko_file        <- file.path(input_dir, "RI_KO_pred_metagenome_unstrat.tsv")
c_ko_file         <- file.path(input_dir, "C_KO_pred_metagenome_unstrat.tsv")
metadata_file     <- file.path(input_dir, "metadata_3group.xlsx")
ko2path_file      <- file.path(input_dir, "ko2pathway.txt")
path_name_file    <- file.path(input_dir, "pathway_list_ko.txt")
brite_file        <- file.path(input_dir, "br08901.txt")

# 元数据列名
sample_col <- "SampleID"
group_col  <- "Group"
sheep_col  <- "SheepID"

# 组别顺序
group_levels <- c("Rum", "Ile", "Col")

# 过滤参数
use_relative_abundance <- TRUE
prev_cut_ko            <- 0.10
mean_cut_ko            <- 1e-6
prev_cut_pathway       <- 0.10
mean_cut_pathway       <- 1e-6
prev_cut_l2            <- 0
mean_cut_l2            <- 0
prev_cut_l1            <- 0
mean_cut_l1            <- 0

# 统计参数
p_adj_method <- "BH"
sig_cut      <- 0.001

# 伪计数（预留）
pseudo <- 1e-12

# =========================================================
# 2. 通用函数
# =========================================================
clean_group <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("Rum", "Rumen", "rum", "rumen")] <- "Rum"
  x[x %in% c("Ile", "Ileum", "ile", "ileum")] <- "Ile"
  x[x %in% c("Col", "Colon", "colon", "col")] <- "Col"
  x
}

clean_ko <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^ko:", "", x, ignore.case = TRUE)
  x <- sub("^KO:", "", x, ignore.case = TRUE)
  x <- sub("^path:", "", x, ignore.case = TRUE)
  x <- sub("^map", "", x, ignore.case = TRUE)
  x <- sub("^ko", "", x, ignore.case = TRUE)
  x <- paste0("K", x)
  x <- sub("^KK", "K", x)
  x
}

clean_pathway <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^path:", "", x, ignore.case = TRUE)
  x <- sub("^map", "", x, ignore.case = TRUE)
  x <- sub("^ko", "", x, ignore.case = TRUE)
  paste0("ko", x)
}

safe_numeric_matrix <- function(df) {
  m <- as.matrix(df)
  mode(m) <- "numeric"
  m[is.na(m)] <- 0
  m
}

calc_prev_mean <- function(mat) {
  data.frame(
    feature = rownames(mat),
    prevalence = rowMeans(mat > 0),
    mean_abundance = rowMeans(mat),
    stringsAsFactors = FALSE
  )
}

filter_feature_matrix <- function(mat, prev_cut = 0, mean_cut = 0) {
  stat <- calc_prev_mean(mat)
  keep <- stat$prevalence >= prev_cut & stat$mean_abundance >= mean_cut
  mat[keep, , drop = FALSE]
}

median_by_group <- function(vec, meta) {
  tapply(vec, meta$Group, median, na.rm = TRUE)
}

mean_by_group <- function(vec, meta) {
  tapply(vec, meta$Group, mean, na.rm = TRUE)
}

paired_test_one_feature <- function(vec, meta) {
  df <- data.frame(
    SheepID = meta$SheepID,
    Group   = meta$Group,
    value   = as.numeric(vec),
    stringsAsFactors = FALSE
  )
  
  wide <- df %>%
    tidyr::pivot_wider(names_from = Group, values_from = value) %>%
    dplyr::select(all_of(c("SheepID", group_levels)))
  
  wide_complete <- wide %>% tidyr::drop_na(all_of(group_levels))
  
  overall_p <- NA_real_
  if (nrow(wide_complete) >= 3) {
    overall_p <- tryCatch(
      stats::friedman.test(as.matrix(wide_complete[, group_levels]))$p.value,
      error = function(e) NA_real_
    )
  }
  
  pair_p <- function(g1, g2) {
    tmp <- wide %>%
      dplyr::select(SheepID, all_of(c(g1, g2))) %>%
      tidyr::drop_na(all_of(c(g1, g2)))
    if (nrow(tmp) < 3) return(NA_real_)
    tryCatch(
      stats::wilcox.test(tmp[[g1]], tmp[[g2]], paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  }
  
  med <- median_by_group(vec, meta)
  mn  <- mean_by_group(vec, meta)
  
  data.frame(
    overall_p = overall_p,
    p_Rum_Ile = pair_p("Rum", "Ile"),
    p_Rum_Col = pair_p("Rum", "Col"),
    p_Ile_Col = pair_p("Ile", "Col"),
    median_Rum = unname(ifelse("Rum" %in% names(med), med["Rum"], NA)),
    median_Ile = unname(ifelse("Ile" %in% names(med), med["Ile"], NA)),
    median_Col = unname(ifelse("Col" %in% names(med), med["Col"], NA)),
    mean_Rum = unname(ifelse("Rum" %in% names(mn), mn["Rum"], NA)),
    mean_Ile = unname(ifelse("Ile" %in% names(mn), mn["Ile"], NA)),
    mean_Col = unname(ifelse("Col" %in% names(mn), mn["Col"], NA)),
    stringsAsFactors = FALSE
  )
}

classify_pattern <- function(mR, mI, mC, qRI, qRC, qIC, sig_cut = 0.001) {
  if (any(is.na(c(mR, mI, mC)))) return("Other")
  
  sigRI <- !is.na(qRI) && qRI < sig_cut
  sigRC <- !is.na(qRC) && qRC < sig_cut
  sigIC <- !is.na(qIC) && qIC < sig_cut
  
  if ((!sigRC) && sigRI && sigIC) {
    if ((mR > mI && mC > mI) || (mR < mI && mC < mI)) return("FHC")
  }
  
  if (sigRI && sigIC) {
    if ((mI > mR && mI > mC) || (mI < mR && mI < mC)) return("IT")
  }
  
  if ((mR > mI && mI > mC && sigRI && sigIC) ||
      (mR < mI && mI < mC && sigRI && sigIC)) {
    return("Grad")
  }
  
  if (mR > mI && mR > mC && sigRI && sigRC) return("FGH")
  if (mC > mR && mC > mI && sigRC && sigIC) return("HGH")
  
  "Other"
}

run_feature_stats <- function(mat, meta, feature_info = NULL, level_name = "feature") {
  res_list <- lapply(seq_len(nrow(mat)), function(i) {
    paired_test_one_feature(mat[i, ], meta)
  })
  
  res <- bind_rows(res_list)
  res <- cbind(feature = rownames(mat), res, stringsAsFactors = FALSE)
  
  res$overall_q <- p.adjust(res$overall_p, method = p_adj_method)
  res$q_Rum_Ile <- p.adjust(res$p_Rum_Ile, method = p_adj_method)
  res$q_Rum_Col <- p.adjust(res$p_Rum_Col, method = p_adj_method)
  res$q_Ile_Col <- p.adjust(res$p_Ile_Col, method = p_adj_method)
  
  res$pattern_class <- purrr::pmap_chr(
    list(res$median_Rum, res$median_Ile, res$median_Col,
         res$q_Rum_Ile, res$q_Rum_Col, res$q_Ile_Col),
    ~ classify_pattern(..1, ..2, ..3, ..4, ..5, ..6, sig_cut = sig_cut)
  )
  
  stat_basic <- calc_prev_mean(mat)
  colnames(stat_basic) <- c("feature", "prevalence", "mean_abundance_all")
  res <- left_join(res, stat_basic, by = "feature")
  
  if (!is.null(feature_info)) {
    res <- left_join(feature_info, res, by = "feature")
  }
  
  res$level <- level_name
  
  res <- res %>%
    dplyr::relocate(level, .before = 1) %>%
    arrange(overall_q, q_Rum_Col, q_Ile_Col, q_Rum_Ile, feature)
  
  res
}

top_n_each_pattern <- function(df, n = 20) {
  df %>%
    filter(pattern_class != "Other") %>%
    group_by(pattern_class) %>%
    arrange(overall_q, .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup()
}

write_multi_sheet_xlsx <- function(file, sheet_list) {
  wb <- createWorkbook()
  for (nm in names(sheet_list)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheet_list[[nm]])
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

# =========================================================
# 3. 读取 metadata
# =========================================================
meta <- readxl::read_excel(metadata_file)
meta <- as.data.frame(meta)

need_cols <- c(sample_col, group_col, sheep_col)
if (!all(need_cols %in% colnames(meta))) {
  stop("metadata 缺少必要列：", paste(setdiff(need_cols, colnames(meta)), collapse = ", "))
}

meta <- meta %>%
  dplyr::rename(
    SampleID = all_of(sample_col),
    Group    = all_of(group_col),
    SheepID  = all_of(sheep_col)
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    Group    = clean_group(Group),
    SheepID  = as.character(SheepID)
  ) %>%
  filter(Group %in% group_levels)

meta$Group <- factor(meta$Group, levels = group_levels)

sheep_keep <- meta %>%
  count(SheepID, Group) %>%
  tidyr::pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Rum >= 1, Ile >= 1, Col >= 1) %>%
  pull(SheepID)

meta <- meta %>%
  filter(SheepID %in% sheep_keep) %>%
  arrange(SheepID, Group)

# =========================================================
# 4. 读取 KO 表
# =========================================================
read_ko_table <- function(file, meta_sample_ids) {
  dat <- data.table::fread(file, sep = "\t", header = TRUE, check.names = FALSE)
  dat <- as.data.frame(dat)
  
  colnames(dat)[1] <- "feature_raw"
  dat$feature <- clean_ko(dat$feature_raw)
  
  sample_cols <- intersect(colnames(dat), meta_sample_ids)
  if (length(sample_cols) == 0) {
    stop("文件中未找到与 metadata 匹配的样本列：", basename(file))
  }
  
  dat2 <- dat[, c("feature", sample_cols), drop = FALSE]
  dat2 <- dat2 %>%
    group_by(feature) %>%
    summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  
  dat2
}

ri_ko <- read_ko_table(ri_ko_file, meta$SampleID)
c_ko  <- read_ko_table(c_ko_file,  meta$SampleID)

ko_all <- full_join(ri_ko, c_ko, by = "feature")
ko_all[is.na(ko_all)] <- 0

sample_cols_final <- meta$SampleID
sample_cols_final <- sample_cols_final[sample_cols_final %in% colnames(ko_all)]
ko_all <- ko_all[, c("feature", sample_cols_final), drop = FALSE]

ko_mat <- safe_numeric_matrix(ko_all[, -1, drop = FALSE])
rownames(ko_mat) <- ko_all$feature

if (use_relative_abundance) {
  col_sums <- colSums(ko_mat, na.rm = TRUE)
  col_sums[col_sums == 0] <- 1
  ko_mat <- sweep(ko_mat, 2, col_sums, "/")
}

ko_mat_filt <- filter_feature_matrix(ko_mat, prev_cut = prev_cut_ko, mean_cut = mean_cut_ko)

# =========================================================
# 5. 注释与映射
# =========================================================
ko2path <- data.table::fread(ko2path_file, sep = "\t", header = FALSE)
ko2path <- as.data.frame(ko2path)
colnames(ko2path)[1:2] <- c("ko_raw", "path_raw")

ko2path <- ko2path %>%
  mutate(
    feature = clean_ko(ko_raw),
    pathway = clean_pathway(path_raw)
  ) %>%
  select(feature, pathway) %>%
  distinct()

path_name <- data.table::fread(path_name_file, sep = "\t", header = FALSE)
path_name <- as.data.frame(path_name)
colnames(path_name)[1:2] <- c("pathway", "pathway_name")

path_name <- path_name %>%
  mutate(pathway = clean_pathway(pathway)) %>%
  distinct()

brite_lines <- readLines(brite_file, warn = FALSE, encoding = "UTF-8")
brite_lines <- gsub("\t", " ", brite_lines)

current_A <- NA_character_
current_B <- NA_character_
brite_map <- list()

for (ln in brite_lines) {
  ln2 <- trimws(ln)
  if (ln2 == "" || ln2 %in% c("!", "#")) next
  
  if (grepl("^A", ln2)) {
    current_A <- trimws(sub("^A", "", ln2))
  } else if (grepl("^B", ln2)) {
    current_B <- trimws(sub("^B", "", ln2))
  } else if (grepl("^C\\s+[0-9]{5}", ln2)) {
    tmp <- sub("^C\\s+([0-9]{5})\\s+", "\\1\t", ln2)
    sp <- strsplit(tmp, "\t")[[1]]
    pid <- paste0("ko", sp[1])
    pname <- ifelse(length(sp) >= 2, sp[2], NA_character_)
    brite_map[[length(brite_map) + 1]] <- data.frame(
      pathway = pid,
      pathway_name_brite = pname,
      L1 = current_A,
      L2 = current_B,
      stringsAsFactors = FALSE
    )
  }
}

brite_map <- bind_rows(brite_map) %>% distinct()

path_info <- full_join(path_name, brite_map, by = "pathway") %>%
  mutate(pathway_name = coalesce(pathway_name, pathway_name_brite)) %>%
  select(pathway, pathway_name, L1, L2) %>%
  distinct()

ko_info_full <- ko2path %>%
  left_join(path_info, by = "pathway") %>%
  group_by(feature) %>%
  summarise(
    pathway_n = n_distinct(pathway),
    pathway = paste(sort(unique(pathway)), collapse = "; "),
    pathway_name = paste(sort(unique(na.omit(pathway_name))), collapse = "; "),
    L1 = paste(sort(unique(na.omit(L1))), collapse = "; "),
    L2 = paste(sort(unique(na.omit(L2))), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    in_filtered_KO = feature %in% rownames(ko_mat_filt)
  ) %>%
  arrange(feature)

# =========================================================
# 6. KO 聚合到 pathway / L2 / L1
# =========================================================
ko2path_use <- ko2path %>%
  filter(feature %in% rownames(ko_mat_filt)) %>%
  distinct()

path_list <- split(ko2path_use$feature, ko2path_use$pathway)
path_list <- path_list[lengths(path_list) > 0]

path_mat <- t(sapply(path_list, function(kos) {
  colSums(ko_mat_filt[intersect(kos, rownames(ko_mat_filt)), , drop = FALSE], na.rm = TRUE)
}))
path_mat <- as.matrix(path_mat)
mode(path_mat) <- "numeric"

if (is.null(rownames(path_mat))) {
  rownames(path_mat) <- names(path_list)
}

path_mat_filt <- filter_feature_matrix(path_mat, prev_cut = prev_cut_pathway, mean_cut = mean_cut_pathway)

path_info_use <- path_info %>%
  filter(pathway %in% rownames(path_mat_filt), !is.na(L2), L2 != "")

l2_list <- split(path_info_use$pathway, path_info_use$L2)
l2_list <- l2_list[lengths(l2_list) > 0]

l2_mat <- t(sapply(l2_list, function(pths) {
  colSums(path_mat_filt[intersect(pths, rownames(path_mat_filt)), , drop = FALSE], na.rm = TRUE)
}))
l2_mat <- as.matrix(l2_mat)
mode(l2_mat) <- "numeric"
l2_mat_filt <- filter_feature_matrix(l2_mat, prev_cut = prev_cut_l2, mean_cut = mean_cut_l2)

path_info_use2 <- path_info %>%
  filter(pathway %in% rownames(path_mat_filt), !is.na(L1), L1 != "")

l1_list <- split(path_info_use2$pathway, path_info_use2$L1)
l1_list <- l1_list[lengths(l1_list) > 0]

l1_mat <- t(sapply(l1_list, function(pths) {
  colSums(path_mat_filt[intersect(pths, rownames(path_mat_filt)), , drop = FALSE], na.rm = TRUE)
}))
l1_mat <- as.matrix(l1_mat)
mode(l1_mat) <- "numeric"
l1_mat_filt <- filter_feature_matrix(l1_mat, prev_cut = prev_cut_l1, mean_cut = mean_cut_l1)

# =========================================================
# 7. 统一 metadata 顺序
# =========================================================
meta_use <- meta %>% filter(SampleID %in% colnames(ko_mat_filt))
meta_use <- meta_use[match(colnames(ko_mat_filt), meta_use$SampleID), ]
stopifnot(all(meta_use$SampleID == colnames(ko_mat_filt)))

# =========================================================
# 8. 各层级统计
# =========================================================
ko_info_used <- ko_info_full %>% filter(feature %in% rownames(ko_mat_filt))

ko_res <- run_feature_stats(
  mat = ko_mat_filt,
  meta = meta_use,
  feature_info = ko_info_used,
  level_name = "KO"
)

path_feature_info <- path_info %>% rename(feature = pathway)

path_res <- run_feature_stats(
  mat = path_mat_filt,
  meta = meta_use,
  feature_info = path_feature_info,
  level_name = "Pathway"
)

l2_info <- data.frame(
  feature = rownames(l2_mat_filt),
  L2 = rownames(l2_mat_filt),
  stringsAsFactors = FALSE
)

l2_res <- run_feature_stats(
  mat = l2_mat_filt,
  meta = meta_use,
  feature_info = l2_info,
  level_name = "L2"
)

l1_info <- data.frame(
  feature = rownames(l1_mat_filt),
  L1 = rownames(l1_mat_filt),
  stringsAsFactors = FALSE
)

l1_res <- run_feature_stats(
  mat = l1_mat_filt,
  meta = meta_use,
  feature_info = l1_info,
  level_name = "L1"
)

# 显著结果
ko_sig   <- ko_res   %>% filter(overall_q < sig_cut)
path_sig <- path_res %>% filter(overall_q < sig_cut)
l2_sig   <- l2_res   %>% filter(overall_q < sig_cut)
l1_sig   <- l1_res   %>% filter(overall_q < sig_cut)

# 每类模式代表条目
ko_top_pattern   <- top_n_each_pattern(ko_sig, n = 20)
path_top_pattern <- top_n_each_pattern(path_sig, n = 20)
l2_top_pattern   <- top_n_each_pattern(l2_sig, n = 20)
l1_top_pattern   <- top_n_each_pattern(l1_sig, n = 20)

# =========================================================
# 9. 汇总表
# =========================================================
pattern_summary <- bind_rows(
  ko_res   %>% count(level, pattern_class, name = "n"),
  path_res %>% count(level, pattern_class, name = "n"),
  l2_res   %>% count(level, pattern_class, name = "n"),
  l1_res   %>% count(level, pattern_class, name = "n")
) %>%
  arrange(level, match(pattern_class, c("FGH", "FHC", "Grad", "HGH", "IT", "Other")))

sig_summary <- bind_rows(
  ko_res   %>% summarise(level = "KO",      n_total = n(), n_sig = sum(overall_q < sig_cut, na.rm = TRUE)),
  path_res %>% summarise(level = "Pathway", n_total = n(), n_sig = sum(overall_q < sig_cut, na.rm = TRUE)),
  l2_res   %>% summarise(level = "L2",      n_total = n(), n_sig = sum(overall_q < sig_cut, na.rm = TRUE)),
  l1_res   %>% summarise(level = "L1",      n_total = n(), n_sig = sum(overall_q < sig_cut, na.rm = TRUE))
)

ko_median_table <- ko_res %>%
  select(
    feature, pathway_n, pathway, pathway_name, L1, L2,
    median_Rum, median_Ile, median_Col,
    mean_Rum, mean_Ile, mean_Col,
    overall_p, overall_q, p_Rum_Ile, q_Rum_Ile, p_Rum_Col, q_Rum_Col, p_Ile_Col, q_Ile_Col,
    prevalence, mean_abundance_all, pattern_class
  ) %>%
  arrange(overall_q, pattern_class, feature)

path_median_table <- path_res %>%
  select(
    feature, pathway_name, L1, L2,
    median_Rum, median_Ile, median_Col,
    mean_Rum, mean_Ile, mean_Col,
    overall_p, overall_q, p_Rum_Ile, q_Rum_Ile, p_Rum_Col, q_Rum_Col, p_Ile_Col, q_Ile_Col,
    prevalence, mean_abundance_all, pattern_class
  ) %>%
  arrange(overall_q, pattern_class, feature)

# =========================================================
# 10. 检查表
# =========================================================
sample_check <- meta_use %>%
  count(Group, name = "n_samples") %>%
  mutate(n_sheep = length(unique(meta_use$SheepID)))

file_check <- data.frame(
  item = c(
    "metadata_used_samples",
    "unique_sheep",
    "KO_before_filter",
    "KO_after_filter",
    "pathway_after_filter",
    "L2_after_filter",
    "L1_after_filter"
  ),
  value = c(
    nrow(meta_use),
    length(unique(meta_use$SheepID)),
    nrow(ko_mat),
    nrow(ko_mat_filt),
    nrow(path_mat_filt),
    nrow(l2_mat_filt),
    nrow(l1_mat_filt)
  ),
  stringsAsFactors = FALSE
)

filter_params <- data.frame(
  parameter = c(
    "use_relative_abundance",
    "prev_cut_ko", "mean_cut_ko",
    "prev_cut_pathway", "mean_cut_pathway",
    "prev_cut_l2", "mean_cut_l2",
    "prev_cut_l1", "mean_cut_l1",
    "p_adj_method", "sig_cut", "pseudo"
  ),
  value = c(
    as.character(use_relative_abundance),
    prev_cut_ko, mean_cut_ko,
    prev_cut_pathway, mean_cut_pathway,
    prev_cut_l2, mean_cut_l2,
    prev_cut_l1, mean_cut_l1,
    p_adj_method, sig_cut, pseudo
  ),
  stringsAsFactors = FALSE
)

# =========================================================
# 11. 输出：过程表全部塞进一个 Excel
# =========================================================
process_workbook <- file.path(out_dir, "PICRUSt2_过程表_合并.xlsx")

write_multi_sheet_xlsx(
  process_workbook,
  list(
    metadata_used = meta_use,
    check_summary = file_check,
    sample_summary = sample_check,
    parameter_summary = filter_params,
    KO_annotation_full = ko_info_full,
    pathway_annotation = path_info,
    KO_abundance_filtered = data.frame(feature = rownames(ko_mat_filt), ko_mat_filt, check.names = FALSE),
    pathway_abundance_filtered = data.frame(feature = rownames(path_mat_filt), path_mat_filt, check.names = FALSE),
    L2_abundance_filtered = data.frame(feature = rownames(l2_mat_filt), l2_mat_filt, check.names = FALSE),
    L1_abundance_filtered = data.frame(feature = rownames(l1_mat_filt), l1_mat_filt, check.names = FALSE)
  )
)

# =========================================================
# 12. 输出：有用结果全部塞进一个 Excel
# =========================================================
result_workbook <- file.path(out_dir, "PICRUSt2_核心结果_合并.xlsx")

write_multi_sheet_xlsx(
  result_workbook,
  list(
    summary_significance = sig_summary,
    summary_pattern = pattern_summary,
    
    KO_all_results = ko_res,
    KO_sig_results = ko_sig,
    KO_top_each_pattern = ko_top_pattern,
    KO_group_median_table = ko_median_table,
    
    pathway_all_results = path_res,
    pathway_sig_results = path_sig,
    pathway_top_each_pattern = path_top_pattern,
    pathway_group_median_table = path_median_table,
    
    L2_all_results = l2_res,
    L2_sig_results = l2_sig,
    L2_top_each_pattern = l2_top_pattern,
    
    L1_all_results = l1_res,
    L1_sig_results = l1_sig,
    L1_top_each_pattern = l1_top_pattern
  )
)

# =========================================================
# 13. 结束
# =========================================================
cat("\n基础数据处理完成（精简输出版）。\n")
cat("输出目录：", out_dir, "\n")
cat("过程表：", process_workbook, "\n")
cat("核心结果：", result_workbook, "\n")