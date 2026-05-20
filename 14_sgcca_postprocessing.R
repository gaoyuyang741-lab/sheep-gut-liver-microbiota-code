###############################################################
###############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)

# ---- avoid namespace masking ----
select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
arrange    <- dplyr::arrange
rename     <- dplyr::rename
summarise  <- dplyr::summarise
distinct   <- dplyr::distinct
pull       <- dplyr::pull
bind_rows  <- dplyr::bind_rows
bind_cols  <- dplyr::bind_cols
left_join  <- dplyr::left_join
inner_join <- dplyr::inner_join

# ============================================================
# 一、路径参数
# ============================================================
# Input and output directories
# This script reads sGCCA output files from "results/sgcca/sgcca_5block_8traits_out".
# Refined sGCCA results will be saved in "results/sgcca/refined_results".

input_root  <- file.path("results", "sgcca", "sgcca_5block_8traits_out")
output_root <- file.path("results", "sgcca", "refined_results")

dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
# ============================================================
# 二、后处理标准参数区
# ============================================================
params <- list(
  # ---------- 第一层：模型是否成立 ----------
  ave_inner_a = 0.40,
  ave_inner_b = 0.30,
  
  corr_trait_liver = 0.40,
  corr_other_pairs = 0.35,
  min_other_pairs_n = 2,
  
  # ---------- 第二层：主链由谁组成 ----------
  feature_method = "top_n",
  
  top_n_rum   = 10,
  top_n_ile   = 10,
  top_n_col   = 10,
  top_n_liver = 10,
  top_n_trait = 1,
  
  quantile_rum   = 0.80,
  quantile_ile   = 0.80,
  quantile_col   = 0.80,
  quantile_liver = 0.80,
  quantile_trait = 0.00,
  
  min_abs_loading_rum   = 0,
  min_abs_loading_ile   = 0,
  min_abs_loading_col   = 0,
  min_abs_loading_liver = 0,
  min_abs_loading_trait = 0,
  
  # ---------- 第三层：链条结构是什么 ----------
  edge_source_liver = 0.40,
  edge_source_liver_candidate = 0.30,
  edge_liver_trait = 0.40,
  edge_source_trait = 0.00,
  
  liver_hub_min_links = 2,
  max_chain_per_source = 2,
  max_total_chains = 10
)

# ============================================================
# 三、辅助函数
# ============================================================
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  suppressMessages(readr::read_csv(fp, show_col_types = FALSE))
}

safe_read_lines <- function(fp) {
  if (!file.exists(fp)) return(character(0))
  readLines(fp, warn = FALSE, encoding = "UTF-8")
}

extract_first_number <- function(text, pattern) {
  hit <- stringr::str_match(text, pattern)[, 2]
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NA_real_)
  as.numeric(hit[1])
}

parse_model_diagnostics <- function(txt_lines) {
  txt <- paste(txt_lines, collapse = "\n")
  trait <- stringr::str_match(txt, 'Trait:\\s*\\n\\[1\\]\\s*"([^"]+)"')[, 2]
  if (is.na(trait)) {
    trait <- stringr::str_match(txt, 'Trait:\\s*([A-Za-z0-9_\\.]+)')[, 2]
  }
  
  ave_inner <- extract_first_number(txt, '\\$AVE_inner\\s*\\[1\\]\\s*([0-9eE\\.\\-]+)')
  ave_outer <- extract_first_number(txt, '\\$AVE_outer\\s*\\[1\\]\\s*([0-9eE\\.\\-]+)')
  
  tibble(
    Trait = ifelse(is.na(trait), NA_character_, trait),
    AVE_outer = ave_outer,
    AVE_inner = ave_inner
  )
}

prep_corr_long <- function(corr_df) {
  stopifnot(!is.null(corr_df))
  first_col <- names(corr_df)[1]
  
  corr_df %>%
    rename(row_block = all_of(first_col)) %>%
    pivot_longer(-row_block, names_to = "col_block", values_to = "correlation") %>%
    mutate(
      block_a = pmin(row_block, col_block),
      block_b = pmax(row_block, col_block)
    ) %>%
    filter(block_a != block_b) %>%
    distinct(block_a, block_b, .keep_all = TRUE) %>%
    mutate(abs_corr = abs(correlation))
}

judge_model <- function(diag_df, corr_long, params) {
  ave_inner <- diag_df$AVE_inner[1]
  
  trait_liver_ok <- corr_long %>%
    filter((block_a == "Liver_comp1" & block_b == "Trait_comp1") |
             (block_a == "Trait_comp1" & block_b == "Liver_comp1")) %>%
    summarise(
      ok = any(abs_corr >= params$corr_trait_liver),
      corr = dplyr::first(correlation)
    )
  
  other_ok_n <- corr_long %>%
    filter(!(block_a == "Liver_comp1" & block_b == "Trait_comp1")) %>%
    summarise(n = sum(abs_corr >= params$corr_other_pairs)) %>%
    pull(n)
  
  model_level <- dplyr::case_when(
    !is.na(ave_inner) && ave_inner >= params$ave_inner_a &&
      isTRUE(trait_liver_ok$ok) && other_ok_n >= params$min_other_pairs_n ~ "A_模型成立_可作为主体解释",
    !is.na(ave_inner) && ave_inner >= params$ave_inner_b ~ "B_模型基本成立_可辅助解释",
    TRUE ~ "C_模型较弱_不建议深讲"
  )
  
  tibble(
    trait_liver_corr = trait_liver_ok$corr,
    trait_liver_pass = isTRUE(trait_liver_ok$ok),
    other_pair_pass_n = other_ok_n,
    model_level = model_level,
    model_pass = model_level != "C_模型较弱_不建议深讲"
  )
}

pick_core_features_one_block <- function(df_block, block_name, params) {
  if (nrow(df_block) == 0) return(df_block)
  
  min_abs_loading <- switch(
    block_name,
    Rum   = params$min_abs_loading_rum,
    Ile   = params$min_abs_loading_ile,
    Col   = params$min_abs_loading_col,
    Liver = params$min_abs_loading_liver,
    Trait = params$min_abs_loading_trait,
    0
  )
  
  out <- df_block %>%
    filter(abs_comp1 >= min_abs_loading)
  
  if (params$feature_method == "top_n") {
    n_keep <- switch(
      block_name,
      Rum   = params$top_n_rum,
      Ile   = params$top_n_ile,
      Col   = params$top_n_col,
      Liver = params$top_n_liver,
      Trait = params$top_n_trait,
      10
    )
    n_keep2 <- min(n_keep, nrow(out))
    out <- out %>%
      arrange(desc(abs_comp1)) %>%
      slice_head(n = n_keep2)
    
  } else if (params$feature_method == "quantile") {
    q_cut <- switch(
      block_name,
      Rum   = params$quantile_rum,
      Ile   = params$quantile_ile,
      Col   = params$quantile_col,
      Liver = params$quantile_liver,
      Trait = params$quantile_trait,
      0.8
    )
    thr <- stats::quantile(out$abs_comp1, probs = q_cut, na.rm = TRUE)
    out <- out %>% filter(abs_comp1 >= thr)
  }
  
  out
}

build_core_features <- function(load_df, trait_name, params) {
  stopifnot(!is.null(load_df))
  
  load_df <- load_df %>%
    mutate(
      direction_comp1 = if_else(comp1 >= 0, "positive", "negative"),
      sign_comp1 = if_else(comp1 >= 0, 1L, -1L)
    )
  
  trait_row <- load_df %>% filter(Block == "Trait") %>% arrange(desc(abs_comp1)) %>% slice(1)
  trait_sign <- if (nrow(trait_row) == 0) 1L else trait_row$sign_comp1[1]
  trait_feature <- if (nrow(trait_row) == 0) trait_name else trait_row$Feature[1]
  
  core_df <- load_df %>%
    group_split(Block) %>%
    map_dfr(~ pick_core_features_one_block(.x, unique(.x$Block), params)) %>%
    mutate(
      Trait = trait_name,
      trait_feature = trait_feature,
      trait_sign_comp1 = trait_sign,
      module_vs_trait = if_else(
        sign_comp1 == trait_sign,
        "same_direction_as_trait",
        "opposite_direction_to_trait"
      ),
      core_feature = TRUE
    ) %>%
    arrange(factor(Block, levels = c("Rum", "Ile", "Col", "Liver", "Trait")), desc(abs_comp1))
  
  list(core_df = core_df, trait_sign = trait_sign, trait_feature = trait_feature)
}

# ------------------------------------------------------------
# 新增：block 内 z-score 标准化函数
# ------------------------------------------------------------
standardize_chain_score_within_block <- function(df) {
  if (nrow(df) == 0) return(df)
  
  df %>%
    group_by(source_block) %>%
    mutate(
      block_chain_n = n(),
      block_chain_score_mean = mean(chain_score_raw, na.rm = TRUE),
      block_chain_score_sd   = stats::sd(chain_score_raw, na.rm = TRUE),
      chain_score_z = dplyr::case_when(
        is.na(block_chain_score_sd) ~ 0,
        block_chain_score_sd == 0 ~ 0,
        TRUE ~ (chain_score_raw - block_chain_score_mean) / block_chain_score_sd
      )
    ) %>%
    ungroup()
}

build_chain_table <- function(edges_df, core_features_df, trait_name, params) {
  stopifnot(!is.null(edges_df))
  stopifnot(!is.null(core_features_df))
  
  source_core <- core_features_df %>%
    filter(Block %in% c("Rum", "Ile", "Col")) %>%
    select(
      Block, Feature,
      source_loading = comp1,
      source_abs_loading = abs_comp1,
      source_direction = direction_comp1
    )
  
  liver_core <- core_features_df %>%
    filter(Block == "Liver") %>%
    select(
      Feature,
      liver_loading = comp1,
      liver_abs_loading = abs_comp1,
      liver_direction = direction_comp1
    )
  
  trait_core <- core_features_df %>%
    filter(Block == "Trait") %>%
    slice(1) %>%
    transmute(
      trait_feature = Feature,
      trait_loading = comp1,
      trait_abs_loading = abs_comp1,
      trait_direction = direction_comp1
    )
  
  source_liver_edges <- edges_df %>%
    filter(pair %in% c("Rum_Liver", "Ile_Liver", "Col_Liver")) %>%
    rename(
      source_block = block1, source_feature = var1,
      liver_block = block2, liver_feature = var2,
      corr_source_liver = correlation,
      abs_corr_source_liver = abs_corr
    ) %>%
    filter(abs_corr_source_liver >= params$edge_source_liver)
  
  liver_trait_edges <- edges_df %>%
    filter(pair == "Liver_Trait") %>%
    rename(
      liver_block = block1, liver_feature = var1,
      trait_block = block2, trait_feature = var2,
      corr_liver_trait = correlation,
      abs_corr_liver_trait = abs_corr
    ) %>%
    filter(abs_corr_liver_trait >= params$edge_liver_trait)
  
  source_trait_edges <- edges_df %>%
    filter(pair %in% c("Rum_Trait", "Ile_Trait", "Col_Trait")) %>%
    rename(
      source_block = block1, source_feature = var1,
      trait_block = block2, trait_feature = var2,
      corr_source_trait = correlation,
      abs_corr_source_trait = abs_corr
    )
  
  if (params$edge_source_trait > 0) {
    source_trait_edges <- source_trait_edges %>%
      filter(abs_corr_source_trait >= params$edge_source_trait)
  }
  
  chain_df <- source_liver_edges %>%
    inner_join(source_core, by = c("source_block" = "Block", "source_feature" = "Feature")) %>%
    inner_join(liver_core, by = c("liver_feature" = "Feature")) %>%
    inner_join(
      liver_trait_edges %>% select(liver_feature, trait_feature, corr_liver_trait, abs_corr_liver_trait),
      by = "liver_feature"
    ) %>%
    left_join(
      source_trait_edges %>% select(source_block, source_feature, trait_feature,
                                    corr_source_trait, abs_corr_source_trait),
      by = c("source_block", "source_feature", "trait_feature")
    ) %>%
    mutate(
      Trait = trait_name,
      source_to_liver_direction = if_else(corr_source_liver >= 0, "positive", "negative"),
      liver_to_trait_direction  = if_else(corr_liver_trait >= 0, "positive", "negative"),
      source_to_trait_direction = case_when(
        is.na(corr_source_trait) ~ NA_character_,
        corr_source_trait >= 0 ~ "positive",
        TRUE ~ "negative"
      ),
      chain_text = paste0(
        source_block, ":", source_feature,
        if_else(corr_source_liver >= 0, " ↑→ ", " ↓→ "),
        "Liver:", liver_feature,
        if_else(corr_liver_trait >= 0, " ↑→ ", " ↓→ "),
        trait_feature
      ),
      chain_score_raw = source_abs_loading * liver_abs_loading *
        abs_corr_source_liver * abs_corr_liver_trait,
      chain_type = case_when(
        corr_source_liver >= 0 & corr_liver_trait >= 0 ~ "同向-同向链",
        corr_source_liver >= 0 & corr_liver_trait <  0 ~ "同向-反向链",
        corr_source_liver <  0 & corr_liver_trait >= 0 ~ "反向-同向链",
        TRUE ~ "反向-反向链"
      )
    ) %>%
    standardize_chain_score_within_block() %>%
    group_by(source_block, source_feature) %>%
    arrange(desc(chain_score_z), desc(chain_score_raw), .by_group = TRUE) %>%
    slice_head(n = params$max_chain_per_source) %>%
    ungroup() %>%
    arrange(desc(chain_score_z), desc(chain_score_raw)) %>%
    slice_head(n = params$max_total_chains)
  
  candidate_edges <- edges_df %>%
    filter(pair %in% c("Rum_Liver", "Ile_Liver", "Col_Liver")) %>%
    rename(
      source_block = block1, source_feature = var1,
      liver_feature = var2,
      corr_source_liver = correlation,
      abs_corr_source_liver = abs_corr
    ) %>%
    filter(
      abs_corr_source_liver >= params$edge_source_liver_candidate,
      abs_corr_source_liver < params$edge_source_liver
    ) %>%
    inner_join(source_core, by = c("source_block" = "Block", "source_feature" = "Feature")) %>%
    inner_join(liver_core, by = c("liver_feature" = "Feature")) %>%
    mutate(
      Trait = trait_name,
      candidate_edge_score_raw = source_abs_loading * liver_abs_loading * abs_corr_source_liver
    ) %>%
    group_by(source_block) %>%
    mutate(
      candidate_block_n = n(),
      candidate_edge_score_mean = mean(candidate_edge_score_raw, na.rm = TRUE),
      candidate_edge_score_sd = stats::sd(candidate_edge_score_raw, na.rm = TRUE),
      candidate_edge_score_z = dplyr::case_when(
        is.na(candidate_edge_score_sd) ~ 0,
        candidate_edge_score_sd == 0 ~ 0,
        TRUE ~ (candidate_edge_score_raw - candidate_edge_score_mean) / candidate_edge_score_sd
      )
    ) %>%
    ungroup() %>%
    arrange(desc(candidate_edge_score_z), desc(candidate_edge_score_raw))
  
  list(chain_df = chain_df, candidate_edges = candidate_edges)
}

build_liver_hub_table <- function(chain_df, params) {
  if (is.null(chain_df) || nrow(chain_df) == 0) {
    return(tibble())
  }
  
  trait_name <- unique(chain_df$Trait)
  if (length(trait_name) == 0) trait_name <- NA_character_
  
  chain_df %>%
    group_by(liver_feature) %>%
    summarise(
      n_source_links = n_distinct(paste(source_block, source_feature, sep = "::")),
      n_rum_links = n_distinct(source_feature[source_block == "Rum"]),
      n_ile_links = n_distinct(source_feature[source_block == "Ile"]),
      n_col_links = n_distinct(source_feature[source_block == "Col"]),
      mean_abs_corr_source_liver = mean(abs_corr_source_liver, na.rm = TRUE),
      abs_corr_liver_trait = mean(abs_corr_liver_trait, na.rm = TRUE),
      mean_chain_score_raw = mean(chain_score_raw, na.rm = TRUE),
      mean_chain_score_z = mean(chain_score_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Trait = trait_name[1],
      hub_flag = n_source_links >= params$liver_hub_min_links
    ) %>%
    select(Trait, everything()) %>%
    arrange(desc(hub_flag), desc(n_source_links), desc(mean_chain_score_z), desc(mean_chain_score_raw))
}

build_module_summary <- function(core_features_df, chain_df, trait_name) {
  feature_module <- core_features_df %>%
    filter(Block %in% c("Rum", "Ile", "Col", "Liver")) %>%
    group_by(Block, module_vs_trait) %>%
    summarise(
      n_feature = n(),
      features = paste(Feature, collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(Trait = trait_name) %>%
    select(Trait, everything())
  
  chain_module <- if (!is.null(chain_df) && nrow(chain_df) > 0) {
    chain_df %>%
      group_by(chain_type) %>%
      summarise(
        n_chain = n(),
        chains = paste(chain_text, collapse = " | "),
        .groups = "drop"
      ) %>%
      mutate(
        Trait = trait_name,
        mean_chain_score_raw = mean(chain_df$chain_score_raw, na.rm = TRUE),
        mean_chain_score_z = mean(chain_df$chain_score_z, na.rm = TRUE)
      ) %>%
      select(Trait, everything())
  } else {
    tibble(
      Trait = trait_name,
      chain_type = character(),
      n_chain = integer(),
      chains = character(),
      mean_chain_score_raw = numeric(),
      mean_chain_score_z = numeric()
    )
  }
  
  list(feature_module = feature_module, chain_module = chain_module)
}

add_sheet_with_style <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  hs <- createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom",
    halign = "center"
  )
  addStyle(wb, sheet_name, hs, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}

process_one_trait <- function(trait_dir, output_root, params) {
  trait_folder <- basename(trait_dir)
  message("\n==============================")
  message("Processing: ", trait_folder)
  message("==============================")
  
  out_dir <- file.path(output_root, trait_folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fp_diag   <- file.path(trait_dir, "model_diagnostics.txt")
  fp_corr   <- file.path(trait_dir, "block_score_correlations.csv")
  fp_load   <- file.path(trait_dir, "loadings_nonzero_all_blocks.csv")
  fp_edges  <- file.path(trait_dir, "network_edges.csv")
  fp_scores <- file.path(trait_dir, "component1_scores_wide.csv")
  fp_dim    <- file.path(trait_dir, "block_dimensions.csv")
  
  diag_df  <- parse_model_diagnostics(safe_read_lines(fp_diag))
  corr_df  <- safe_read_csv(fp_corr)
  load_df  <- safe_read_csv(fp_load)
  edges_df <- safe_read_csv(fp_edges)
  score_df <- safe_read_csv(fp_scores)
  dim_df   <- safe_read_csv(fp_dim)
  
  if (is.null(corr_df) || is.null(load_df) || is.null(edges_df)) {
    warning("缺少关键文件，跳过：", trait_folder)
    return(NULL)
  }
  
  trait_name <- diag_df$Trait[1]
  if (is.na(trait_name) || is.null(trait_name) || trait_name == "") {
    tmp_trait <- load_df %>% filter(Block == "Trait") %>% slice(1) %>% pull(Feature)
    trait_name <- ifelse(length(tmp_trait) == 0, stringr::str_remove(trait_folder, "^sGCCA_"), tmp_trait[1])
    diag_df$Trait <- trait_name
  }
  
  corr_long   <- prep_corr_long(corr_df)
  model_judge <- judge_model(diag_df, corr_long, params) %>% mutate(Trait = trait_name)
  layer1_df   <- diag_df %>% bind_cols(model_judge %>% select(-Trait))
  
  core_obj <- build_core_features(load_df, trait_name, params)
  core_df  <- core_obj$core_df
  
  chain_obj <- build_chain_table(edges_df, core_df, trait_name, params)
  chain_df  <- chain_obj$chain_df
  cand_df   <- chain_obj$candidate_edges
  
  liver_hub_df <- build_liver_hub_table(chain_df, params)
  module_obj   <- build_module_summary(core_df, chain_df, trait_name)
  
  layer1_pairs_df <- corr_long %>%
    mutate(Trait = trait_name) %>%
    select(Trait, block_a, block_b, correlation, abs_corr)
  
  write_csv(layer1_df, file.path(out_dir, "01_model_assessment.csv"))
  write_csv(layer1_pairs_df, file.path(out_dir, "01b_block_score_correlations_long.csv"))
  if (!is.null(dim_df))   write_csv(dim_df,   file.path(out_dir, "01c_block_dimensions.csv"))
  if (!is.null(score_df)) write_csv(score_df, file.path(out_dir, "01d_component1_scores_wide.csv"))
  write_csv(core_df,      file.path(out_dir, "02_core_features.csv"))
  write_csv(chain_df,     file.path(out_dir, "03_core_chains.csv"))
  write_csv(cand_df,      file.path(out_dir, "03b_candidate_source_liver_edges.csv"))
  write_csv(liver_hub_df, file.path(out_dir, "04_liver_hubs.csv"))
  write_csv(module_obj$feature_module, file.path(out_dir, "05_feature_module_summary.csv"))
  write_csv(module_obj$chain_module,   file.path(out_dir, "05b_chain_module_summary.csv"))
  
  wb <- createWorkbook()
  add_sheet_with_style(wb, "01_model_assessment", layer1_df)
  add_sheet_with_style(wb, "01b_block_corr_long", layer1_pairs_df)
  if (!is.null(dim_df))   add_sheet_with_style(wb, "01c_block_dimensions", dim_df)
  if (!is.null(score_df)) add_sheet_with_style(wb, "01d_component1_scores", score_df)
  add_sheet_with_style(wb, "02_core_features", core_df)
  add_sheet_with_style(wb, "03_core_chains", chain_df)
  add_sheet_with_style(wb, "03b_candidate_edges", cand_df)
  add_sheet_with_style(wb, "04_liver_hubs", liver_hub_df)
  add_sheet_with_style(wb, "05_feature_modules", module_obj$feature_module)
  add_sheet_with_style(wb, "05b_chain_modules", module_obj$chain_module)
  saveWorkbook(wb, file.path(out_dir, paste0(trait_folder, "_postprocessed.xlsx")), overwrite = TRUE)
  
  list(
    summary = layer1_df %>%
      select(Trait, AVE_outer, AVE_inner, trait_liver_corr,
             trait_liver_pass, other_pair_pass_n, model_level, model_pass),
    core_features = core_df,
    core_chains = chain_df,
    liver_hubs = liver_hub_df,
    feature_module = module_obj$feature_module,
    chain_module = module_obj$chain_module
  )
}

# ============================================================
# 四、批量处理八个表型
# ============================================================
trait_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
trait_dirs <- trait_dirs[grepl("^sGCCA_", basename(trait_dirs))]

if (length(trait_dirs) == 0) {
  stop("在 input_root 下未找到 sGCCA_* 文件夹，请检查路径。")
}

res_list <- purrr::map(trait_dirs, process_one_trait, output_root = output_root, params = params)
res_list <- res_list[!vapply(res_list, is.null, logical(1))]

# ============================================================
# 五、跨表型总汇总
# ============================================================
all_summary        <- bind_rows(purrr::map(res_list, "summary"))
all_core_features  <- bind_rows(purrr::map(res_list, "core_features"))
all_core_chains    <- bind_rows(purrr::map(res_list, "core_chains"))
all_liver_hubs     <- bind_rows(purrr::map(res_list, "liver_hubs"))
all_feature_module <- bind_rows(purrr::map(res_list, "feature_module"))
all_chain_module   <- bind_rows(purrr::map(res_list, "chain_module"))

write_csv(all_summary,        file.path(output_root, "sGCCA_postprocess_summary_all_traits.csv"))
write_csv(all_core_features,  file.path(output_root, "sGCCA_core_features_all_traits.csv"))
write_csv(all_core_chains,    file.path(output_root, "sGCCA_core_chains_all_traits.csv"))
write_csv(all_liver_hubs,     file.path(output_root, "sGCCA_liver_hubs_all_traits.csv"))
write_csv(all_feature_module, file.path(output_root, "sGCCA_feature_modules_all_traits.csv"))
write_csv(all_chain_module,   file.path(output_root, "sGCCA_chain_modules_all_traits.csv"))

wb_all <- createWorkbook()
add_sheet_with_style(wb_all, "summary_all_traits", all_summary)
add_sheet_with_style(wb_all, "core_features_all", all_core_features)
add_sheet_with_style(wb_all, "core_chains_all", all_core_chains)
add_sheet_with_style(wb_all, "liver_hubs_all", all_liver_hubs)
add_sheet_with_style(wb_all, "feature_modules_all", all_feature_module)
add_sheet_with_style(wb_all, "chain_modules_all", all_chain_module)
saveWorkbook(wb_all, file.path(output_root, "sGCCA_postprocess_all_traits.xlsx"), overwrite = TRUE)

message("\n全部完成。输出目录：", output_root)