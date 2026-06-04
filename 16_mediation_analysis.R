############################
## 0. Environment and packages
############################
rm(list = ls())
gc()

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "zCompositions", "compositions", "mediation",
  "openxlsx", "tibble", "parallel", "speedglm"
)

to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(zCompositions)
  library(compositions)
  library(mediation)
  library(openxlsx)
  library(tibble)
  library(parallel)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(123)

############################
## 1. Paths
############################
# Input and output directories
# Please place the required input files in "data/mediation/input".
# Output files will be saved in "results/mediation".

base_dir <- file.path("data", "mediation", "input")
out_dir  <- file.path("results", "mediation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  col          = file.path(base_dir, "Colon_genus_abundance.xlsx"),
  liver        = file.path(base_dir, "肝脏筛选基因.xlsx"),
  blood        = file.path(base_dir, "blood_phenotype_sGCCA_input.xlsx"),
  tail         = file.path(base_dir, "尾脂数据_16s编号转换.xlsx"),
  old_combined = file.path(base_dir, "combined_mediation_results.csv")
)

############################
## 2. Sheet configuration
############################
sheet_cfg <- list(
  col_sheet   = 1,
  liver_sheet = "formal_log2TPM",
  blood_sheet = "selected_raw_table",
  tail_sheet  = "尾脂_16s编号"
)

############################
## 3. Parameters
############################
prevalence_cut <- 0.50
min_nonzero_n  <- NULL

boot_sims <- 1000
save_every_n <- 5000

n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
chunk_size <- 5000

############################
## 4. File-name configuration
############################
checkpoint_files <- list(
  col_partial_csv = file.path(out_dir, "col_results_partial.csv"),
  col_partial_rds = file.path(out_dir, "col_results_partial.rds"),
  
  col_final_csv   = file.path(out_dir, "col_mediation_results.csv"),
  col_final_rds   = file.path(out_dir, "col_mediation_results.rds"),
  
  combined_csv    = file.path(out_dir, "combined_mediation_results.csv"),
  combined_rds    = file.path(out_dir, "combined_mediation_results.rds"),
  
  summary_rds     = file.path(out_dir, "analysis_summary.rds"),
  excel_out       = file.path(out_dir, "mediation_full_analysis_results.xlsx")
)

############################
## 5. Helper functions
############################

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x
}

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  suppressWarnings(as.numeric(x))
}

zscore_df <- function(df, id_col = "SampleID") {
  stopifnot(id_col %in% colnames(df))
  
  out <- as.data.frame(df, check.names = FALSE)
  num_cols <- setdiff(colnames(out), id_col)
  
  for (cc in num_cols) {
    x <- safe_numeric(out[[cc]])
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      out[[cc]] <- 0
    } else {
      out[[cc]] <- (x - m) / s
    }
  }
  out
}

make_safe_names_with_map <- function(x, prefix = "V") {
  orig <- as.character(x)
  safe <- make.names(orig, unique = TRUE)
  tibble::tibble(original = orig, safe = safe, prefix = prefix)
}

## 5.1 Read the genus abundance table: the first column is Genus and the remaining columns are samples
read_genus_table <- function(fp, sheet = 1) {
  df <- readxl::read_excel(fp, sheet = sheet)
  df <- as.data.frame(df, check.names = FALSE)
  
  if (!"Genus" %in% colnames(df)) {
    stop("The 'Genus' column was not found in the file: ", fp)
  }
  
  df <- df[!is.na(df$Genus) & df$Genus != "", , drop = FALSE]
  
  sample_cols <- setdiff(colnames(df), "Genus")
  sample_cols_clean <- clean_id(sample_cols)
  colnames(df) <- c("Genus", sample_cols_clean)
  
  for (cc in sample_cols_clean) {
    df[[cc]] <- safe_numeric(df[[cc]])
  }
  df[is.na(df)] <- 0
  
  ## Merge duplicated genus names
  df <- df %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.matrix(df[, sample_cols_clean, drop = FALSE])
  rownames(mat) <- df$Genus
  mode(mat) <- "numeric"
  
  mat
}

## 5.2 Prevalence filtering based on the proportion of samples with abundance > 0
prevalence_filter <- function(mat, prevalence_cut = 0.50, min_nonzero_n = NULL) {
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

## 5.3 CZM zero replacement and CLR transformation (sample x genus)
clr_transform_czm <- function(mat_genus_by_sample) {
  if (nrow(mat_genus_by_sample) < 2) {
    stop("The number of features after filtering is < 2; CLR transformation cannot be performed. Please relax prevalence_cut.")
  }
  
  x <- t(mat_genus_by_sample)  ## rows = samples, columns = genera
  
  x_nozero <- zCompositions::cmultRepl(
    X = x,
    method = "CZM",
    output = "p-counts",
    label = 0
  )
  
  x_clr <- compositions::clr(x_nozero, base = exp(1))
  x_clr <- as.matrix(x_clr)
  
  x_clr
}

## 5.4 Prepare the microbial block
prep_microbe_block_pairwise_style <- function(filepath, sheet = 1, block_name = "Col",
                                              prevalence_cut = 0.50, min_nonzero_n = NULL) {
  message("Reading microbial data: ", block_name)
  
  raw_mat <- read_genus_table(filepath, sheet = sheet)
  
  prev_obj <- prevalence_filter(
    raw_mat,
    prevalence_cut = prevalence_cut,
    min_nonzero_n = min_nonzero_n
  )
  
  message(block_name, " features retained after prevalence filtering: ", nrow(prev_obj$mat))
  
  clr_mat <- clr_transform_czm(prev_obj$mat)
  
  clr_z_df <- zscore_df(
    tibble::rownames_to_column(as.data.frame(clr_mat, check.names = FALSE), "SampleID"),
    id_col = "SampleID"
  )
  clr_z_df$SampleID <- clean_id(clr_z_df$SampleID)
  
  feat_map <- make_safe_names_with_map(colnames(clr_z_df)[-1], prefix = block_name)
  safe_names <- feat_map$safe
  names(safe_names) <- feat_map$original
  colnames(clr_z_df) <- c("SampleID", safe_names[colnames(clr_z_df)[-1]])
  
  list(
    data = clr_z_df,
    map  = feat_map,
    prevalence_info = prev_obj$info,
    n_feature_before = nrow(raw_mat),
    n_feature_after  = ncol(clr_z_df) - 1,
    block = block_name
  )
}

## 5.5 liver expression
prep_liver_expression <- function(filepath, sheet = "formal_log2TPM") {
  message("Reading the liver gene expression matrix...")
  
  raw_df <- readxl::read_excel(filepath, sheet = sheet, skip = 1)
  raw_df <- as.data.frame(raw_df, check.names = FALSE)
  
  if (!"Gene" %in% colnames(raw_df)) {
    stop("The 'Gene' column was not found in the liver expression table.")
  }
  
  raw_df <- raw_df[!is.na(raw_df$Gene) & raw_df$Gene != "", , drop = FALSE]
  
  sample_cols <- setdiff(colnames(raw_df), "Gene")
  sample_cols_clean <- clean_id(sample_cols)
  colnames(raw_df) <- c("Gene", sample_cols_clean)
  
  for (cc in sample_cols_clean) {
    raw_df[[cc]] <- safe_numeric(raw_df[[cc]])
  }
  
  raw_df <- raw_df %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.matrix(raw_df[, sample_cols_clean, drop = FALSE])
  rownames(mat) <- raw_df$Gene
  mode(mat) <- "numeric"
  
  sf <- as.data.frame(t(mat), check.names = FALSE)
  sf <- tibble::rownames_to_column(sf, var = "SampleID")
  sf$SampleID <- clean_id(sf$SampleID)
  
  sf_z <- zscore_df(sf, id_col = "SampleID")
  
  gene_map <- make_safe_names_with_map(colnames(sf_z)[-1], prefix = "Liver")
  safe_names <- gene_map$safe
  names(safe_names) <- gene_map$original
  colnames(sf_z) <- c("SampleID", safe_names[colnames(sf_z)[-1]])
  
  list(
    data = sf_z,
    map  = gene_map,
    n_gene_after = ncol(sf_z) - 1
  )
}

## 5.6 blood traits
prep_blood_traits <- function(filepath, sheet = "selected_raw_table") {
  message("Reading blood phenotypes (5 traits)...")
  
  raw_df <- readxl::read_excel(filepath, sheet = sheet)
  raw_df <- as.data.frame(raw_df, check.names = FALSE)
  
  need_cols <- c("Sample16S", "TC", "TG", "LDL", "TBA", "GLU")
  miss_cols <- setdiff(need_cols, colnames(raw_df))
  if (length(miss_cols) > 0) {
    stop("The blood phenotype table is missing columns: ", paste(miss_cols, collapse = ", "))
  }
  
  out <- raw_df %>%
    dplyr::transmute(
      SampleID = clean_id(Sample16S),
      TC  = safe_numeric(TC),
      TG  = safe_numeric(TG),
      LDL = safe_numeric(LDL),
      TBA = safe_numeric(TBA),
      GLU = safe_numeric(GLU)
    ) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out
}

## 5.7 tail traits
prep_tail_traits <- function(filepath, sheet = "尾脂_16s编号") {
  message("Reading tail-fat phenotypes (3 traits)...")
  
  raw_df <- readxl::read_excel(filepath, sheet = sheet, skip = 4)
  raw_df <- as.data.frame(raw_df, check.names = FALSE)
  
  need_cols <- c("瘤胃和肠道16s编号", "尾脂g", "尾脂/胴体重g/kg", "尾脂/宰前活重g/kg")
  miss_cols <- setdiff(need_cols, colnames(raw_df))
  if (length(miss_cols) > 0) {
    stop("The tail-fat phenotype table is missing columns: ", paste(miss_cols, collapse = ", "))
  }
  
  out <- raw_df %>%
    dplyr::transmute(
      SampleID = clean_id(`瘤胃和肠道16s编号`),
      TailFat_g = safe_numeric(`尾脂g`),
      TailFat_Carcass_gkg = safe_numeric(`尾脂/胴体重g/kg`),
      TailFat_PreBW_gkg = safe_numeric(`尾脂/宰前活重g/kg`)
    ) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out
}

## 5.8 Merge traits and apply z-score standardization
prep_all_traits <- function(blood_df, tail_df) {
  trait_df <- dplyr::full_join(blood_df, tail_df, by = "SampleID")
  trait_z <- zscore_df(trait_df, id_col = "SampleID")
  
  trait_map <- tibble::tibble(
    original = c("TC", "TG", "LDL", "TBA", "GLU",
                 "TailFat_g", "TailFat_Carcass_gkg", "TailFat_PreBW_gkg"),
    safe = make.names(c("TC", "TG", "LDL", "TBA", "GLU",
                        "TailFat_g", "TailFat_Carcass_gkg", "TailFat_PreBW_gkg"), unique = TRUE),
    prefix = "Trait"
  )
  
  safe_names <- trait_map$safe
  names(safe_names) <- trait_map$original
  colnames(trait_z) <- c("SampleID", safe_names[colnames(trait_z)[-1]])
  
  list(
    data = trait_z,
    map  = trait_map
  )
}

## 5.9 Align the three data blocks
align_three_blocks <- function(source_df, mediator_df, trait_df) {
  common_ids <- Reduce(intersect, list(
    source_df$SampleID,
    mediator_df$SampleID,
    trait_df$SampleID
  ))
  
  common_ids <- sort(unique(common_ids))
  
  if (length(common_ids) < 5) {
    stop("Too few shared samples across the three data blocks (<5). Please check whether sample IDs are consistent.")
  }
  
  source_aln <- source_df %>%
    dplyr::filter(SampleID %in% common_ids) %>%
    dplyr::arrange(match(SampleID, common_ids))
  
  mediator_aln <- mediator_df %>%
    dplyr::filter(SampleID %in% common_ids) %>%
    dplyr::arrange(match(SampleID, common_ids))
  
  trait_aln <- trait_df %>%
    dplyr::filter(SampleID %in% common_ids) %>%
    dplyr::arrange(match(SampleID, common_ids))
  
  stopifnot(identical(source_aln$SampleID, mediator_aln$SampleID))
  stopifnot(identical(source_aln$SampleID, trait_aln$SampleID))
  
  list(
    source   = source_aln,
    mediator = mediator_aln,
    trait    = trait_aln,
    n_sample = length(common_ids),
    ids      = common_ids
  )
}

## 5.10 Build all candidate chains without pre-screening
build_all_candidates <- function(source_df, mediator_df, trait_df,
                                 source_map, mediator_map, trait_map,
                                 source_block = "Col") {
  message("Building all candidate chains for: ", source_block)
  
  source_names   <- setdiff(colnames(source_df), "SampleID")
  mediator_names <- setdiff(colnames(mediator_df), "SampleID")
  trait_names    <- setdiff(colnames(trait_df), "SampleID")
  
  cand_df <- tidyr::expand_grid(
    source_safe = source_names,
    mediator_safe = mediator_names,
    trait_safe = trait_names
  ) %>%
    dplyr::mutate(source_block = source_block)
  
  source_map2 <- source_map %>%
    dplyr::select(source_safe = safe, source_feature = original)
  mediator_map2 <- mediator_map %>%
    dplyr::select(mediator_safe = safe, mediator_feature = original)
  trait_map2 <- trait_map %>%
    dplyr::select(trait_safe = safe, trait_feature = original)
  
  cand_df <- cand_df %>%
    dplyr::left_join(source_map2, by = "source_safe") %>%
    dplyr::left_join(mediator_map2, by = "mediator_safe") %>%
    dplyr::left_join(trait_map2, by = "trait_safe") %>%
    dplyr::distinct(source_block, source_safe, mediator_safe, trait_safe, .keep_all = TRUE)
  
  cand_df
}

extract_lm_coef <- function(fit, term) {
  sm <- summary(fit)
  cf <- sm$coefficients
  if (!term %in% rownames(cf)) {
    return(c(est = NA_real_, p = NA_real_))
  }
  c(est = cf[term, "Estimate"], p = cf[term, "Pr(>|t|)"])
}

get_med_item <- function(obj, candidates) {
  for (nm in candidates) {
    if (!is.null(obj[[nm]]) && length(obj[[nm]]) > 0) return(obj[[nm]][1])
  }
  NA_real_
}

make_chain_id <- function(df) {
  paste(df$source_block, df$source_safe, df$mediator_safe, df$trait_safe, sep = "||")
}

empty_med_row <- function(cid, source_block, source_safe, mediator_safe, trait_safe, n_use, status_txt) {
  tibble::tibble(
    chain_id = cid,
    source_block = source_block,
    source_safe = source_safe,
    mediator_safe = mediator_safe,
    trait_safe = trait_safe,
    n = n_use,
    
    a_est = NA_real_,
    a_p   = NA_real_,
    b_est = NA_real_,
    b_p   = NA_real_,
    cprime_est = NA_real_,
    cprime_p   = NA_real_,
    total_est  = NA_real_,
    total_p    = NA_real_,
    
    acme = NA_real_,
    acme_p = NA_real_,
    ade = NA_real_,
    ade_p = NA_real_,
    prop_med = NA_real_,
    prop_med_p = NA_real_,
    med_total = NA_real_,
    med_total_p = NA_real_,
    
    status = status_txt
  )
}

save_partial_results <- function(res_df, csv_path, rds_path) {
  utils::write.csv(res_df, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(res_df, rds_path)
}

run_one_mediation_chain <- function(i, cand_df, source_mat, mediator_mat, trait_mat, sims = 1000) {
  s_idx <- cand_df$source_idx[i]
  m_idx <- cand_df$mediator_idx[i]
  t_idx <- cand_df$trait_idx[i]
  cid   <- cand_df$chain_id[i]
  
  source_vec <- source_mat[, s_idx]
  medi_vec   <- mediator_mat[, m_idx]
  out_vec    <- trait_mat[, t_idx]
  ok <- is.finite(source_vec) & is.finite(medi_vec) & is.finite(out_vec)
  n_use <- sum(ok)
  
  if (n_use < 5) {
    return(empty_med_row(
      cid = cid,
      source_block = cand_df$source_block[i],
      source_safe = cand_df$source_safe[i],
      mediator_safe = cand_df$mediator_safe[i],
      trait_safe = cand_df$trait_safe[i],
      n_use = n_use,
      status_txt = "skip_n_too_small"
    ))
  }
  
  source_vec <- source_vec[ok]
  medi_vec   <- medi_vec[ok]
  out_vec    <- out_vec[ok]
  
  if (stats::sd(source_vec) == 0 || stats::sd(medi_vec) == 0 || stats::sd(out_vec) == 0) {
    return(empty_med_row(
      cid = cid,
      source_block = cand_df$source_block[i],
      source_safe = cand_df$source_safe[i],
      mediator_safe = cand_df$mediator_safe[i],
      trait_safe = cand_df$trait_safe[i],
      n_use = n_use,
      status_txt = "skip_zero_variance"
    ))
  }
  
  df <- data.frame(source = source_vec, mediator = medi_vec, outcome = out_vec)
  
  tryCatch({
    fit_m <- stats::lm(mediator ~ source, data = df, model = FALSE, x = FALSE, y = FALSE)
    fit_y <- stats::lm(outcome ~ source + mediator, data = df, model = FALSE, x = FALSE, y = FALSE)
    fit_total <- stats::lm(outcome ~ source, data = df, model = FALSE, x = FALSE, y = FALSE)
    
    med_obj <- mediation::mediate(
      model.m = fit_m,
      model.y = fit_y,
      treat = "source",
      mediator = "mediator",
      boot = TRUE,
      sims = sims,
      long = FALSE,
      use_speed = FALSE
    )
    
    a_info      <- extract_lm_coef(fit_m, "source")
    b_info      <- extract_lm_coef(fit_y, "mediator")
    cprime_info <- extract_lm_coef(fit_y, "source")
    total_info  <- extract_lm_coef(fit_total, "source")
    
    tibble::tibble(
      chain_id = cid,
      source_block = cand_df$source_block[i],
      source_safe = cand_df$source_safe[i],
      mediator_safe = cand_df$mediator_safe[i],
      trait_safe = cand_df$trait_safe[i],
      n = n_use,
      
      a_est = unname(a_info["est"]),
      a_p   = unname(a_info["p"]),
      b_est = unname(b_info["est"]),
      b_p   = unname(b_info["p"]),
      cprime_est = unname(cprime_info["est"]),
      cprime_p   = unname(cprime_info["p"]),
      total_est  = unname(total_info["est"]),
      total_p    = unname(total_info["p"]),
      
      acme = get_med_item(med_obj, c("d.avg", "d0", "d1")),
      acme_p = get_med_item(med_obj, c("d.avg.p", "d0.p", "d1.p")),
      
      ade = get_med_item(med_obj, c("z.avg", "z0", "z1")),
      ade_p = get_med_item(med_obj, c("z.avg.p", "z0.p", "z1.p")),
      
      prop_med = get_med_item(med_obj, c("n.avg", "n0", "n1")),
      prop_med_p = get_med_item(med_obj, c("n.avg.p", "n0.p", "n1.p")),
      
      med_total = get_med_item(med_obj, c("tau.coef")),
      med_total_p = get_med_item(med_obj, c("tau.p")),
      
      status = "ok"
    )
  }, error = function(e) {
    empty_med_row(
      cid = cid,
      source_block = cand_df$source_block[i],
      source_safe = cand_df$source_safe[i],
      mediator_safe = cand_df$mediator_safe[i],
      trait_safe = cand_df$trait_safe[i],
      n_use = n_use,
      status_txt = paste0("error: ", conditionMessage(e))
    )
  })
}

run_mediation_candidates_resume <- function(aligned_source, aligned_mediator, aligned_trait,
                                            cand_df, source_block = "Col",
                                            sims = 1000,
                                            partial_csv,
                                            partial_rds,
                                            save_every_n = 200,
                                            n_cores = 1,
                                            chunk_size = NULL) {
  if (nrow(cand_df) == 0) return(tibble::tibble())
  
  source_feat_names <- setdiff(colnames(aligned_source), "SampleID")
  mediator_feat_names <- setdiff(colnames(aligned_mediator), "SampleID")
  trait_feat_names <- setdiff(colnames(aligned_trait), "SampleID")
  
  cand_df <- cand_df %>%
    dplyr::mutate(
      chain_id = make_chain_id(.),
      source_idx = match(source_safe, source_feat_names),
      mediator_idx = match(mediator_safe, mediator_feat_names),
      trait_idx = match(trait_safe, trait_feat_names)
    )
  
  bad_idx_n <- sum(!is.finite(cand_df$source_idx) | !is.finite(cand_df$mediator_idx) | !is.finite(cand_df$trait_idx))
  if (bad_idx_n > 0) {
    warning(source_block, " candidate chains failed during index mapping and have been removed automatically: ", bad_idx_n)
    cand_df <- cand_df %>%
      dplyr::filter(is.finite(source_idx), is.finite(mediator_idx), is.finite(trait_idx))
  }
  
  existing_res <- NULL
  if (file.exists(partial_rds)) {
    message(source_block, " existing partial results detected; resuming from checkpoint: ", partial_rds)
    existing_res <- readRDS(partial_rds)
    
    if (!"chain_id" %in% colnames(existing_res)) {
      existing_res <- existing_res %>%
        dplyr::mutate(chain_id = paste(source_block, source_safe, mediator_safe, trait_safe, sep = "||"))
    }
  } else {
    existing_res <- tibble::tibble()
  }
  
  done_ids <- if (nrow(existing_res) > 0) unique(existing_res$chain_id) else character(0)
  
  cand_todo <- cand_df %>%
    dplyr::filter(!chain_id %in% done_ids)
  
  message(source_block, " total candidate chains: ", nrow(cand_df))
  message(source_block, " completed chains: ", length(done_ids))
  message(source_block, " chains to run in this session: ", nrow(cand_todo))
  
  if (nrow(cand_todo) == 0) {
    message(source_block, " no additional run is required; all chains have already been completed.")
    return(existing_res)
  }
  
  source_mat <- as.matrix(aligned_source[, setdiff(colnames(aligned_source), "SampleID"), drop = FALSE])
  mediator_mat <- as.matrix(aligned_mediator[, setdiff(colnames(aligned_mediator), "SampleID"), drop = FALSE])
  trait_mat <- as.matrix(aligned_trait[, setdiff(colnames(aligned_trait), "SampleID"), drop = FALSE])
  
  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size < 1) {
    chunk_size <- max(save_every_n, n_cores * 5)
  }
  
  idx_groups <- split(seq_len(nrow(cand_todo)), ceiling(seq_len(nrow(cand_todo)) / chunk_size))
  merged_res <- existing_res
  
  message(source_block, " parallel cores: ", n_cores, "; chains per chunk: ", chunk_size, "; total chunks: ", length(idx_groups))
  
  cl <- NULL
  if (.Platform$OS.type == "windows" && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(mediation)
        library(tibble)
        library(speedglm)
      })
      NULL
    })
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "source_mat", "mediator_mat", "trait_mat", "cand_todo",
        "extract_lm_coef", "get_med_item", "empty_med_row", "run_one_mediation_chain",
        "sims"
      ),
      envir = environment()
    )
  }
  
  done_now <- 0L
  rows_since_save <- 0L
  for (g in seq_along(idx_groups)) {
    t0 <- Sys.time()
    idx <- idx_groups[[g]]
    
    batch_res <- if (!is.null(cl) && n_cores > 1) {
      parallel::parLapplyLB(
        cl,
        as.list(idx),
        function(i) run_one_mediation_chain(i, cand_todo, source_mat, mediator_mat, trait_mat, sims = sims)
      )
    } else if (.Platform$OS.type != "windows" && n_cores > 1) {
      parallel::mclapply(
        idx,
        function(i) run_one_mediation_chain(i, cand_todo, source_mat, mediator_mat, trait_mat, sims = sims),
        mc.cores = n_cores,
        mc.preschedule = FALSE
      )
    } else {
      lapply(idx, function(i) run_one_mediation_chain(i, cand_todo, source_mat, mediator_mat, trait_mat, sims = sims))
    }
    
    batch_res <- dplyr::bind_rows(batch_res)
    merged_res <- dplyr::bind_rows(merged_res, batch_res) %>%
      dplyr::distinct(chain_id, .keep_all = TRUE)
    
    done_now <- done_now + length(idx)
    rows_since_save <- rows_since_save + nrow(batch_res)
    
    if (rows_since_save >= save_every_n || g == length(idx_groups)) {
      save_partial_results(merged_res, partial_csv, partial_rds)
      rows_since_save <- 0L
    }
    
    dt_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(
      source_block, " completed chunk ", g, " / ", length(idx_groups), "; cumulative completed chains: ",
      length(unique(merged_res$chain_id)), " / ", nrow(cand_df),
      " (newly completed in this session: ", done_now, " / ", nrow(cand_todo), ")",
      "; elapsed time for this chunk: ", round(dt_min, 2), " min"
    )
    gc(verbose = FALSE)
  }
  
  if (file.exists(partial_rds)) {
    final_res <- readRDS(partial_rds)
  } else {
    final_res <- merged_res
  }
  final_res
}

attach_original_names <- function(res_df, cand_df) {
  if (nrow(res_df) == 0) return(res_df)
  
  key_df <- cand_df %>%
    dplyr::mutate(chain_id = make_chain_id(.)) %>%
    dplyr::select(
      chain_id,
      source_block,
      source_safe,
      mediator_safe,
      trait_safe,
      source_feature,
      mediator_feature,
      trait_feature
    ) %>%
    dplyr::distinct()
  
  res_df %>%
    dplyr::left_join(key_df, by = c("chain_id", "source_block", "source_safe", "mediator_safe", "trait_safe")) %>%
    dplyr::relocate(source_block, source_feature, mediator_feature, trait_feature, .before = 1)
}

## 5.11 Harmonize and complete key columns
standardize_result_columns <- function(df) {
  need_cols <- c(
    "chain_id",
    "source_block", "source_feature", "mediator_feature", "trait_feature",
    "source_safe", "mediator_safe", "trait_safe",
    "n",
    "a_est", "a_p",
    "b_est", "b_p",
    "cprime_est", "cprime_p",
    "total_est", "total_p",
    "acme", "acme_p",
    "ade", "ade_p",
    "prop_med", "prop_med_p",
    "med_total", "med_total_p",
    "status",
    "final_sig"
  )
  
  for (cc in need_cols) {
    if (!cc %in% colnames(df)) {
      if (cc %in% c("chain_id", "source_block", "source_feature", "mediator_feature", "trait_feature",
                    "source_safe", "mediator_safe", "trait_safe", "status", "final_sig")) {
        df[[cc]] <- NA_character_
      } else {
        df[[cc]] <- NA_real_
      }
    }
  }
  
  df <- df[, need_cols, drop = FALSE]
  df
}

## 5.12 Recalculate final_sig
add_final_sig <- function(df) {
  if (nrow(df) == 0) {
    df$final_sig <- character(0)
    return(df)
  }
  
  df %>%
    dplyr::mutate(
      final_sig = ifelse(
        status == "ok" &
          is.finite(a_p) & a_p < 0.05 &
          is.finite(b_p) & b_p < 0.05 &
          is.finite(acme_p) & acme_p < 0.05 &
          is.finite(prop_med) & prop_med > 0 & prop_med < 1 &
          is.finite(prop_med_p) & prop_med_p < 0.05,
        "yes", "no"
      )
    )
}

############################
## 6. Read and preprocess data
############################
col_obj <- prep_microbe_block_pairwise_style(
  filepath = paths$col,
  sheet = sheet_cfg$col_sheet,
  block_name = "Col",
  prevalence_cut = prevalence_cut,
  min_nonzero_n = min_nonzero_n
)

liver_obj <- prep_liver_expression(
  filepath = paths$liver,
  sheet = sheet_cfg$liver_sheet
)

blood_df <- prep_blood_traits(
  filepath = paths$blood,
  sheet = sheet_cfg$blood_sheet
)

tail_df <- prep_tail_traits(
  filepath = paths$tail,
  sheet = sheet_cfg$tail_sheet
)

trait_obj <- prep_all_traits(
  blood_df = blood_df,
  tail_df  = tail_df
)

############################
## 7. Align samples
############################
col_aln <- align_three_blocks(
  source_df   = col_obj$data,
  mediator_df = liver_obj$data,
  trait_df    = trait_obj$data
)

############################
## 8. Build all candidate chains for Col
############################
col_cand <- build_all_candidates(
  source_df    = col_aln$source,
  mediator_df  = col_aln$mediator,
  trait_df     = col_aln$trait,
  source_map   = col_obj$map,
  mediator_map = liver_obj$map,
  trait_map    = trait_obj$map,
  source_block = "Col"
)

utils::write.csv(col_cand, file.path(out_dir, "col_all_candidates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(col_cand, file.path(out_dir, "col_all_candidates.rds"))

############################
## 9. Formal mediation analysis for Col with checkpoint-based resumption
############################
col_res <- run_mediation_candidates_resume(
  aligned_source   = col_aln$source,
  aligned_mediator = liver_obj$data %>%
    dplyr::filter(SampleID %in% col_aln$ids) %>%
    dplyr::arrange(match(SampleID, col_aln$ids)),
  aligned_trait    = trait_obj$data %>%
    dplyr::filter(SampleID %in% col_aln$ids) %>%
    dplyr::arrange(match(SampleID, col_aln$ids)),
  cand_df = col_cand,
  source_block = "Col",
  sims = boot_sims,
  partial_csv = checkpoint_files$col_partial_csv,
  partial_rds = checkpoint_files$col_partial_rds,
  save_every_n = save_every_n,
  n_cores = n_cores,
  chunk_size = chunk_size
)

utils::write.csv(col_res, checkpoint_files$col_final_csv, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(col_res, checkpoint_files$col_final_rds)

############################
## 10. Restore original feature names
############################
col_res2 <- attach_original_names(col_res, col_cand)
col_res2 <- add_final_sig(col_res2)
col_res2 <- standardize_result_columns(col_res2)

############################
## 11. Read the previous combined results (Rum + Ile)
############################
if (!file.exists(paths$old_combined)) {
  stop("The previous combined file was not found: ", paths$old_combined)
}

old_combined <- utils::read.csv(paths$old_combined, check.names = FALSE, stringsAsFactors = FALSE)
old_combined <- standardize_result_columns(old_combined)
old_combined <- add_final_sig(old_combined)

############################
## 12. Merge previous results and Col results
############################
combined_res <- dplyr::bind_rows(old_combined, col_res2)

combined_res <- combined_res %>%
  dplyr::distinct(chain_id, .keep_all = TRUE) %>%
  add_final_sig() %>%
  dplyr::arrange(dplyr::desc(final_sig == "yes"), acme_p, source_block)

utils::write.csv(combined_res, checkpoint_files$combined_csv, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(combined_res, checkpoint_files$combined_rds)

############################
## 13. Summary information table
############################
summary_info <- tibble::tibble(
  item = c(
    "prevalence_cut",
    "min_nonzero_n",
    "candidate_method",
    "boot_sims",
    "save_every_n",
    "n_cores",
    "chunk_size",
    "Col_features_after",
    "Liver_genes_after",
    "Trait_count_after",
    "Col_common_samples",
    "Col_all_candidates",
    "Col_mediation_rows",
    "Old_combined_rows",
    "New_combined_rows",
    "New_combined_final_sig_rows"
  ),
  value = c(
    prevalence_cut,
    ifelse(is.null(min_nonzero_n), NA, min_nonzero_n),
    "No prescreen; all source-mediator-trait combinations entered mediation directly",
    boot_sims,
    save_every_n,
    n_cores,
    chunk_size,
    col_obj$n_feature_after,
    liver_obj$n_gene_after,
    ncol(trait_obj$data) - 1,
    col_aln$n_sample,
    nrow(col_cand),
    nrow(col_res2),
    nrow(old_combined),
    nrow(combined_res),
    ifelse(nrow(combined_res) == 0, 0, sum(combined_res$final_sig == "yes", na.rm = TRUE))
  )
)

saveRDS(summary_info, checkpoint_files$summary_rds)

############################
## 14. Export Excel file
############################
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "README")
readme_text <- data.frame(
  Section = c(
    "Analysis framework",
    "Microbial preprocessing",
    "Liver gene processing",
    "Phenotype processing",
    "Candidate-chain inclusion rule",
    "Formal mediation analysis",
    "Parallel acceleration",
    "Final significance criteria",
    "Checkpoint-resumption note",
    "Result integration note"
  ),
  Detail = c(
    "This script adds Col -> liver gene -> 8 traits and merges the completed results with the existing Rum + Ile combined results",
    paste0("Consistent with the original script: prevalence filtering (proportion threshold ", prevalence_cut, ") -> CZM zero replacement -> CLR -> z-score"),
    "formal_log2TPM -> z-score",
    "Five blood traits and three tail-fat traits are merged and then z-score standardized",
    "No pre-screening is applied; all source x mediator x trait combinations enter the formal mediation analysis",
    paste0("mediate(boot=TRUE, sims=", boot_sims, ")"),
    paste0("Chunk-wise parallelization; n_cores=", n_cores, "; chunk_size=", chunk_size),
    "a_p<0.05, b_p<0.05, acme_p<0.05, 0<prop_med<1, prop_med_p<0.05",
    paste0("Partial results are saved automatically every ", save_every_n, " chains; rerunning the script will resume from the checkpoint"),
    "Final combined table = previous Rum + Ile combined_mediation_results.csv + new Col results"
  ),
  stringsAsFactors = FALSE
)
openxlsx::writeData(wb, "README", readme_text)

openxlsx::addWorksheet(wb, "summary")
openxlsx::writeData(wb, "summary", summary_info)

openxlsx::addWorksheet(wb, "col_feature_map")
openxlsx::writeData(wb, "col_feature_map", col_obj$map)

openxlsx::addWorksheet(wb, "liver_gene_map")
openxlsx::writeData(wb, "liver_gene_map", liver_obj$map)

openxlsx::addWorksheet(wb, "trait_map")
openxlsx::writeData(wb, "trait_map", trait_obj$map)

openxlsx::addWorksheet(wb, "col_prevalence_info")
openxlsx::writeData(wb, "col_prevalence_info", col_obj$prevalence_info)

openxlsx::addWorksheet(wb, "col_all_candidates")
openxlsx::writeData(wb, "col_all_candidates", col_cand)

openxlsx::addWorksheet(wb, "col_results")
openxlsx::writeData(wb, "col_results", col_res2)

openxlsx::addWorksheet(wb, "old_combined_results")
openxlsx::writeData(wb, "old_combined_results", old_combined)

openxlsx::addWorksheet(wb, "combined_results")
openxlsx::writeData(wb, "combined_results", combined_res)

sig_only <- combined_res %>% dplyr::filter(final_sig == "yes")
openxlsx::addWorksheet(wb, "combined_sig_only")
openxlsx::writeData(wb, "combined_sig_only", sig_only)

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

openxlsx::saveWorkbook(wb, checkpoint_files$excel_out, overwrite = TRUE)

############################
## 15. Console messages
############################
cat("\n=============================\n")
cat("Colon mediation analysis completed using the original logic with the added Col block.\n")
cat("=============================\n")
cat("Result directory: ", out_dir, "\n")
cat("Excel summary file: ", checkpoint_files$excel_out, "\n\n")

cat("Key statistics:\n")
cat("Parallel cores: ", n_cores, "\n")
cat("Chains per chunk: ", chunk_size, "\n")
cat("Col features after prevalence filtering: ", col_obj$n_feature_after, "\n")
cat("Number of liver genes: ", liver_obj$n_gene_after, "\n")
cat("Total Col candidate chains: ", nrow(col_cand), "\n")
cat("Col formal result rows: ", nrow(col_res2), "\n")
cat("Previous combined result rows: ", nrow(old_combined), "\n")
cat("New combined result rows: ", nrow(combined_res), "\n")
cat("Final significant chains in the new combined results: ", sum(combined_res$final_sig == "yes", na.rm = TRUE), "\n")
cat("\n")