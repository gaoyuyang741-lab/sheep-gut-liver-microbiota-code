############################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################
## 0. Install / load packages
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
## 1. Path settings
############################
# Input and output directories
# Please place the required cross-compartment coordinated genera input files in "data/coordinated_genera/input".
# Output files will be saved in "results/coordinated_genera".

base_dir <- file.path("data", "coordinated_genera", "input")
out_dir  <- file.path("results", "coordinated_genera")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rum_fp <- file.path(base_dir, "Rum_genus_abundance.xlsx")
ile_fp <- file.path(base_dir, "Ile_genus_abundance.xlsx")
col_fp <- file.path(base_dir, "Colon_genus_abundance.xlsx")

process_xlsx <- file.path(out_dir, "01_continuity_structure_analysis_process_file.xlsx")
result_xlsx  <- file.path(out_dir, "02_continuity_structure_analysis_result_file.xlsx")

############################
## 2. Parameter settings
############################
# Baseline threshold: prevalence >= 0.5 in at least one compartment
prevalence_core_cut <- 0.50

# Traceability threshold for the other end: prevalence >= 0.2 or nonzero sample count >= 6
prevalence_trace_cut <- 0.20
nonzero_trace_cut    <- 6

# Significance threshold
q_cut_main <- 0.05

# Strong-correlation threshold
abs_r_cut_strict <- 0.30

# Supporting-evidence threshold for the third pair
abs_r_cut_support <- 0.20

############################
## 3. Basic functions
############################

clean_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x
}

read_genus_table <- function(fp, sheet = 1) {
  df <- readxl::read_excel(fp, sheet = sheet)
  df <- as.data.frame(df)
  
  if (!"Genus" %in% colnames(df)) {
    stop("No 'Genus' column was found in file: ", fp)
  }
  
  df <- df[!is.na(df$Genus) & df$Genus != "", , drop = FALSE]
  
  sample_cols <- setdiff(colnames(df), "Genus")
  sample_cols_clean <- clean_sample_id(sample_cols)
  colnames(df) <- c("Genus", sample_cols_clean)
  
  for (cc in sample_cols_clean) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }
  df[is.na(df)] <- 0
  
  # Merge duplicated genus names
  df <- df |>
    dplyr::group_by(Genus) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.matrix(df[, sample_cols_clean, drop = FALSE])
  rownames(mat) <- df$Genus
  mode(mat) <- "numeric"
  
  return(mat)  # genus x sample
}

calc_genus_stats <- function(mat) {
  data.frame(
    Genus = rownames(mat),
    prevalence = rowMeans(mat > 0, na.rm = TRUE),
    nonzero_n = rowSums(mat > 0, na.rm = TRUE),
    mean_abundance = rowMeans(mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

filter_single_site_core <- function(mat, prevalence_core_cut = 0.5) {
  st <- calc_genus_stats(mat)
  st$kept_single_site_core <- st$prevalence >= prevalence_core_cut
  list(
    mat = mat[st$kept_single_site_core, , drop = FALSE],
    info = st
  )
}

clr_transform_czm <- function(mat_genus_by_sample) {
  if (nrow(mat_genus_by_sample) < 2) {
    stop("Fewer than two genera remained after filtering; CLR cannot be performed.")
  }
  
  x <- t(mat_genus_by_sample)  # sample x genus
  
  x_nozero <- zCompositions::cmultRepl(
    X = x,
    method = "CZM",
    output = "p-counts",
    label = 0
  )
  
  x_clr <- compositions::clr(x_nozero, base = exp(1))
  x_clr <- as.matrix(x_clr)
  return(x_clr)  # sample x genus
}

clr_mat_to_df <- function(clr_mat) {
  out <- as.data.frame(clr_mat, check.names = FALSE)
  out <- tibble::rownames_to_column(out, "SampleID")
  out$SampleID <- clean_sample_id(out$SampleID)
  out
}

align_two_blocks <- function(df_x, df_y) {
  df_x$SampleID <- clean_sample_id(df_x$SampleID)
  df_y$SampleID <- clean_sample_id(df_y$SampleID)
  
  common_ids <- intersect(df_x$SampleID, df_y$SampleID)
  if (length(common_ids) < 3) {
    stop("Too few common samples for analysis.")
  }
  
  df_x2 <- df_x |>
    dplyr::filter(SampleID %in% common_ids) |>
    dplyr::arrange(SampleID)
  
  df_y2 <- df_y |>
    dplyr::filter(SampleID %in% common_ids) |>
    dplyr::arrange(SampleID)
  
  if (!identical(df_x2$SampleID, df_y2$SampleID)) {
    stop("Sample alignment failed.")
  }
  
  list(x = df_x2, y = df_y2, common_ids = common_ids)
}

get_pair_entry_table <- function(stats1, stats2,
                                 site1, site2,
                                 prevalence_core_cut = 0.5,
                                 prevalence_trace_cut = 0.2,
                                 nonzero_trace_cut = 6) {
  s1 <- stats1 |>
    dplyr::rename(
      prevalence_1 = prevalence,
      nonzero_n_1 = nonzero_n,
      mean_abundance_1 = mean_abundance
    )
  
  s2 <- stats2 |>
    dplyr::rename(
      prevalence_2 = prevalence,
      nonzero_n_2 = nonzero_n,
      mean_abundance_2 = mean_abundance
    )
  
  out <- dplyr::full_join(s1, s2, by = "Genus")
  
  out$prevalence_1[is.na(out$prevalence_1)] <- 0
  out$prevalence_2[is.na(out$prevalence_2)] <- 0
  out$nonzero_n_1[is.na(out$nonzero_n_1)] <- 0
  out$nonzero_n_2[is.na(out$nonzero_n_2)] <- 0
  out$mean_abundance_1[is.na(out$mean_abundance_1)] <- 0
  out$mean_abundance_2[is.na(out$mean_abundance_2)] <- 0
  
  out$core_1  <- out$prevalence_1 >= prevalence_core_cut
  out$core_2  <- out$prevalence_2 >= prevalence_core_cut
  out$trace_1 <- (out$prevalence_1 >= prevalence_trace_cut) | (out$nonzero_n_1 >= nonzero_trace_cut)
  out$trace_2 <- (out$prevalence_2 >= prevalence_trace_cut) | (out$nonzero_n_2 >= nonzero_trace_cut)
  
  # At least one end is core and the other end is traceable
  out$pass_pair_entry <- (out$core_1 & out$trace_2) | (out$core_2 & out$trace_1)
  
  out$entry_type <- dplyr::case_when(
    out$core_1 & out$core_2 ~ paste0(site1, "_core + ", site2, "_core"),
    out$core_1 & !out$core_2 & out$trace_2 ~ paste0(site1, "_core + ", site2, "_traceable"),
    !out$core_1 & out$core_2 & out$trace_1 ~ paste0(site2, "_core + ", site1, "_traceable"),
    TRUE ~ "not included in analysis"
  )
  
  out
}

calc_pair_correlations <- function(df_x, df_y, entry_df, pair_name,
                                   q_cut_main = 0.05,
                                   abs_r_cut_strict = 0.30,
                                   method_main = "pearson",
                                   method_sens = "spearman") {
  aligned <- align_two_blocks(df_x, df_y)
  x <- aligned$x
  y <- aligned$y
  
  genus_x <- setdiff(colnames(x), "SampleID")
  genus_y <- setdiff(colnames(y), "SampleID")
  common_genus <- intersect(genus_x, genus_y)
  
  entry_df2 <- entry_df |>
    dplyr::filter(pass_pair_entry) |>
    dplyr::filter(Genus %in% common_genus)
  
  if (nrow(entry_df2) == 0) {
    stop(pair_name, ": no analyzable genus remained after applying the new inclusion rules.")
  }
  
  res_list <- vector("list", nrow(entry_df2))
  
  for (i in seq_len(nrow(entry_df2))) {
    g <- entry_df2$Genus[i]
    
    vx <- suppressWarnings(as.numeric(x[[g]]))
    vy <- suppressWarnings(as.numeric(y[[g]]))
    
    ok <- is.finite(vx) & is.finite(vy)
    n_ok <- sum(ok)
    
    r_main <- NA_real_
    p_main <- NA_real_
    r_sens <- NA_real_
    p_sens <- NA_real_
    
    if (n_ok >= 3) {
      ct_main <- suppressWarnings(
        cor.test(vx[ok], vy[ok], method = method_main)
      )
      r_main <- unname(ct_main$estimate)
      p_main <- ct_main$p.value
      
      ct_sens <- suppressWarnings(
        cor.test(vx[ok], vy[ok], method = method_sens, exact = FALSE)
      )
      r_sens <- unname(ct_sens$estimate)
      p_sens <- ct_sens$p.value
    }
    
    res_list[[i]] <- data.frame(
      Genus = g,
      n = n_ok,
      r_main = r_main,
      p_main = p_main,
      r_sens = r_sens,
      p_sens = p_sens,
      pair = pair_name,
      stringsAsFactors = FALSE
    )
  }
  
  res <- dplyr::bind_rows(res_list)
  res$q_main <- p.adjust(res$p_main, method = "BH")
  res$q_sens <- p.adjust(res$p_sens, method = "BH")
  res$abs_r_main <- abs(res$r_main)
  res$abs_r_sens <- abs(res$r_sens)
  
  res <- res |>
    dplyr::left_join(entry_df2, by = "Genus")
  
  # Main-analysis direction
  res$main_direction <- ifelse(
    is.na(res$r_main), NA,
    ifelse(res$r_main > 0, "positive",
           ifelse(res$r_main < 0, "substitution", "zero correlation"))
  )
  
  # Sensitivity-analysis direction
  res$sens_direction <- ifelse(
    is.na(res$r_sens), NA,
    ifelse(res$r_sens > 0, "positive",
           ifelse(res$r_sens < 0, "substitution", "zero correlation"))
  )
  
  # Strict significance: Pearson q < 0.05
  res$main_sig <- !is.na(res$q_main) & (res$q_main < q_cut_main)
  
  # Strict strong evidence: Pearson q < 0.05 and |r| >= 0.3
  res$main_sig_strict <- !is.na(res$q_main) &
    (res$q_main < q_cut_main) &
    (abs(res$r_main) >= abs_r_cut_strict)
  
  # Sensitivity support: Spearman direction is consistent and |r| >= 0.2
  res$sens_support <- !is.na(res$r_sens) & !is.na(res$r_main) &
    (sign(res$r_sens) == sign(res$r_main)) &
    (abs(res$r_sens) >= abs_r_cut_support)
  
  # Evidence grading for each pair
  res$pair_evidence_level <- dplyr::case_when(
    res$main_sig_strict & res$sens_support ~ "strong evidence",
    res$main_sig & res$sens_support ~ "relatively strong evidence",
    res$main_sig & !res$sens_support ~ "main-analysis evidence",
    !res$main_sig & res$sens_support & !is.na(res$r_main) & (abs(res$r_main) >= abs_r_cut_support) ~ "trend support",
    TRUE ~ "no evidence"
  )
  
  # Continuity type for each pair
  res$continuity_type <- dplyr::case_when(
    res$pair_evidence_level != "no evidence" & res$r_main > 0 ~ "positive-continuity type",
    res$pair_evidence_level != "no evidence" & res$r_main < 0 ~ "substitution-continuity type",
    TRUE ~ "undetermined"
  )
  
  res <- res |>
    dplyr::arrange(q_main, dplyr::desc(abs_r_main))
  
  res
}

make_pair_summary <- function(res_df, pair_name) {
  sig_all <- res_df |>
    dplyr::filter(main_sig)
  
  sig_strict <- res_df |>
    dplyr::filter(main_sig_strict)
  
  data.frame(
    pair = pair_name,
    total_tested_genus = nrow(res_df),
    strong_evidence = sum(res_df$pair_evidence_level == "strong evidence", na.rm = TRUE),
    relatively_strong_evidence = sum(res_df$pair_evidence_level == "relatively strong evidence", na.rm = TRUE),
    main_only_evidence = sum(res_df$pair_evidence_level == "main-analysis evidence", na.rm = TRUE),
    trend_support = sum(res_df$pair_evidence_level == "trend support", na.rm = TRUE),
    main_sig_total = nrow(sig_all),
    main_sig_strict = nrow(sig_strict),
    same_direction_main = sum(sig_all$r_main > 0, na.rm = TRUE),
    replacement_main = sum(sig_all$r_main < 0, na.rm = TRUE),
    median_abs_r_main_sig = ifelse(nrow(sig_all) > 0, median(abs(sig_all$r_main), na.rm = TRUE), NA_real_),
    stringsAsFactors = FALSE
  )
}

classify_structure_across_three <- function(pattern_df) {
  pattern_df$triple_strict_support <- with(
    pattern_df,
    sig_RI_strict & sig_RC_strict & sig_IC_strict
  )
  
  pattern_df$triple_main_support <- with(
    pattern_df,
    sig_RI_main & sig_RC_main & sig_IC_main
  )
  
  pattern_df$two_main_one_trend <- with(
    pattern_df,
    (
      sig_RI_main & sig_RC_main & (!sig_IC_main) & trend_IC
    ) |
      (
        sig_RI_main & sig_IC_main & (!sig_RC_main) & trend_RC
      ) |
      (
        sig_RC_main & sig_IC_main & (!sig_RI_main) & trend_RI
      )
  )
  
  pattern_df$any_two_main <- with(
    pattern_df,
    (sig_RI_main + sig_RC_main + sig_IC_main) >= 2
  )
  
  pattern_df$any_one_strict <- with(
    pattern_df,
    sig_RI_strict | sig_RC_strict | sig_IC_strict
  )
  
  pattern_df$all_same_sign_main <- with(
    pattern_df,
    !is.na(sign_RI) & !is.na(sign_RC) & !is.na(sign_IC) &
      (sign_RI == sign_RC) & (sign_RC == sign_IC) &
      sign_RI != 0
  )
  
  pattern_df$all_same_sign_available <- with(
    pattern_df,
    {
      s1 <- sign_RI
      s2 <- sign_RC
      s3 <- sign_IC
      all_non_na <- !is.na(s1) & !is.na(s2) & !is.na(s3)
      all_non_na & (s1 == s2) & (s2 == s3) & s1 != 0
    }
  )
  
  pattern_df$structure_grade <- dplyr::case_when(
    pattern_df$triple_strict_support & pattern_df$all_same_sign_available ~ "three-compartment continuity: strong evidence",
    pattern_df$triple_main_support & pattern_df$all_same_sign_available ~ "three-compartment continuity: main evidence",
    pattern_df$two_main_one_trend ~ "three-compartment continuity: supporting evidence",
    pattern_df$sig_RI_main & pattern_df$sig_RC_main & !pattern_df$sig_IC_main ~ "rumen dual-link evidence",
    pattern_df$sig_RI_main & pattern_df$sig_IC_main & !pattern_df$sig_RC_main ~ "ileal bridge evidence",
    pattern_df$sig_RC_main & pattern_df$sig_IC_main & !pattern_df$sig_RI_main ~ "hindgut dual-link evidence",
    pattern_df$sig_RI_strict & !pattern_df$sig_RC_main & !pattern_df$sig_IC_main ~ "anterior continuity: strong evidence",
    pattern_df$sig_RC_strict & !pattern_df$sig_RI_main & !pattern_df$sig_IC_main ~ "foregut-hindgut coupling: strong evidence",
    pattern_df$sig_IC_strict & !pattern_df$sig_RI_main & !pattern_df$sig_RC_main ~ "posterior continuity: strong evidence",
    pattern_df$sig_RI_main & !pattern_df$sig_RC_main & !pattern_df$sig_IC_main ~ "anterior continuity: main evidence",
    pattern_df$sig_RC_main & !pattern_df$sig_RI_main & !pattern_df$sig_IC_main ~ "foregut-hindgut coupling: main evidence",
    pattern_df$sig_IC_main & !pattern_df$sig_RI_main & !pattern_df$sig_RC_main ~ "posterior continuity: main evidence",
    TRUE ~ "no clear structural evidence"
  )
  
  pattern_df
}

make_pattern_table <- function(res_RI, res_RC, res_IC,
                               rum_stats, ile_stats, col_stats) {
  all_genus <- sort(unique(c(res_RI$Genus, res_RC$Genus, res_IC$Genus)))
  out <- data.frame(Genus = all_genus, stringsAsFactors = FALSE)
  
  tmp_RI <- res_RI |>
    dplyr::select(
      Genus,
      r_RI = r_main, q_RI = q_main, abs_r_RI = abs_r_main,
      r_RI_sens = r_sens, q_RI_sens = q_sens,
      type_RI = continuity_type,
      evidence_RI = pair_evidence_level,
      sig_RI_main = main_sig,
      sig_RI_strict = main_sig_strict
    )
  
  tmp_RC <- res_RC |>
    dplyr::select(
      Genus,
      r_RC = r_main, q_RC = q_main, abs_r_RC = abs_r_main,
      r_RC_sens = r_sens, q_RC_sens = q_sens,
      type_RC = continuity_type,
      evidence_RC = pair_evidence_level,
      sig_RC_main = main_sig,
      sig_RC_strict = main_sig_strict
    )
  
  tmp_IC <- res_IC |>
    dplyr::select(
      Genus,
      r_IC = r_main, q_IC = q_main, abs_r_IC = abs_r_main,
      r_IC_sens = r_sens, q_IC_sens = q_sens,
      type_IC = continuity_type,
      evidence_IC = pair_evidence_level,
      sig_IC_main = main_sig,
      sig_IC_strict = main_sig_strict
    )
  
  out <- out |>
    dplyr::left_join(tmp_RI, by = "Genus") |>
    dplyr::left_join(tmp_RC, by = "Genus") |>
    dplyr::left_join(tmp_IC, by = "Genus")
  
  for (cc in c("sig_RI_main", "sig_RI_strict", "sig_RC_main", "sig_RC_strict", "sig_IC_main", "sig_IC_strict")) {
    out[[cc]][is.na(out[[cc]])] <- FALSE
  }
  
  out$trend_RI <- !is.na(out$r_RI) & !out$sig_RI_main & (abs(out$r_RI) >= abs_r_cut_support)
  out$trend_RC <- !is.na(out$r_RC) & !out$sig_RC_main & (abs(out$r_RC) >= abs_r_cut_support)
  out$trend_IC <- !is.na(out$r_IC) & !out$sig_IC_main & (abs(out$r_IC) >= abs_r_cut_support)
  
  out$sign_RI <- ifelse(is.na(out$r_RI), NA, sign(out$r_RI))
  out$sign_RC <- ifelse(is.na(out$r_RC), NA, sign(out$r_RC))
  out$sign_IC <- ifelse(is.na(out$r_IC), NA, sign(out$r_IC))
  
  out$structure_code_main <- paste0(
    as.integer(out$sig_RI_main),
    as.integer(out$sig_RC_main),
    as.integer(out$sig_IC_main)
  )
  
  out$structure_code_strict <- paste0(
    as.integer(out$sig_RI_strict),
    as.integer(out$sig_RC_strict),
    as.integer(out$sig_IC_strict)
  )
  
  out$pair_signature <- paste(
    ifelse(is.na(out$evidence_RI), "RI:NA", paste0("RI:", out$evidence_RI)),
    ifelse(is.na(out$evidence_RC), "RC:NA", paste0("RC:", out$evidence_RC)),
    ifelse(is.na(out$evidence_IC), "IC:NA", paste0("IC:", out$evidence_IC)),
    sep = " | "
  )
  
  sR <- rum_stats |>
    dplyr::rename(
      R_prevalence = prevalence,
      R_nonzero_n = nonzero_n,
      R_mean_abundance = mean_abundance
    )
  
  sI <- ile_stats |>
    dplyr::rename(
      I_prevalence = prevalence,
      I_nonzero_n = nonzero_n,
      I_mean_abundance = mean_abundance
    )
  
  sC <- col_stats |>
    dplyr::rename(
      C_prevalence = prevalence,
      C_nonzero_n = nonzero_n,
      C_mean_abundance = mean_abundance
    )
  
  out <- out |>
    dplyr::left_join(sR, by = "Genus") |>
    dplyr::left_join(sI, by = "Genus") |>
    dplyr::left_join(sC, by = "Genus")
  
  out$abundance_pattern <- apply(out, 1, function(z) {
    mean_R <- suppressWarnings(as.numeric(z["R_mean_abundance"]))
    mean_I <- suppressWarnings(as.numeric(z["I_mean_abundance"]))
    mean_C <- suppressWarnings(as.numeric(z["C_mean_abundance"]))
    
    vals <- c(R = mean_R, I = mean_I, C = mean_C)
    if (any(is.na(vals))) return(NA_character_)
    
    ord <- names(sort(vals, decreasing = TRUE))
    tol_base <- max(vals, na.rm = TRUE)
    tol <- ifelse(is.finite(tol_base) && tol_base > 0, 0.2 * tol_base, 0)
    
    if (abs(mean_R - mean_C) <= tol && mean_I < mean_R && mean_I < mean_C) {
      return("Rum ≈ Col higher than Ile")
    }
    if (abs(mean_R - mean_C) <= tol && mean_I > mean_R && mean_I > mean_C) {
      return("Ile higher than Rum ≈ Col")
    }
    
    paste(ord, collapse = " > ")
  })
  
  out <- classify_structure_across_three(out)
  out
}

extract_clr_df_by_genus <- function(clr_df, genus_vec) {
  keep <- intersect(genus_vec, setdiff(colnames(clr_df), "SampleID"))
  clr_df[, c("SampleID", keep), drop = FALSE]
}

############################
## 4. Read data
############################
cat("======================================\n")
cat("Reading data...\n")
cat("======================================\n")

rum_raw_mat <- read_genus_table(rum_fp, sheet = 1)
ile_raw_mat <- read_genus_table(ile_fp, sheet = 1)
col_raw_mat <- read_genus_table(col_fp, sheet = 1)

############################
## 5. Sample check
############################
cat("======================================\n")
cat("Checking whether sample columns are consistent...\n")
cat("======================================\n")

rum_samples <- clean_sample_id(colnames(rum_raw_mat))
ile_samples <- clean_sample_id(colnames(ile_raw_mat))
col_samples <- clean_sample_id(colnames(col_raw_mat))

if (!identical(sort(rum_samples), sort(ile_samples)) ||
    !identical(sort(rum_samples), sort(col_samples))) {
  stop("Sample columns are inconsistent across the three compartments; please check first.")
}

sample_order <- sort(rum_samples)
rum_raw_mat <- rum_raw_mat[, sample_order, drop = FALSE]
ile_raw_mat <- ile_raw_mat[, sample_order, drop = FALSE]
col_raw_mat <- col_raw_mat[, sample_order, drop = FALSE]

############################
## 6. Compartment-specific core-genus filtering (still based on prevalence 0.5)
############################
cat("======================================\n")
cat("Filtering compartment-specific core genera (prevalence >= 0.5)...\n")
cat("======================================\n")

rum_core <- filter_single_site_core(rum_raw_mat, prevalence_core_cut = prevalence_core_cut)
ile_core <- filter_single_site_core(ile_raw_mat, prevalence_core_cut = prevalence_core_cut)
col_core <- filter_single_site_core(col_raw_mat, prevalence_core_cut = prevalence_core_cut)

rum_raw_core <- rum_core$mat
ile_raw_core <- ile_core$mat
col_raw_core <- col_core$mat

rum_stats_all <- calc_genus_stats(rum_raw_mat)
ile_stats_all <- calc_genus_stats(ile_raw_mat)
col_stats_all <- calc_genus_stats(col_raw_mat)

############################
## 7. CLR transformation
## Note: CLR is still performed using compartment-specific core-genus matrices
############################
cat("======================================\n")
cat("CZM zero replacement + CLR...\n")
cat("======================================\n")

rum_clr <- clr_transform_czm(rum_raw_core)
ile_clr <- clr_transform_czm(ile_raw_core)
col_clr <- clr_transform_czm(col_raw_core)

rum_clr_df <- clr_mat_to_df(rum_clr)
ile_clr_df <- clr_mat_to_df(ile_clr)
col_clr_df <- clr_mat_to_df(col_clr)

############################
## 8. Build three inclusion tables using the "one-end core + other-end traceable" rule
############################
cat("======================================\n")
cat("Building paired inclusion-rule tables...\n")
cat("======================================\n")

entry_RI <- get_pair_entry_table(
  stats1 = rum_stats_all,
  stats2 = ile_stats_all,
  site1 = "Rumen",
  site2 = "Ileum",
  prevalence_core_cut = prevalence_core_cut,
  prevalence_trace_cut = prevalence_trace_cut,
  nonzero_trace_cut = nonzero_trace_cut
)

entry_RC <- get_pair_entry_table(
  stats1 = rum_stats_all,
  stats2 = col_stats_all,
  site1 = "Rumen",
  site2 = "Colon",
  prevalence_core_cut = prevalence_core_cut,
  prevalence_trace_cut = prevalence_trace_cut,
  nonzero_trace_cut = nonzero_trace_cut
)

entry_IC <- get_pair_entry_table(
  stats1 = ile_stats_all,
  stats2 = col_stats_all,
  site1 = "Ileum",
  site2 = "Colon",
  prevalence_core_cut = prevalence_core_cut,
  prevalence_trace_cut = prevalence_trace_cut,
  nonzero_trace_cut = nonzero_trace_cut
)

############################
## 9. Paired correlation analysis
## Main analysis: Pearson
## Sensitivity analysis: Spearman
############################
cat("======================================\n")
cat("Calculating paired continuity correlations...\n")
cat("======================================\n")

res_RI <- calc_pair_correlations(
  df_x = rum_clr_df,
  df_y = ile_clr_df,
  entry_df = entry_RI,
  pair_name = "Rumen_vs_Ileum",
  q_cut_main = q_cut_main,
  abs_r_cut_strict = abs_r_cut_strict,
  method_main = "pearson",
  method_sens = "spearman"
)

res_RC <- calc_pair_correlations(
  df_x = rum_clr_df,
  df_y = col_clr_df,
  entry_df = entry_RC,
  pair_name = "Rumen_vs_Colon",
  q_cut_main = q_cut_main,
  abs_r_cut_strict = abs_r_cut_strict,
  method_main = "pearson",
  method_sens = "spearman"
)

res_IC <- calc_pair_correlations(
  df_x = ile_clr_df,
  df_y = col_clr_df,
  entry_df = entry_IC,
  pair_name = "Ileum_vs_Colon",
  q_cut_main = q_cut_main,
  abs_r_cut_strict = abs_r_cut_strict,
  method_main = "pearson",
  method_sens = "spearman"
)

############################
## 10. Extract key results
############################
sig_RI_main <- res_RI |>
  dplyr::filter(main_sig)

sig_RC_main <- res_RC |>
  dplyr::filter(main_sig)

sig_IC_main <- res_IC |>
  dplyr::filter(main_sig)

sig_RI_strict <- res_RI |>
  dplyr::filter(main_sig_strict)

sig_RC_strict <- res_RC |>
  dplyr::filter(main_sig_strict)

sig_IC_strict <- res_IC |>
  dplyr::filter(main_sig_strict)

RI_same <- sig_RI_main |>
  dplyr::filter(r_main > 0)

RI_replace <- sig_RI_main |>
  dplyr::filter(r_main < 0)

RC_same <- sig_RC_main |>
  dplyr::filter(r_main > 0)

RC_replace <- sig_RC_main |>
  dplyr::filter(r_main < 0)

IC_same <- sig_IC_main |>
  dplyr::filter(r_main > 0)

IC_replace <- sig_IC_main |>
  dplyr::filter(r_main < 0)

############################
## 11. Three-compartment structural pattern table
############################
cat("======================================\n")
cat("Building the three-compartment structural pattern table...\n")
cat("======================================\n")

pattern_df <- make_pattern_table(
  res_RI = res_RI,
  res_RC = res_RC,
  res_IC = res_IC,
  rum_stats = rum_stats_all,
  ile_stats = ile_stats_all,
  col_stats = col_stats_all
)

core_triple_strong <- pattern_df |>
  dplyr::filter(structure_grade == "three-compartment continuity: strong evidence") |>
  dplyr::arrange(q_RI, q_RC, q_IC, dplyr::desc(abs_r_RI + abs_r_RC + abs_r_IC))

core_triple_main <- pattern_df |>
  dplyr::filter(structure_grade == "three-compartment continuity: main evidence") |>
  dplyr::arrange(q_RI, q_RC, q_IC, dplyr::desc(abs_r_RI + abs_r_RC + abs_r_IC))

core_triple_support <- pattern_df |>
  dplyr::filter(structure_grade == "three-compartment continuity: supporting evidence") |>
  dplyr::arrange(q_RI, q_RC, q_IC, dplyr::desc(abs_r_RI + abs_r_RC + abs_r_IC))

############################
## 12. Summary statistics
############################
summary_pairs <- dplyr::bind_rows(
  make_pair_summary(res_RI, "Rumen_vs_Ileum"),
  make_pair_summary(res_RC, "Rumen_vs_Colon"),
  make_pair_summary(res_IC, "Ileum_vs_Colon")
)

summary_structure_grade <- pattern_df |>
  dplyr::count(structure_grade, sort = TRUE, name = "n")

summary_pair_signature <- pattern_df |>
  dplyr::count(pair_signature, sort = TRUE, name = "n")

############################
## 13. Export CLR matrices for the three categories of three-compartment core genera
############################
triple_strong_genus  <- core_triple_strong$Genus
triple_main_genus    <- core_triple_main$Genus
triple_support_genus <- core_triple_support$Genus

rum_triple_strong_clr  <- extract_clr_df_by_genus(rum_clr_df, triple_strong_genus)
ile_triple_strong_clr  <- extract_clr_df_by_genus(ile_clr_df, triple_strong_genus)
col_triple_strong_clr  <- extract_clr_df_by_genus(col_clr_df, triple_strong_genus)

rum_triple_main_clr    <- extract_clr_df_by_genus(rum_clr_df, triple_main_genus)
ile_triple_main_clr    <- extract_clr_df_by_genus(ile_clr_df, triple_main_genus)
col_triple_main_clr    <- extract_clr_df_by_genus(col_clr_df, triple_main_genus)

rum_triple_support_clr <- extract_clr_df_by_genus(rum_clr_df, triple_support_genus)
ile_triple_support_clr <- extract_clr_df_by_genus(ile_clr_df, triple_support_genus)
col_triple_support_clr <- extract_clr_df_by_genus(col_clr_df, triple_support_genus)

############################
## 14. Write process file
############################
cat("======================================\n")
cat("Writing process file...\n")
cat("======================================\n")

process_list <- list(
  Rumen_raw_count = data.frame(Genus = rownames(rum_raw_mat), rum_raw_mat, check.names = FALSE),
  Ileum_raw_count = data.frame(Genus = rownames(ile_raw_mat), ile_raw_mat, check.names = FALSE),
  Colon_raw_count = data.frame(Genus = rownames(col_raw_mat), col_raw_mat, check.names = FALSE),
  
  Rumen_all_stats = rum_stats_all,
  Ileum_all_stats = ile_stats_all,
  Colon_all_stats = col_stats_all,
  
  Rumen_single_site_core_info = rum_core$info,
  Ileum_single_site_core_info = ile_core$info,
  Colon_single_site_core_info = col_core$info,
  
  Rumen_core_count = data.frame(Genus = rownames(rum_raw_core), rum_raw_core, check.names = FALSE),
  Ileum_core_count = data.frame(Genus = rownames(ile_raw_core), ile_raw_core, check.names = FALSE),
  Colon_core_count = data.frame(Genus = rownames(col_raw_core), col_raw_core, check.names = FALSE),
  
  Rumen_CLR = rum_clr_df,
  Ileum_CLR = ile_clr_df,
  Colon_CLR = col_clr_df,
  
  Entry_Rumen_vs_Ileum = entry_RI,
  Entry_Rumen_vs_Colon = entry_RC,
  Entry_Ileum_vs_Colon = entry_IC
)

openxlsx::write.xlsx(
  process_list,
  file = process_xlsx,
  overwrite = TRUE
)

############################
## 15. Write result file
############################
cat("======================================\n")
cat("Writing result file...\n")
cat("======================================\n")

result_list <- list(
  Pair_summary = summary_pairs,
  
  Rumen_vs_Ileum_all = res_RI,
  Rumen_vs_Colon_all = res_RC,
  Ileum_vs_Colon_all = res_IC,
  
  Rumen_vs_Ileum_main_sig = sig_RI_main,
  Rumen_vs_Ileum_strict = sig_RI_strict,
  Rumen_vs_Ileum_same_direction = RI_same,
  Rumen_vs_Ileum_replacement = RI_replace,
  
  Rumen_vs_Colon_main_sig = sig_RC_main,
  Rumen_vs_Colon_strict = sig_RC_strict,
  Rumen_vs_Colon_same_direction = RC_same,
  Rumen_vs_Colon_replacement = RC_replace,
  
  Ileum_vs_Colon_main_sig = sig_IC_main,
  Ileum_vs_Colon_strict = sig_IC_strict,
  Ileum_vs_Colon_same_direction = IC_same,
  Ileum_vs_Colon_replacement = IC_replace,
  
  Pattern_table = pattern_df,
  Structure_grade_summary = summary_structure_grade,
  Pair_signature_summary = summary_pair_signature,
  
  Triple_strong = core_triple_strong,
  Triple_main = core_triple_main,
  Triple_support = core_triple_support,
  
  Rumen_triple_strong_CLR = rum_triple_strong_clr,
  Ileum_triple_strong_CLR = ile_triple_strong_clr,
  Colon_triple_strong_CLR = col_triple_strong_clr,
  
  Rumen_triple_main_CLR = rum_triple_main_clr,
  Ileum_triple_main_CLR = ile_triple_main_clr,
  Colon_triple_main_CLR = col_triple_main_clr,
  
  Rumen_triple_support_CLR = rum_triple_support_clr,
  Ileum_triple_support_CLR = ile_triple_support_clr,
  Colon_triple_support_CLR = col_triple_support_clr
)

openxlsx::write.xlsx(
  result_list,
  file = result_xlsx,
  overwrite = TRUE
)

############################
## 16. Console summary
############################
cat("======================================\n")
cat("Analysis completed.\n")
cat("Process file: ", process_xlsx, "\n")
cat("Result file: ", result_xlsx, "\n")
cat("======================================\n")

cat("Number of compartment-specific core genera (prevalence >= ", prevalence_core_cut, "）：\n", sep = "")
cat("Rumen :", nrow(rum_raw_core), "\n")
cat("Ileum :", nrow(ile_raw_core), "\n")
cat("Colon :", nrow(col_raw_core), "\n\n")

cat("Number of genera included in the three paired analyses:\n")
cat("Rumen_vs_Ileum :", sum(entry_RI$pass_pair_entry, na.rm = TRUE), "\n")
cat("Rumen_vs_Colon :", sum(entry_RC$pass_pair_entry, na.rm = TRUE), "\n")
cat("Ileum_vs_Colon :", sum(entry_IC$pass_pair_entry, na.rm = TRUE), "\n\n")

cat("Number of genera significant in the Pearson main analysis (q < ", q_cut_main, "）：\n", sep = "")
cat("Rumen_vs_Ileum :", nrow(sig_RI_main), "\n")
cat("Rumen_vs_Colon :", nrow(sig_RC_main), "\n")
cat("Ileum_vs_Colon :", nrow(sig_IC_main), "\n\n")

cat("Number of strictly significant genera (q < ", q_cut_main, " and |r| >= ", abs_r_cut_strict, "）：\n", sep = "")
cat("Rumen_vs_Ileum :", nrow(sig_RI_strict), "\n")
cat("Rumen_vs_Colon :", nrow(sig_RC_strict), "\n")
cat("Ileum_vs_Colon :", nrow(sig_IC_strict), "\n\n")

cat("Three-compartment structural grades:\n")
print(summary_structure_grade)
cat("======================================\n")