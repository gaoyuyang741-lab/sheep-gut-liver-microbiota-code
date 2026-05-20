# ============================================================
# 5-block unsupervised sGCCA for rumen 16S + ileum 16S + colon 16S
# + liver transcriptome + 1 phenotype
#
# Traits run separately:
#   GLU, TC, TG, TBA, LDL,
#   TailFat_g, TailFat_Carcass_g_per_kg, TailFat_PreLive_g_per_kg
#
# Main preprocessing strategy:
#   · Rumen / Ileum / Colon microbiome:
#       abundance(count) -> prevalence filter
#       -> CZM zero replacement -> CLR -> z-score
#   · Liver transcriptome: use formal_log2TPM sheet directly -> z-score
#   · Phenotype: single trait block -> z-score
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  BiocManager::install("mixOmics", ask = FALSE, update = FALSE)
}
if (!requireNamespace("zCompositions", quietly = TRUE)) {
  install.packages("zCompositions", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
if (!requireNamespace("compositions", quietly = TRUE)) {
  install.packages("compositions", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(mixOmics)
  library(zCompositions)
  library(compositions)
})

# ---- avoid namespace masking ----
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
bind_rows <- dplyr::bind_rows
full_join <- dplyr::full_join
rename    <- dplyr::rename

# ============================================================
# Parameters
# ============================================================
# Input and output directories
# Please place the required sGCCA input files in "data/sgcca/input".
# Output files will be saved in "results/sgcca/sgcca_5block_8traits_out".

workdir <- file.path("data", "sgcca", "input")

file_rum     <- file.path(workdir, "Rum_genus_abundance.xlsx")
file_ile     <- file.path(workdir, "Ile_genus_abundance.xlsx")
file_col     <- file.path(workdir, "Col_genus_abundance.xlsx")
file_liver   <- file.path(workdir, "肝脏筛选基因.xlsx")
file_blood   <- file.path(workdir, "blood_phenotype_sGCCA_input.xlsx")
file_tailfat <- file.path(workdir, "尾脂数据_16s编号转换.xlsx")

liver_sheet <- "formal_log2TPM"
blood_sheet <- "sgcca_blood_input"

traits_to_run <- c(
  "GLU",
  "TC",
  "TG",
  "TBA",
  "LDL",
  "TailFat_g",
  "TailFat_Carcass_g_per_kg",
  "TailFat_PreLive_g_per_kg"
)

tailfat_traits <- c(
  "TailFat_g",
  "TailFat_Carcass_g_per_kg",
  "TailFat_PreLive_g_per_kg"
)

# Microbiome preprocessing
prevalence_cut_microbe <- 0.50

# sGCCA parameters
ncomp          <- 1
seed           <- 123
sparsity_level <- 0.30
network_cutoff <- 0.30
topN           <- 30

outdir <- file.path("results", "sgcca", "sgcca_5block_8traits_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# ============================================================
# Helpers
# ============================================================
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, useBytes = TRUE)
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("-", "_", x)
  x <- gsub("[()（）\\[\\]]", "", x)
  x
}

safe_numeric_df <- function(df) {
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  df
}

read_microbe_first_sheet <- function(path, prevalence_cut = 0.10) {
  df <- readxl::read_excel(path, sheet = 1, col_names = TRUE)
  df <- df[, !grepl("^Unnamed", colnames(df)), drop = FALSE]
  
  all_na_cols <- vapply(df, function(col) all(is.na(col)), logical(1))
  if (any(all_na_cols)) df <- df[, !all_na_cols, drop = FALSE]
  
  feat_col <- colnames(df)[1]
  feats <- make.names(as.character(df[[feat_col]]), unique = TRUE)
  
  mat_df <- df[, -1, drop = FALSE]
  colnames(mat_df) <- normalize_ids(colnames(mat_df))
  mat_df <- safe_numeric_df(mat_df)
  
  X <- t(as.matrix(mat_df))
  rownames(X) <- colnames(mat_df)
  colnames(X) <- feats
  
  keep_var <- apply(X, 2, function(v) {
    vv <- var(v, na.rm = TRUE)
    is.finite(vv) && vv > 0
  })
  X <- X[, keep_var, drop = FALSE]
  
  prev <- colMeans(X > 0, na.rm = TRUE)
  keep_prev <- prev >= prevalence_cut
  X <- X[, keep_prev, drop = FALSE]
  
  if (ncol(X) == 0) {
    stop(sprintf("No microbiome feature remains after prevalence filtering: %s", path))
  }
  X
}

read_liver_strict_log2tpm <- function(path, sheet = "formal_log2TPM") {
  df <- readxl::read_excel(path, sheet = sheet, skip = 1, col_names = TRUE)
  df <- df[, !grepl("^Unnamed", colnames(df)), drop = FALSE]
  
  genes <- make.names(as.character(df[[1]]), unique = TRUE)
  
  mat_df <- df[, -1, drop = FALSE]
  colnames(mat_df) <- normalize_ids(colnames(mat_df))
  mat_df <- safe_numeric_df(mat_df)
  
  X <- t(as.matrix(mat_df))
  rownames(X) <- colnames(mat_df)
  colnames(X) <- genes
  
  keep <- apply(X, 2, function(v) {
    vv <- var(v, na.rm = TRUE)
    is.finite(vv) && vv > 0
  })
  X <- X[, keep, drop = FALSE]
  
  if (ncol(X) == 0) stop("No liver feature remains after filtering.")
  X
}

read_blood_trait_block <- function(path, trait, sheet = "sgcca_blood_input") {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE)
  
  feat_row <- which(as.character(raw[[1]]) == "Feature")
  if (length(feat_row) == 0) {
    stop("Could not find the 'Feature' row in blood input file.")
  }
  feat_row <- feat_row[1]
  
  sample_ids <- unlist(raw[feat_row, -1], use.names = FALSE)
  sample_ids <- normalize_ids(sample_ids)
  
  dat <- raw[(feat_row + 1):nrow(raw), , drop = FALSE]
  colnames(dat) <- c("Feature", sample_ids)
  dat$Feature <- as.character(dat$Feature)
  
  if (!(trait %in% dat$Feature)) {
    stop(sprintf("Trait '%s' not found in blood input sheet.", trait))
  }
  
  one <- dat[dat$Feature == trait, , drop = FALSE]
  vals <- suppressWarnings(as.numeric(one[1, -1]))
  
  X <- matrix(vals, ncol = 1)
  rownames(X) <- sample_ids
  colnames(X) <- trait
  X
}

read_tailfat_trait_block <- function(path, trait, sheet = "尾脂_16s编号") {
  # 前4行是说明，第5行才是真正表头
  df <- readxl::read_excel(path, sheet = sheet, skip = 4, col_names = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- df[, !grepl("^Unnamed", colnames(df)), drop = FALSE]
  
  need_cols <- c("瘤胃和肠道16s编号", "尾脂g", "尾脂/胴体重g/kg", "尾脂/宰前活重g/kg")
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop(sprintf(
      "Tail fat file missing columns: %s",
      paste(miss_cols, collapse = ", ")
    ))
  }
  
  col_map <- c(
    "TailFat_g"                = "尾脂g",
    "TailFat_Carcass_g_per_kg" = "尾脂/胴体重g/kg",
    "TailFat_PreLive_g_per_kg" = "尾脂/宰前活重g/kg"
  )
  
  src_col <- col_map[trait]
  if (is.na(src_col)) {
    stop(sprintf("Unknown tail fat trait: %s", trait))
  }
  
  df <- df[, c("瘤胃和肠道16s编号", src_col), drop = FALSE]
  colnames(df) <- c("SampleID", trait)
  
  df$SampleID <- normalize_ids(df$SampleID)
  df[[trait]] <- suppressWarnings(as.numeric(df[[trait]]))
  df <- df[!is.na(df$SampleID) & !is.na(df[[trait]]), , drop = FALSE]
  
  X <- as.matrix(df[, trait, drop = FALSE])
  rownames(X) <- df$SampleID
  colnames(X) <- trait
  X
}

microbe_preprocess <- function(X) {
  X <- as.matrix(X)
  mode(X) <- "numeric"
  
  rs <- rowSums(X, na.rm = TRUE)
  if (any(!is.finite(rs) | rs <= 0)) {
    stop("Microbiome block contains sample(s) with non-positive total abundance.")
  }
  
  X_nozero <- zCompositions::cmultRepl(
    X,
    method = "CZM",
    output = "p-counts"
  )
  
  X_clr <- compositions::clr(X_nozero, base = exp(1))
  X_clr <- as.matrix(X_clr)
  
  X_z <- scale(X_clr)
  keep <- apply(X_z, 2, function(v) all(is.finite(v)))
  X_z[, keep, drop = FALSE]
}

z_preprocess <- function(X) {
  X <- as.matrix(X)
  Xs <- scale(X)
  keep <- apply(Xs, 2, function(v) all(is.finite(v)))
  Xs[, keep, drop = FALSE]
}

auto_penalty <- function(p, sparsity = 0.3) {
  min_c1 <- 1 / sqrt(p)
  max_c1 <- 1
  min_c1 + sparsity * (max_c1 - min_c1)
}

safe_get_block_scores <- function(mb, blocks, ncomp = 1) {
  cand <- c("Y", "variates", "scores")
  for (slot in cand) {
    if (!is.null(mb[[slot]]) && is.list(mb[[slot]])) {
      out <- lapply(mb[[slot]], function(M) {
        M <- as.matrix(M)
        M <- M[, seq_len(min(ncomp, ncol(M))), drop = FALSE]
        colnames(M) <- paste0("comp", seq_len(ncol(M)))
        M
      })
      return(out)
    }
  }
  
  if (is.null(mb$a)) stop("No score-related slot found and no loadings slot available.")
  
  out <- list()
  for (b in names(blocks)) {
    A <- mb$a[[b]]
    if (is.null(A)) next
    X <- as.matrix(blocks[[b]])
    Xcs <- scale(X, center = TRUE, scale = TRUE)
    p <- min(ncol(Xcs), nrow(A))
    S <- Xcs[, seq_len(p), drop = FALSE] %*%
      as.matrix(A[seq_len(p), seq_len(min(ncomp, ncol(A))), drop = FALSE])
    rownames(S) <- rownames(X)
    colnames(S) <- paste0("comp", seq_len(ncol(S)))
    out[[b]] <- S
  }
  out
}

safe_get_block_loadings <- function(mb, blocks, scores_blocks, ncomp = 1) {
  if (!is.null(mb$a)) return(mb$a)
  if (!is.null(mb$astar)) return(mb$astar)
  
  out <- list()
  for (b in names(blocks)) {
    X <- as.matrix(blocks[[b]])
    if (is.null(scores_blocks[[b]])) next
    S <- as.matrix(scores_blocks[[b]])[, seq_len(min(ncomp, ncol(scores_blocks[[b]]))), drop = FALSE]
    Xs <- scale(X, center = TRUE, scale = TRUE)
    L <- suppressWarnings(cor(Xs, S, use = "pairwise.complete.obs"))
    L <- as.matrix(L)
    rownames(L) <- colnames(X)
    colnames(L) <- paste0("comp", seq_len(ncol(L)))
    out[[b]] <- L
  }
  if (length(out) == 0) stop("Failed to obtain block loadings.")
  out
}

build_network_edges <- function(blocks, loadings_list, comp = 1, cutoff = 0.3) {
  blks <- names(blocks)
  out  <- list()
  
  sel_feats <- lapply(blks, function(b) {
    Lb <- as.matrix(loadings_list[[b]])
    if (is.null(Lb)) character(0) else rownames(Lb)[abs(Lb[, comp, drop = TRUE]) > 0]
  })
  names(sel_feats) <- blks
  
  if (length(blks) < 2) return(data.frame())
  
  for (i in 1:(length(blks) - 1)) {
    for (j in (i + 1):length(blks)) {
      b1 <- blks[i]
      b2 <- blks[j]
      f1 <- sel_feats[[b1]]
      f2 <- sel_feats[[b2]]
      
      if (length(f1) == 0 || length(f2) == 0) next
      
      X <- as.matrix(blocks[[b1]][, f1, drop = FALSE])
      Y <- as.matrix(blocks[[b2]][, f2, drop = FALSE])
      
      C <- suppressWarnings(cor(X, Y, use = "pairwise.complete.obs"))
      if (!is.matrix(C)) next
      
      df <- as.data.frame(as.table(C))
      colnames(df) <- c("var1", "var2", "correlation")
      
      df <- df %>%
        filter(is.finite(correlation)) %>%
        mutate(
          block1   = b1,
          block2   = b2,
          abs_corr = abs(correlation)
        ) %>%
        filter(abs_corr >= cutoff) %>%
        select(block1, var1, block2, var2, correlation, abs_corr) %>%
        arrange(desc(abs_corr))
      
      if (nrow(df) > 0) out[[paste(b1, b2, sep = "_")]] <- df
    }
  }
  
  if (length(out) == 0) data.frame() else bind_rows(out, .id = "pair")
}

write_loading_outputs <- function(loadings_list, outdir_trait, topN = 30) {
  all_nonzero <- list()
  top_tbls    <- list()
  
  for (b in names(loadings_list)) {
    L <- as.matrix(loadings_list[[b]])
    if (is.null(L) || nrow(L) == 0) next
    if (!"comp1" %in% colnames(L)) colnames(L)[1] <- "comp1"
    
    df <- data.frame(
      Block     = b,
      Feature   = rownames(L),
      comp1     = L[, "comp1", drop = TRUE],
      abs_comp1 = abs(L[, "comp1", drop = TRUE]),
      stringsAsFactors = FALSE
    ) %>% arrange(desc(abs_comp1))
    
    df_nonzero <- df %>% filter(abs_comp1 > 0)
    df_top     <- head(df_nonzero, topN)
    
    all_nonzero[[b]] <- df_nonzero
    top_tbls[[b]]    <- df_top
    
    if (nrow(df_top) > 0) {
      p <- ggplot(df_top, aes(x = reorder(Feature, comp1), y = comp1)) +
        geom_col() +
        coord_flip() +
        labs(
          title = paste0(b, " top loadings (comp1)"),
          x = NULL,
          y = "Loading"
        ) +
        theme_bw(base_size = 11)
      
      ggsave(
        file.path(outdir_trait, paste0("top_loadings_", b, ".pdf")),
        p, width = 7, height = 6
      )
    }
  }
  
  readr::write_csv(bind_rows(all_nonzero), file.path(outdir_trait, "loadings_nonzero_all_blocks.csv"))
  readr::write_csv(bind_rows(top_tbls), file.path(outdir_trait, "loadings_topN_all_blocks.csv"))
}

# ============================================================
# Read fixed blocks once
# ============================================================
message(">> Reading fixed blocks ...")
rum_raw   <- read_microbe_first_sheet(file_rum, prevalence_cut = prevalence_cut_microbe)
ile_raw   <- read_microbe_first_sheet(file_ile, prevalence_cut = prevalence_cut_microbe)
col_raw   <- read_microbe_first_sheet(file_col, prevalence_cut = prevalence_cut_microbe)
liver_raw <- read_liver_strict_log2tpm(file_liver, sheet = liver_sheet)

message(">> Dimensions before per-trait alignment:")
print(list(
  Rum   = dim(rum_raw),
  Ile   = dim(ile_raw),
  Col   = dim(col_raw),
  Liver = dim(liver_raw)
))

# ============================================================
# Run trait by trait
# ============================================================
set.seed(seed)
summary_list <- list()

for (trait in traits_to_run) {
  message("\n==============================")
  message(">> Running sGCCA for trait: ", trait)
  message("==============================")
  
  trait_dir <- file.path(outdir, paste0("sGCCA_", trait))
  dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (trait %in% tailfat_traits) {
    pheno_raw <- read_tailfat_trait_block(file_tailfat, trait = trait, sheet = "尾脂_16s编号")
  } else {
    pheno_raw <- read_blood_trait_block(file_blood, trait = trait, sheet = blood_sheet)
  }
  
  common_ids <- Reduce(intersect, list(
    rownames(rum_raw),
    rownames(ile_raw),
    rownames(col_raw),
    rownames(liver_raw),
    rownames(pheno_raw)
  ))
  common_ids <- sort(common_ids)
  
  if (length(common_ids) < 10) {
    stop(sprintf("Too few common samples for trait %s: %d", trait, length(common_ids)))
  }
  
  message(">> Common samples used: ", length(common_ids))
  
  blocks_raw <- list(
    Rum   = rum_raw[common_ids, , drop = FALSE],
    Ile   = ile_raw[common_ids, , drop = FALSE],
    Col   = col_raw[common_ids, , drop = FALSE],
    Liver = liver_raw[common_ids, , drop = FALSE],
    Trait = pheno_raw[common_ids, , drop = FALSE]
  )
  
  blocks <- list(
    Rum   = microbe_preprocess(blocks_raw$Rum),
    Ile   = microbe_preprocess(blocks_raw$Ile),
    Col   = microbe_preprocess(blocks_raw$Col),
    Liver = z_preprocess(blocks_raw$Liver),
    Trait = z_preprocess(blocks_raw$Trait)
  )
  
  blocks <- lapply(blocks, function(x) {
    x[, apply(x, 2, function(v) all(is.finite(v))), drop = FALSE]
  })
  
  if (any(vapply(blocks, ncol, integer(1)) == 0)) {
    stop(sprintf("At least one block has zero columns after preprocessing for trait %s.", trait))
  }
  
  readr::write_csv(
    data.frame(SampleID = common_ids),
    file.path(trait_dir, "aligned_samples.csv")
  )
  
  readr::write_csv(
    data.frame(
      Block    = names(blocks),
      nSample  = vapply(blocks, nrow, integer(1)),
      nFeature = vapply(blocks, ncol, integer(1))
    ),
    file.path(trait_dir, "block_dimensions.csv")
  )
  
  B <- length(blocks)
  design <- matrix(1, nrow = B, ncol = B, dimnames = list(names(blocks), names(blocks)))
  diag(design) <- 0
  
  penalty_vec <- sapply(blocks, function(m) auto_penalty(ncol(m), sparsity = sparsity_level))
  names(penalty_vec) <- names(blocks)
  
  message(">> penalty (c1):")
  print(penalty_vec)
  
  use_mixomics_sgcca <- exists("wrapper.sgcca", where = asNamespace("mixOmics"), inherits = FALSE)
  
  if (use_mixomics_sgcca) {
    mb <- mixOmics::wrapper.sgcca(
      X       = blocks,
      design  = design,
      penalty = penalty_vec,
      ncomp   = ncomp,
      scale   = TRUE
    )
  } else if (requireNamespace("RGCCA", quietly = TRUE)) {
    mb <- RGCCA::sgcca(
      A       = blocks,
      C       = design,
      c1      = penalty_vec,
      ncomp   = rep(ncomp, B),
      scale   = TRUE,
      verbose = TRUE
    )
  } else {
    stop("Neither mixOmics::wrapper.sgcca nor RGCCA::sgcca is available.")
  }
  
  saveRDS(mb, file.path(trait_dir, paste0("sgcca_model_", trait, ".rds")))
  
  diag_file <- file.path(trait_dir, "model_diagnostics.txt")
  sink(diag_file)
  cat("===== sGCCA Model Diagnostics =====\n\n")
  cat("Trait:\n")
  print(trait); cat("\n")
  cat("Penalty (c1) per block:\n")
  print(penalty_vec); cat("\n")
  if (!is.null(mb$AVE)) {
    cat("Average Variance Explained (AVE):\n")
    print(mb$AVE); cat("\n")
  } else {
    cat("No AVE slot in model.\n\n")
  }
  sink()
  
  scores_blocks   <- safe_get_block_scores(mb, blocks, ncomp = ncomp)
  loadings_blocks <- safe_get_block_loadings(mb, blocks, scores_blocks, ncomp = ncomp)
  
  scores_long <- bind_rows(lapply(names(scores_blocks), function(b) {
    S <- as.data.frame(scores_blocks[[b]])
    S$SampleID <- rownames(scores_blocks[[b]])
    S$Block <- b
    S
  }))
  readr::write_csv(scores_long, file.path(trait_dir, "block_scores.csv"))
  
  write_loading_outputs(loadings_blocks, trait_dir, topN = topN)
  
  score_mat <- do.call(cbind, lapply(names(scores_blocks), function(b) {
    x <- as.matrix(scores_blocks[[b]])
    colnames(x) <- paste0(b, "_", colnames(x))
    x
  }))
  
  block_cor <- suppressWarnings(cor(score_mat, use = "pairwise.complete.obs"))
  write.csv(block_cor, file.path(trait_dir, "block_score_correlations.csv"), row.names = TRUE)
  
  pdf(file.path(trait_dir, "block_score_correlations_heatmap.pdf"), width = 7, height = 6)
  pheatmap::pheatmap(block_cor, main = paste0("Score correlation heatmap - ", trait))
  dev.off()
  
  comp1_tbl <- Reduce(function(x, y) full_join(x, y, by = "SampleID"),
                      lapply(names(scores_blocks), function(b) {
                        data.frame(
                          SampleID = rownames(scores_blocks[[b]]),
                          score    = scores_blocks[[b]][, 1, drop = TRUE],
                          stringsAsFactors = FALSE
                        ) %>%
                          rename(!!paste0(b, "_comp1") := score)
                      }))
  readr::write_csv(comp1_tbl, file.path(trait_dir, "component1_scores_wide.csv"))
  
  edges <- build_network_edges(blocks, loadings_blocks, comp = 1, cutoff = network_cutoff)
  if (nrow(edges) > 0) {
    readr::write_csv(edges, file.path(trait_dir, "network_edges.csv"))
  } else {
    readr::write_csv(data.frame(), file.path(trait_dir, "network_edges.csv"))
  }
  
  summary_list[[trait]] <- data.frame(
    Trait         = trait,
    nSample       = length(common_ids),
    Rum_feature   = ncol(blocks$Rum),
    Ile_feature   = ncol(blocks$Ile),
    Col_feature   = ncol(blocks$Col),
    Liver_feature = ncol(blocks$Liver),
    Trait_feature = ncol(blocks$Trait),
    stringsAsFactors = FALSE
  )
}

readr::write_csv(bind_rows(summary_list), file.path(outdir, "sgcca_run_summary.csv"))
message("\nAll done. Outputs are saved in: ", normalizePath(outdir, winslash = "/", mustWork = FALSE))