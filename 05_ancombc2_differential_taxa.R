# =========================================================
# Three-site genus-level spatial stratification analysis script
# Design rationale:
# 1. Use the genus count table for ANCOM-BC2
# 2. Use the Friedman test as the overall paired three-group entry test
# 3. Use paired Wilcoxon tests as pairwise auxiliary tests
# 4. Classify spatial stratification patterns based on median relative abundance across the three sites
# =========================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

# =========================================================
# 0. Load packages
# =========================================================
need_cran <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tibble",
  "stringr",
  "purrr",
  "tidyr",
  "lme4"
)

for (pkg in need_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", ask = FALSE, update = FALSE)
}
if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  BiocManager::install("ANCOMBC", ask = FALSE, update = FALSE)
}

library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(tidyr)
library(lme4)
library(phyloseq)
library(ANCOMBC)

# =========================================================
# 1. Parameter settings
# =========================================================
# Input and output directories
# Please place the required input files in "data/feature_taxa/input".
# Output files will be saved in "results/feature_taxa/q0.001".

in_dir  <- file.path("data", "feature_taxa", "input")
out_dir <- file.path("results", "feature_taxa", "q0.001")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

abund_file <- file.path(in_dir, "genus_abundance_3group_merged.xlsx")
meta_file  <- file.path(in_dir, "metadata_3group.xlsx")
tax_file   <- file.path(in_dir, "taxonomy_genus.xlsx")

# Grouping and pairing
group_var    <- "Group"
subject_var  <- "SheepID"
group_levels <- c("Rum", "Ile", "Col")

# ANCOM-BC2 parameters
prv_cut      <- 0.10
lib_cut      <- 1000
pseudo_sens  <- TRUE
struc_zero   <- FALSE
neg_lb       <- FALSE
alpha_main   <- 0.001
p_adj_method <- "BH"

# Friedman / pairwise / pattern-classification thresholds
overall_q_cut  <- 0.001
pairwise_q_cut <- 0.001

# Number of representative genera to output
top_n_each_subclass <- 15

# =========================================================
# 2. Read data
# =========================================================
abund_df <- read_excel(abund_file, sheet = 1)
meta_df  <- read_excel(meta_file, sheet = 1)
tax_df   <- read_excel(tax_file, sheet = 1)

colnames(abund_df)[1] <- "Genus"
colnames(meta_df) <- c("SampleID", "Group", "SheepID")
colnames(tax_df)[1] <- "Genus"

abund_df <- abund_df %>%
  filter(!is.na(Genus), Genus != "")

# Deduplicate genera in the abundance table: sum rows if duplicated
if (anyDuplicated(abund_df$Genus) > 0) {
  abund_df <- abund_df %>%
    group_by(Genus) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Deduplicate genera in taxonomy: keep the first record
if (anyDuplicated(tax_df$Genus) > 0) {
  tax_df <- tax_df %>%
    distinct(Genus, .keep_all = TRUE)
}

meta_df <- meta_df %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = factor(Group, levels = group_levels),
    SheepID = as.character(SheepID)
  )

# =========================================================
# 3. Input consistency check
# =========================================================
sample_cols <- colnames(abund_df)[-1]

missing_in_meta  <- setdiff(sample_cols, meta_df$SampleID)
missing_in_abund <- setdiff(meta_df$SampleID, sample_cols)

if (length(missing_in_meta) > 0) {
  stop("The following samples in the abundance table are absent from metadata:\n", paste(missing_in_meta, collapse = ", "))
}
if (length(missing_in_abund) > 0) {
  stop("The following samples in metadata are absent from the abundance table:\n", paste(missing_in_abund, collapse = ", "))
}

# Reorder abundance table according to metadata
abund_df <- abund_df %>%
  select(Genus, all_of(meta_df$SampleID))

group_count <- table(meta_df$Group)
sheep_count <- table(meta_df$SheepID)

qc_sample_alignment <- data.frame(
  check_item = c(
    "Abundance-table sample count",
    "metadata sample count",
    "Shared sample count",
    "Rum sample count",
    "Ile sample count",
    "Col sample count",
    "SheepID count",
    "Whether each SheepID has exactly three samples"
  ),
  value = c(
    length(sample_cols),
    nrow(meta_df),
    length(intersect(sample_cols, meta_df$SampleID)),
    unname(group_count["Rum"]),
    unname(group_count["Ile"]),
    unname(group_count["Col"]),
    length(unique(meta_df$SheepID)),
    ifelse(all(sheep_count == 3), "Yes", "No")
  )
)

if (!all(sheep_count == 3)) {
  warning("Not all SheepID values correspond to exactly three site samples; please check the paired structure.")
}

# =========================================================
# 4. Construct phyloseq object
# =========================================================
otu_mat <- abund_df %>%
  column_to_rownames("Genus") %>%
  as.matrix()

storage.mode(otu_mat) <- "numeric"

meta_use <- meta_df %>%
  as.data.frame()
rownames(meta_use) <- meta_use$SampleID

tax_use <- tax_df %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "Unclassified", as.character(.x)))) %>%
  right_join(data.frame(Genus = rownames(otu_mat)), by = "Genus") %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "Unclassified", .x))) %>%
  distinct(Genus, .keep_all = TRUE)

tax_mat <- tax_use %>%
  column_to_rownames("Genus") %>%
  as.matrix()

tax_mat <- tax_mat[rownames(otu_mat), , drop = FALSE]

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(meta_use),
  tax_table(tax_mat)
)

# =========================================================
# 5. Preprocessing statistics
# =========================================================
lib_sizes <- sample_sums(ps)

prev_df <- data.frame(
  Genus = taxa_names(ps),
  prevalence_all = apply(otu_table(ps), 1, function(x) mean(x > 0)),
  prevalence_Rum = apply(otu_table(ps)[, sample_data(ps)$Group == "Rum"], 1, function(x) mean(x > 0)),
  prevalence_Ile = apply(otu_table(ps)[, sample_data(ps)$Group == "Ile"], 1, function(x) mean(x > 0)),
  prevalence_Col = apply(otu_table(ps)[, sample_data(ps)$Group == "Col"], 1, function(x) mean(x > 0))
)

qc_preprocess <- data.frame(
  item = c(
    "Raw genus count",
    "sample count",
    "Minimum sequencing depth",
    "Median sequencing depth",
    "Maximum sequencing depth",
    "Number of genera with overall prevalence >= 0.10"
  ),
  value = c(
    ntaxa(ps),
    nsamples(ps),
    min(lib_sizes),
    median(lib_sizes),
    max(lib_sizes),
    sum(prev_df$prevalence_all >= prv_cut)
  )
)

# =========================================================
# 6. Run ANCOM-BC2
#    ANCOM-BC2 is retained here, but res_global is no longer used as the overall entry test
# =========================================================
message("Starting ANCOM-BC2 ...")

set.seed(123)

res_ancom <- ancombc2(
  data = ps,
  tax_level = NULL,
  fix_formula = group_var,
  rand_formula = paste0("(1|", subject_var, ")"),
  p_adj_method = p_adj_method,
  pseudo_sens = pseudo_sens,
  prv_cut = prv_cut,
  lib_cut = lib_cut,
  group = group_var,
  struc_zero = struc_zero,
  neg_lb = neg_lb,
  alpha = alpha_main,
  global = TRUE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  ),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  verbose = TRUE
)

message("ANCOM-BC2 completed.")

# =========================================================
# 7. Format main ANCOM-BC2 results
# =========================================================
res_main <- as.data.frame(res_ancom$res)

if ("taxon" %in% colnames(res_main)) {
  res_main$Genus <- res_main$taxon
} else {
  res_main <- res_main %>% rownames_to_column("Genus")
}

res_main <- res_main %>%
  select(Genus, everything())

ancom_main_tbl <- res_main %>%
  transmute(
    Genus,
    lfc_GroupIle,
    lfc_GroupCol,
    p_GroupIle,
    p_GroupCol,
    q_GroupIle,
    q_GroupCol,
    diff_GroupIle,
    diff_GroupCol,
    passed_ss_GroupIle,
    passed_ss_GroupCol,
    diff_robust_GroupIle,
    diff_robust_GroupCol
  ) %>%
  left_join(prev_df, by = "Genus") %>%
  left_join(tax_df, by = "Genus")

# Save raw ANCOM-BC2 results separately
ancom_raw_tbl <- res_main %>%
  left_join(prev_df, by = "Genus") %>%
  left_join(tax_df, by = "Genus")

# =========================================================
# 8. Calculate median relative abundance across the three sites (for pattern classification only)
# =========================================================
rel_mat <- sweep(otu_mat, 2, colSums(otu_mat), "/")
rel_mat[is.na(rel_mat)] <- 0

get_group_median <- function(mat, meta, grp) {
  idx <- meta$SampleID[meta$Group == grp]
  apply(mat[, idx, drop = FALSE], 1, median, na.rm = TRUE)
}

med_R <- get_group_median(rel_mat, meta_df, "Rum")
med_I <- get_group_median(rel_mat, meta_df, "Ile")
med_C <- get_group_median(rel_mat, meta_df, "Col")

make_order_string <- function(r, i, c) {
  vals <- c(Rum = r, Ile = i, Col = c)
  paste(names(sort(vals, decreasing = TRUE)), collapse = " > ")
}

pattern_base_tbl <- data.frame(
  Genus = rownames(rel_mat),
  med_Rum = med_R,
  med_Ile = med_I,
  med_Col = med_C,
  stringsAsFactors = FALSE
)

pattern_base_tbl$order_string <- pmap_chr(
  list(pattern_base_tbl$med_Rum, pattern_base_tbl$med_Ile, pattern_base_tbl$med_Col),
  make_order_string
)

# =========================================================
# 9. Construct long and wide tables
# =========================================================
rel_long <- as.data.frame(rel_mat) %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "RelAbund") %>%
  left_join(meta_df, by = "SampleID")

rel_wide <- rel_long %>%
  select(Genus, SheepID, Group, RelAbund) %>%
  pivot_wider(names_from = Group, values_from = RelAbund)

# =========================================================
# 10. Friedman overall test
#    This is the actual overall entry test for the paired three-site design
# =========================================================
message("Starting the Friedman overall test ...")

safe_friedman <- function(df_one_genus) {
  df_one_genus <- df_one_genus %>%
    select(SheepID, Rum, Ile, Col) %>%
    arrange(SheepID)
  
  mat <- as.matrix(df_one_genus[, c("Rum", "Ile", "Col")])
  
  if (nrow(mat) < 3) return(NA_real_)
  if (all(is.na(mat))) return(NA_real_)
  
  out <- tryCatch(
    friedman.test(mat)$p.value,
    error = function(e) NA_real_
  )
  out
}

friedman_tbl <- map_dfr(unique(rel_wide$Genus), function(g) {
  sub <- rel_wide %>% filter(Genus == g)
  data.frame(
    Genus = g,
    p_friedman = safe_friedman(sub)
  )
})

friedman_tbl <- friedman_tbl %>%
  mutate(
    q_friedman = p.adjust(p_friedman, method = "BH"),
    overall_sig_friedman = ifelse(!is.na(q_friedman) & q_friedman <= overall_q_cut, TRUE, FALSE)
  )

overall_tbl <- friedman_tbl %>%
  left_join(pattern_base_tbl %>% select(Genus, med_Rum, med_Ile, med_Col, order_string), by = "Genus") %>%
  left_join(ancom_main_tbl, by = "Genus") %>%
  arrange(q_friedman, p_friedman)

overall_sig_tbl <- overall_tbl %>%
  filter(overall_sig_friedman) %>%
  arrange(q_friedman, p_friedman)

# =========================================================
# 11. Run paired Wilcoxon tests for genera significant in the overall test
# =========================================================
message("Starting paired Wilcoxon tests for genera significant in the overall test ...")

safe_paired_wilcox <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) return(NA_real_)
  if (all(x == y)) return(1)
  
  out <- tryCatch(
    wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  out
}

if (nrow(overall_sig_tbl) == 0) {
  message("No Friedman overall significant genera were found at the current threshold; pairwise Wilcoxon tests were skipped.")
  pairwise_tbl <- data.frame(
    Genus = character(0),
    p_RI = numeric(0),
    p_RC = numeric(0),
    p_IC = numeric(0),
    q_RI = numeric(0),
    q_RC = numeric(0),
    q_IC = numeric(0),
    sig_RI = logical(0),
    sig_RC = logical(0),
    sig_IC = logical(0)
  )
} else {
  pairwise_tbl <- map_dfr(overall_sig_tbl$Genus, function(g) {
    sub <- rel_wide %>% filter(Genus == g)
    
    p_RI <- safe_paired_wilcox(sub$Rum, sub$Ile)
    p_RC <- safe_paired_wilcox(sub$Rum, sub$Col)
    p_IC <- safe_paired_wilcox(sub$Ile, sub$Col)
    
    data.frame(
      Genus = g,
      p_RI = p_RI,
      p_RC = p_RC,
      p_IC = p_IC
    )
  })
  
  pairwise_tbl <- pairwise_tbl %>%
    mutate(
      q_RI = p.adjust(p_RI, method = "BH"),
      q_RC = p.adjust(p_RC, method = "BH"),
      q_IC = p.adjust(p_IC, method = "BH"),
      sig_RI = ifelse(!is.na(q_RI) & q_RI <= pairwise_q_cut, TRUE, FALSE),
      sig_RC = ifelse(!is.na(q_RC) & q_RC <= pairwise_q_cut, TRUE, FALSE),
      sig_IC = ifelse(!is.na(q_IC) & q_IC <= pairwise_q_cut, TRUE, FALSE)
    )
}

# =========================================================
# 12. Spatial stratification pattern-classification function
# =========================================================
classify_pattern <- function(mR, mI, mC, sig_RI, sig_RC, sig_IC) {
  
  dir_string <- paste(names(sort(c(Rum = mR, Ile = mI, Col = mC), decreasing = TRUE)), collapse = " > ")
  
  # Foregut-hindgut coordinated pattern: Rum ≈ Col, higher or lower than Ile
  if (!sig_RC && sig_RI && sig_IC) {
    if (mR > mI && mC > mI) {
      return(c("Foregut-hindgut coordinated pattern", "Foregut-hindgut coordinated high pattern", "Rum ≈ Col > Ile"))
    }
    if (mR < mI && mC < mI) {
      return(c("Foregut-hindgut coordinated pattern", "Foregut-hindgut coordinated low pattern", "Rum ≈ Col < Ile"))
    }
  }
  
  # Ileum-only high / low
  if (sig_RI && sig_IC) {
    if (mI > mR && mI > mC) {
      return(c("Small-intestine transitional pattern", "Ileum-only high pattern", dir_string))
    }
    if (mI < mR && mI < mC) {
      return(c("Small-intestine transitional pattern", "Ileum-only low pattern", dir_string))
    }
  }
  
  # Decreasing gradient pattern: Rum > Ile > Col
  if (mR > mI && mI > mC) {
    if (sig_RC && (sig_RI || sig_IC)) {
      return(c("Gradient pattern", "Decreasing gradient pattern", "Rum > Ile > Col"))
    }
  }
  
  # Increasing gradient pattern: Col > Ile > Rum
  if (mC > mI && mI > mR) {
    if (sig_RC && (sig_RI || sig_IC)) {
      return(c("Gradient pattern", "Increasing gradient pattern", "Col > Ile > Rum"))
    }
  }
  
  # Foregut-enriched pattern
  if (mR > mI && mR > mC && sig_RI && sig_RC) {
    if (!sig_IC) {
      return(c("Foregut-enriched pattern", "Foregut-high with similar remaining sites", "Rum > Ile ≈ Col"))
    } else {
      return(c("Foregut-enriched pattern", "Foregut-high stratified pattern", dir_string))
    }
  }
  
  # Hindgut-enriched pattern
  if (mC > mR && mC > mI && sig_RC && sig_IC) {
    if (!sig_RI) {
      return(c("Hindgut-enriched pattern", "Hindgut-high with similar remaining sites", "Col > Ile ≈ Rum"))
    } else {
      return(c("Hindgut-enriched pattern", "Hindgut-high stratified pattern", dir_string))
    }
  }
  
  # Complex pattern
  return(c("Complex/undetermined pattern", "Complex/undetermined pattern", dir_string))
}

# =========================================================
# 13. Merge and classify
# =========================================================
pattern_tbl <- overall_sig_tbl %>%
  left_join(pairwise_tbl, by = "Genus") %>%
  mutate(
    pattern_raw = pmap(
      list(med_Rum, med_Ile, med_Col, sig_RI, sig_RC, sig_IC),
      classify_pattern
    ),
    pattern_class = map_chr(pattern_raw, 1),
    pattern_subclass = map_chr(pattern_raw, 2),
    pattern_direction = map_chr(pattern_raw, 3)
  ) %>%
  select(
    Genus,
    p_friedman, q_friedman, overall_sig_friedman,
    prevalence_all, prevalence_Rum, prevalence_Ile, prevalence_Col,
    med_Rum, med_Ile, med_Col,
    p_RI, q_RI, sig_RI,
    p_RC, q_RC, sig_RC,
    p_IC, q_IC, sig_IC,
    lfc_GroupIle, lfc_GroupCol,
    p_GroupIle, p_GroupCol,
    q_GroupIle, q_GroupCol,
    diff_GroupIle, diff_GroupCol,
    passed_ss_GroupIle, passed_ss_GroupCol,
    diff_robust_GroupIle, diff_robust_GroupCol,
    pattern_class, pattern_subclass, pattern_direction,
    Phylum, Class, Order, Family
  ) %>%
  arrange(pattern_class, pattern_subclass, q_friedman, p_friedman)

# =========================================================
# 14. Summary tables
# =========================================================
pattern_summary_subclass_tbl <- pattern_tbl %>%
  count(pattern_class, pattern_subclass, name = "n_genus") %>%
  arrange(desc(n_genus), pattern_class, pattern_subclass)

pattern_summary_class_tbl <- pattern_tbl %>%
  count(pattern_class, name = "n_genus") %>%
  mutate(percent = round(100 * n_genus / sum(n_genus), 2)) %>%
  arrange(desc(n_genus))

pattern_top_tbl <- pattern_tbl %>%
  group_by(pattern_class, pattern_subclass) %>%
  arrange(q_friedman, .by_group = TRUE) %>%
  slice_head(n = top_n_each_subclass) %>%
  ungroup()

# Basic table for all genera
all_genus_desc_tbl <- pattern_base_tbl %>%
  left_join(prev_df, by = "Genus") %>%
  left_join(tax_df, by = "Genus") %>%
  left_join(ancom_main_tbl %>% select(
    Genus, lfc_GroupIle, lfc_GroupCol, p_GroupIle, p_GroupCol,
    q_GroupIle, q_GroupCol, diff_GroupIle, diff_GroupCol
  ), by = "Genus") %>%
  left_join(friedman_tbl, by = "Genus") %>%
  arrange(Genus)

# =========================================================
# 15. Save results
# =========================================================
write.xlsx(
  qc_sample_alignment,
  file = file.path(out_dir, "01_QC_sample_alignment.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  qc_preprocess,
  file = file.path(out_dir, "02_QC_preprocess_summary.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  ancom_raw_tbl,
  file = file.path(out_dir, "03_ANCOMBC2_raw_results_all_genus.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  ancom_main_tbl,
  file = file.path(out_dir, "04_ANCOMBC2_main_results_cleaned.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  friedman_tbl,
  file = file.path(out_dir, "05_Friedman_overall_all_genus.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  overall_sig_tbl,
  file = file.path(out_dir, "06_Friedman_overall_significant_genus.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  pairwise_tbl,
  file = file.path(out_dir, "07_pairwise_paired_wilcoxon_for_overall_sig_genus.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  pattern_tbl,
  file = file.path(out_dir, "08_spatial_pattern_classification_main_table.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  pattern_summary_subclass_tbl,
  file = file.path(out_dir, "09_spatial_pattern_summary_subclass.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  pattern_summary_class_tbl,
  file = file.path(out_dir, "10_spatial_pattern_summary_class.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  pattern_top_tbl,
  file = file.path(out_dir, "11_spatial_pattern_top_genus_by_subclass.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  all_genus_desc_tbl,
  file = file.path(out_dir, "12_all_genus_basic_description.xlsx"),
  rowNames = FALSE
)

# Combined workbook
wb <- createWorkbook()

addWorksheet(wb, "QC_sample_alignment")
writeData(wb, "QC_sample_alignment", qc_sample_alignment)

addWorksheet(wb, "QC_preprocess_summary")
writeData(wb, "QC_preprocess_summary", qc_preprocess)

addWorksheet(wb, "ANCOMBC2_raw_all")
writeData(wb, "ANCOMBC2_raw_all", ancom_raw_tbl)

addWorksheet(wb, "ANCOMBC2_main_cleaned")
writeData(wb, "ANCOMBC2_main_cleaned", ancom_main_tbl)

addWorksheet(wb, "Friedman_all")
writeData(wb, "Friedman_all", friedman_tbl)

addWorksheet(wb, "Friedman_sig")
writeData(wb, "Friedman_sig", overall_sig_tbl)

addWorksheet(wb, "pairwise_wilcoxon")
writeData(wb, "pairwise_wilcoxon", pairwise_tbl)

addWorksheet(wb, "pattern_main_table")
writeData(wb, "pattern_main_table", pattern_tbl)

addWorksheet(wb, "pattern_summary_subclass")
writeData(wb, "pattern_summary_subclass", pattern_summary_subclass_tbl)

addWorksheet(wb, "pattern_summary_class")
writeData(wb, "pattern_summary_class", pattern_summary_class_tbl)

addWorksheet(wb, "pattern_top_by_subclass")
writeData(wb, "pattern_top_by_subclass", pattern_top_tbl)

addWorksheet(wb, "all_genus_description")
writeData(wb, "all_genus_description", all_genus_desc_tbl)

saveWorkbook(
  wb,
  file = file.path(out_dir, "00_spatial_pattern_analysis_all_results.xlsx"),
  overwrite = TRUE
)

# =========================================================
# 16. Console summary output
# =========================================================
cat("\n================ Analysis completed ================\n")
cat("Output directory: ", out_dir, "\n")
cat("Raw genus count: ", ntaxa(ps), "\n")
cat("Number of Friedman overall significant genera (q <= ", overall_q_cut, "): ", nrow(overall_sig_tbl), "\n", sep = "")

cat("\nCounts for major pattern categories:\n")
print(pattern_summary_class_tbl)

cat("\nMain output files:\n")
cat("03_ANCOMBC2_raw_results_all_genus.xlsx\n")
cat("04_ANCOMBC2_main_results_cleaned.xlsx\n")
cat("05_Friedman_overall_all_genus.xlsx\n")
cat("06_Friedman_overall_significant_genus.xlsx\n")
cat("07_pairwise_paired_wilcoxon_for_overall_sig_genus.xlsx\n")
cat("08_spatial_pattern_classification_main_table.xlsx\n")
cat("09_spatial_pattern_summary_subclass.xlsx\n")
cat("10_spatial_pattern_summary_class.xlsx\n")
cat("11_spatial_pattern_top_genus_by_subclass.xlsx\n")
cat("12_all_genus_basic_description.xlsx\n")
cat("00_spatial_pattern_analysis_all_results.xlsx\n")
cat("=========================================\n")