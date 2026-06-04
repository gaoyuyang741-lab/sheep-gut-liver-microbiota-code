############################################################
# Alpha diversity analysis for Rumen / Ileum / Colon
# - Merge rumen+ileum file and colon file
# - Keep complete triplets only
# - Friedman test for overall repeated-measures comparison
# - Paired Wilcoxon post hoc tests
# - BH/FDR correction
# - Pretty triplet boxplots with paired lines
#
############################################################

# 0) Packages ------------------------------------------------
pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "ggplot2", "writexl"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# 1) Paths ---------------------------------------------------
# Input and output directories
# Please place the raw alpha-diversity tables in "data/diversity/alpha/input".
# Output files will be saved in "results/diversity/alpha".
# Figure files will be saved in "figures/diversity/alpha/plots_pretty_fdr".

in_file_ri <- file.path(
  "data", "diversity", "alpha", "input",
  "RI_alpha_diversity.xls"
)

in_file_c <- file.path(
  "data", "diversity", "alpha", "input",
  "C_alpha_diversity.xls"
)

out_dir  <- file.path("results", "diversity", "alpha")
plot_dir <- file.path("figures", "diversity", "alpha", "plots_pretty_fdr")

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# 2) Helper: robust reader ----------------------------------
read_alpha_file <- function(in_file) {
  # ---------- attempt 1: read as real Excel ----------
  df <- tryCatch(
    {
      message("Attempting to read as Excel: ", in_file)
      readxl::read_excel(in_file, sheet = 1)
    },
    error = function(e) NULL
  )
  
  # ---------- attempt 2: read as tab-delimited text ----------
  if (is.null(df)) {
    df <- tryCatch(
      {
        message("Excel reading failed; attempting to read as tab-delimited text: ", in_file)
        read.delim(
          in_file,
          header = TRUE,
          sep = "\t",
          check.names = FALSE,
          quote = "",
          comment.char = "",
          stringsAsFactors = FALSE
        )
      },
      error = function(e) NULL
    )
  }
  
  # ---------- attempt 3: fread auto-detect ----------
  if (is.null(df)) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      install.packages("data.table", dependencies = TRUE)
    }
    df <- tryCatch(
      {
        message("Tab-delimited reading failed; attempting automatic detection with data.table::fread: ", in_file)
        data.table::fread(
          in_file,
          data.table = FALSE,
          check.names = FALSE
        )
      },
      error = function(e) NULL
    )
  }
  
  # ---------- fail ----------
  if (is.null(df)) {
    stop(paste0("Unable to read file: ", in_file,
                "\nPlease check whether the file can be opened normally in Excel, or save it as .xlsx and rerun."))
  }
  
  # ---------- normalize column names ----------
  names(df) <- stringr::str_replace_all(names(df), "\\s+", "_")
  names(df) <- stringr::str_replace_all(names(df), "-", "_")
  
  # ensure Sample_ID exists
  if (!"Sample_ID" %in% names(df)) {
    if ("Sample" %in% names(df)) {
      df <- df %>% rename(Sample_ID = Sample)
    } else if ("SampleID" %in% names(df)) {
      df <- df %>% rename(Sample_ID = SampleID)
    } else {
      names(df)[1] <- "Sample_ID"
    }
  }
  
  # identify PD column variant
  pd_candidates <- c("PD_whole_tree", "PD_whole_", "PD_whole_tree_index", "PD_whole")
  pd_col <- intersect(pd_candidates, names(df))
  if (length(pd_col) == 0) {
    pd_col2 <- grep("^PD", names(df), value = TRUE)
    if (length(pd_col2) == 0) {
      stop(paste0("PD column not found. Current column names: ", paste(names(df), collapse = ", ")))
    } else {
      pd_col <- pd_col2[1]
    }
  } else {
    pd_col <- pd_col[1]
  }
  
  # required columns
  need_cols <- c("Sample_ID", "Feature", "ACE", "Chao1", "Simpson", "Shannon", "Coverage")
  missing_cols <- setdiff(need_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing columns: ", paste(missing_cols, collapse = ", "),
                "\nCurrent column names: ", paste(names(df), collapse = ", ")))
  }
  
  # standardize PD name
  df <- df %>% rename(PD_whole_tree = all_of(pd_col))
  
  # numeric columns
  num_cols <- c("Feature", "ACE", "Chao1", "Simpson", "Shannon", "PD_whole_tree", "Coverage")
  df <- df %>%
    mutate(across(any_of(num_cols), ~ suppressWarnings(as.numeric(.))))
  
  return(df)
}

# 3) Read and merge data ------------------------------------
df_ri <- read_alpha_file(in_file_ri)
df_c  <- read_alpha_file(in_file_c)

df_all <- bind_rows(df_ri, df_c)

# 4) Parse site and Sheep_ID --------------------------------
# Expected sample naming:
#   R-XXXX   -> Rum
#   I-XXXX   -> Ile
#   C-XXXX   -> Colon
df_all <- df_all %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    site_code = stringr::str_extract(Sample_ID, "^[RIC]"),
    Site = case_when(
      site_code == "R" ~ "Rum",
      site_code == "I" ~ "Ile",
      site_code == "C" ~ "Colon",
      TRUE ~ NA_character_
    ),
    Sheep_ID = stringr::str_replace(Sample_ID, "^[RIC]-", "")
  )

# basic checks
cat("Total rows after merge:", nrow(df_all), "\n")
cat("Site counts:\n")
print(table(df_all$Site, useNA = "ifany"))

# 5) Keep complete triplets only -----------------------------
triplet_table <- df_all %>%
  filter(!is.na(Site)) %>%
  count(Sheep_ID, Site) %>%
  tidyr::pivot_wider(names_from = Site, values_from = n, values_fill = 0)

needed_sites <- c("Rum", "Ile", "Colon")
for (s in needed_sites) {
  if (!s %in% names(triplet_table)) {
    triplet_table[[s]] <- 0
  }
}

complete_ids <- triplet_table %>%
  filter(Rum == 1, Ile == 1, Colon == 1) %>%
  pull(Sheep_ID)

df_triplet <- df_all %>%
  filter(Sheep_ID %in% complete_ids, Site %in% c("Rum", "Ile", "Colon")) %>%
  mutate(Site = factor(Site, levels = c("Rum", "Ile", "Colon"))) %>%
  arrange(Sheep_ID, Site)

write.csv(triplet_table, file.path(out_dir, "triplet_check_table.csv"), row.names = FALSE)

cat("Samples kept in complete triplets:", nrow(df_triplet), "\n")
cat("Number of complete sheep triplets :", length(unique(df_triplet$Sheep_ID)), "\n")

# 6) Indices -------------------------------------------------
indices <- c("Feature", "ACE", "Chao1", "Simpson", "Shannon", "PD_whole_tree")

# 7) P-value formatter ---------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) {
    formatC(p, format = "e", digits = 2)
  } else {
    format(signif(p, 3), scientific = FALSE, trim = TRUE)
  }
}

# 8) Friedman overall test ----------------------------------
friedman_one <- function(dat, idx) {
  subdat <- dat %>%
    select(Sheep_ID, Site, all_of(idx)) %>%
    tidyr::pivot_wider(names_from = Site, values_from = all_of(idx))
  
  subdat <- subdat %>%
    tidyr::drop_na(Rum, Ile, Colon)
  
  longdat <- subdat %>%
    tidyr::pivot_longer(cols = c(Rum, Ile, Colon), names_to = "Site", values_to = "value") %>%
    mutate(
      Site = factor(Site, levels = c("Rum", "Ile", "Colon")),
      Sheep_ID = factor(Sheep_ID)
    )
  
  p_friedman <- tryCatch(
    friedman.test(value ~ Site | Sheep_ID, data = longdat)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    index = idx,
    n_sheep = nrow(subdat),
    mean_Rum   = mean(subdat$Rum,   na.rm = TRUE),
    mean_Ile   = mean(subdat$Ile,   na.rm = TRUE),
    mean_Colon = mean(subdat$Colon, na.rm = TRUE),
    median_Rum   = median(subdat$Rum,   na.rm = TRUE),
    median_Ile   = median(subdat$Ile,   na.rm = TRUE),
    median_Colon = median(subdat$Colon, na.rm = TRUE),
    p_friedman = p_friedman
  )
}

res_overall <- purrr::map_dfr(indices, ~ friedman_one(df_triplet, .x)) %>%
  mutate(
    p_friedman_fdr = p.adjust(p_friedman, method = "BH")
  ) %>%
  arrange(p_friedman_fdr)

print(res_overall)

write.csv(res_overall, file.path(out_dir, "alpha_friedman_overall.csv"), row.names = FALSE)
writexl::write_xlsx(res_overall, file.path(out_dir, "alpha_friedman_overall.xlsx"))

# 9) Pairwise post hoc tests --------------------------------
pairwise_wilcox_one <- function(dat, idx, site1, site2) {
  subdat <- dat %>%
    filter(Site %in% c(site1, site2)) %>%
    select(Sheep_ID, Site, all_of(idx)) %>%
    tidyr::pivot_wider(names_from = Site, values_from = all_of(idx))
  
  if (!all(c(site1, site2) %in% names(subdat))) {
    return(tibble(
      index = idx,
      contrast = paste(site1, "vs", site2),
      n_pairs = 0,
      mean_1 = NA_real_,
      mean_2 = NA_real_,
      median_1 = NA_real_,
      median_2 = NA_real_,
      p_raw = NA_real_
    ))
  }
  
  subdat <- subdat %>% tidyr::drop_na(all_of(site1), all_of(site2))
  
  p_raw <- tryCatch(
    wilcox.test(subdat[[site1]], subdat[[site2]], paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    index = idx,
    contrast = paste(site1, "vs", site2),
    n_pairs = nrow(subdat),
    mean_1 = mean(subdat[[site1]], na.rm = TRUE),
    mean_2 = mean(subdat[[site2]], na.rm = TRUE),
    median_1 = median(subdat[[site1]], na.rm = TRUE),
    median_2 = median(subdat[[site2]], na.rm = TRUE),
    p_raw = p_raw
  )
}

contrasts <- list(
  c("Rum", "Ile"),
  c("Rum", "Colon"),
  c("Ile", "Colon")
)

res_pairwise <- purrr::map_dfr(indices, function(idx) {
  purrr::map_dfr(contrasts, function(ct) {
    pairwise_wilcox_one(df_triplet, idx, ct[1], ct[2])
  })
}) %>%
  group_by(index) %>%
  mutate(
    p_fdr_within_index = p.adjust(p_raw, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_fdr_global = p.adjust(p_raw, method = "BH")
  )

print(res_pairwise)

write.csv(res_pairwise, file.path(out_dir, "alpha_pairwise_wilcox.csv"), row.names = FALSE)
writexl::write_xlsx(res_pairwise, file.path(out_dir, "alpha_pairwise_wilcox.xlsx"))

# 10) Plotting -----------------------------------------------
pal <- c(
  "Rum"   = "#5B8FF9",  # cool blue
  "Ile"   = "#F6BD16",  # warm yellow-orange
  "Colon" = "#9270CA"   # muted purple
)

get_overall_label <- function(idx, res_overall) {
  p <- res_overall %>%
    filter(index == idx) %>%
    pull(p_friedman_fdr)
  
  paste0("Friedman p(FDR) = ", fmt_p(p))
}

plot_triplet_box <- function(dat, idx, res_overall) {
  p_lab <- get_overall_label(idx, res_overall)
  
  ymax <- max(dat[[idx]], na.rm = TRUE)
  ymin <- min(dat[[idx]], na.rm = TRUE)
  yrange <- ymax - ymin
  if (yrange == 0) yrange <- abs(ymax) * 0.1 + 1e-6
  
  y_text <- ymax + 0.12 * yrange
  
  ggplot(dat, aes(x = Site, y = .data[[idx]], fill = Site)) +
    geom_line(aes(group = Sheep_ID), color = "grey80", linewidth = 0.4, alpha = 0.8) +
    geom_point(
      aes(group = Sheep_ID),
      color = "black", size = 1.7, alpha = 0.75,
      position = position_jitter(width = 0.05, height = 0)
    ) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.88, linewidth = 0.7) +
    scale_fill_manual(values = pal) +
    annotate("text", x = 2, y = y_text, label = p_lab, size = 4) +
    labs(
      title = idx,
      x = NULL,
      y = idx
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    coord_cartesian(ylim = c(ymin, y_text + 0.08 * yrange))
}

for (idx in indices) {
  p <- plot_triplet_box(df_triplet, idx, res_overall)
  
  ggsave(
    file.path(plot_dir, paste0("triplet_alpha_", idx, ".pdf")),
    p, width = 7.0, height = 5.4
  )
  ggsave(
    file.path(plot_dir, paste0("triplet_alpha_", idx, ".png")),
    p, width = 7.0, height = 5.4, dpi = 300
  )
}

# 11) Optional summary file ----------------------------------
summary_txt <- c(
  paste0("Total merged rows: ", nrow(df_all)),
  paste0("Complete triplet sheep: ", length(unique(df_triplet$Sheep_ID))),
  paste0("Complete triplet rows: ", nrow(df_triplet))
)
writeLines(summary_txt, file.path(out_dir, "run_summary.txt"))

message("All done.")
message("Output folder: ", normalizePath(out_dir))
message("Plots folder : ", normalizePath(plot_dir))