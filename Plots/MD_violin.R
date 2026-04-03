library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(officer)
library(rvg)

# ==============================
# Function 1: CC aggregated violin plot
# ==============================
plot_violin_cc <- function(folder_path) {
  # 1. Read all CSV files
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  df_list <- lapply(csv_files, function(file) {
    df <- read_csv(file, show_col_types = FALSE)
    df$group <- tools::file_path_sans_ext(basename(file))
    df
  })
  
  # 2. Merge and normalize (min-max per group), keep abs values < 1
  df_all <- bind_rows(df_list) %>%
    mutate(value = abs(value)) %>%
    filter(value < 1) %>%
    group_by(group) %>%
    mutate(
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      value_minmax = (value - min_val) / (max_val - min_val)
    ) %>%
    ungroup()
  
  # 3. Trim middle 50% (25%-75%)
  df_trimmed <- df_all %>%
    group_by(group) %>%
    mutate(
      q1 = quantile(value_minmax, 0.25, na.rm = TRUE),
      q3 = quantile(value_minmax, 0.75, na.rm = TRUE)
    ) %>%
    filter(value_minmax >= q1 & value_minmax <= q3) %>%
    ungroup()
  
  # 4. Violin + boxplot
  p <- ggplot(df_trimmed, aes(x = group, y = value_minmax, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
    labs(title = "Violin Plot (Middle 50%, Min-Max Normalized)",
         x = "Group", y = "Value (min-max normalized, trimmed)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 5. Pairwise Wilcoxon tests
  df_trimmed$group <- factor(df_trimmed$group)
  comparisons <- combn(levels(df_trimmed$group), 2, simplify = FALSE)
  p <- p + stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    size = 3,
    format.args = list(digits = 6)
  )
  
  # 6. Export p-values table
  pval_df <- do.call(rbind, lapply(comparisons, function(comp) {
    g1 <- comp[1]; g2 <- comp[2]
    v1 <- df_trimmed$value_minmax[df_trimmed$group == g1]
    v2 <- df_trimmed$value_minmax[df_trimmed$group == g2]
    if (length(v1) > 1 && length(v2) > 1) {
      pv <- wilcox.test(v1, v2)$p.value
      pf <- formatC(pv, format = "e", digits = 6)
    } else {
      pf <- NA
    }
    data.frame(group1 = g1, group2 = g2, p_value = pf)
  }))
  write.csv(pval_df, file.path(folder_path, "pairwise_pvalues_trimmed_minmax.csv"), row.names = FALSE)
  
  # 7. Export editable PPT
  ppt_path <- file.path(folder_path, "violin_trimmed_minmax_editable.pptx")
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, dml(ggobj = p), location = ph_location_type(type = "body"))
  doc <- ph_with(doc, "Violin Plot (Middle 50%, Min-Max Normalized)",
                 location = ph_location_type(type = "title"))
  print(doc, target = ppt_path)
  
  message("Exported editable PPT: ", ppt_path)
  return(list(plot = p, pvals = pval_df))
}


# ==============================
# Function 2: singleAA aggregated violin plot
# ==============================
plot_violin_singleaa <- function(folder_path) {
  # 1. Read all CSV files
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  df_list <- lapply(csv_files, function(file) {
    df <- read_csv(file, show_col_types = FALSE)
    df$group <- tools::file_path_sans_ext(basename(file))
    # Automatically detect numeric column
    cand <- intersect(names(df), c("value","ANM_effectiveness","ANM_sensitivity","ANM_stiffness","ANM_sq"))
    valcol <- if (length(cand) >= 1) cand[1] else {
      num_cols <- names(df)[sapply(df, is.numeric)]
      if (length(num_cols) == 0) stop(paste("No numeric column found in:", basename(file)))
      num_cols[1]
    }
    df$value <- df[[valcol]]
    df[, c("value","group")]
  })
  df_all <- bind_rows(df_list)   # no abs, no normalization
  
  # 2. Trim middle 90% (remove top and bottom 5%)
  df_trimmed <- df_all %>%
    group_by(group) %>%
    mutate(
      q1 = quantile(value, 0.05, na.rm = TRUE),
      q3 = quantile(value, 0.95, na.rm = TRUE)
    ) %>%
    filter(value >= q1 & value <= q3) %>%
    ungroup()
  
  # 3. Violin + boxplot
  p <- ggplot(df_trimmed, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
    labs(title = "Violin Plot (Middle 90% of Data, No Normalization)",
         x = "Group", y = "Value (trimmed)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 4. Pairwise Wilcoxon tests (robust wrapper)
  safe_wilcox <- function(x, y) {
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    out <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
    out
  }
  
  df_trimmed$group <- factor(df_trimmed$group)
  comparisons <- combn(levels(df_trimmed$group), 2, simplify = FALSE)
  p <- p + stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    size = 3,
    format.args = list(digits = 6)
  )
  
  # 5. Export p-values table
  pval_df <- do.call(rbind, lapply(comparisons, function(comp) {
    g1 <- comp[1]; g2 <- comp[2]
    v1 <- df_trimmed$value[df_trimmed$group == g1]
    v2 <- df_trimmed$value[df_trimmed$group == g2]
    pv <- safe_wilcox(v1, v2)
    data.frame(group1 = g1, group2 = g2,
               p_value = ifelse(is.na(pv), NA, formatC(pv, format = "e", digits = 6)))
  }))
  write.csv(pval_df, file.path(folder_path, "pairwise_pvalues_trimmed_no_norm.csv"), row.names = FALSE)
  
  # 6. Export editable PPT
  ppt_path <- file.path(folder_path, "violin_trimmed_no_norm_editable.pptx")
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, dml(ggobj = p), location = ph_location_type(type = "body"))
  doc <- ph_with(doc, "Violin Plot (Middle 90%, No Normalization)",
                 location = ph_location_type(type = "title"))
  print(doc, target = ppt_path)
  
  message("Exported editable PPT: ", ppt_path)
  return(list(plot = p, pvals = pval_df))
}


# ==============================
# Example usage
# ==============================
# Example 1: CC aggregated
# result_cc <- plot_violin_cc("F:/Data_0822/dynamics/data/results/CC_aggregated/anm_top3")

# Example 2: singleAA aggregated
# result_saa <- plot_violin_singleaa("F:/Data_0822/dynamics/data/results/singleAA_aggregated/anm_effectiveness")
