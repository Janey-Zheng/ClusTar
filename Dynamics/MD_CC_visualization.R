#!/usr/bin/env Rscript

# ==== Libraries ====
library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(tools)
library(argparse)

# ==== Argument Parser ====
parser <- ArgumentParser(description = "Generate violin plots with min-max normalization and pairwise Wilcoxon tests.")
parser$add_argument("--folder", required = TRUE, help = "Path to the folder containing CSV files")
parser$add_argument("--min_n", type = "integer", default = 5, help = "Minimum number of samples per group to include (default: 5)")
parser$add_argument("--q_low", type = "double", default = 0.05, help = "Lower quantile for trimming (default: 0.05)")
parser$add_argument("--q_high", type = "double", default = 0.95, help = "Upper quantile for trimming (default: 0.95)")
parser$add_argument("--out_csv", default = "pairwise_pvalues_trimmed_minmax.csv", help = "Output CSV file name for p-values")

args <- parser$parse_args()

folder_path <- args$folder
min_n <- args$min_n
q_low <- args$q_low
q_high <- args$q_high
out_csv <- args$out_csv

csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(csv_files) > 0)

# ==== Read and combine data: take the first column as 'value' ====
df_all <- lapply(csv_files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  tibble::tibble(
    value = abs(as.numeric(df[[1]])),   # force numeric
    group = file_path_sans_ext(basename(file))
  )
}) |> bind_rows() |>
  filter(is.finite(value))

# ==== Min-max normalization (remove zero-variance groups) ====
df_all <- df_all |>
  group_by(group) |>
  mutate(min_val = min(value, na.rm = TRUE),
         max_val = max(value, na.rm = TRUE),
         rng = max_val - min_val) |>
  filter(rng > 0) |>
  mutate(value_minmax = (value - min_val) / rng) |>
  ungroup()

# ==== Trim by quantile range (default: 5%–95%) ====
df_trimmed <- df_all |>
  group_by(group) |>
  mutate(q1 = quantile(value_minmax, q_low, na.rm = TRUE),
         q3 = quantile(value_minmax, q_high, na.rm = TRUE)) |>
  filter(value_minmax >= q1, value_minmax <= q3) |>
  ungroup()

# ==== Keep groups with enough samples ====
df_trimmed <- df_trimmed |>
  add_count(group, name = "n_in_group") |>
  filter(n_in_group >= min_n)

# ==== Violin plot ====
p <- ggplot(df_trimmed, aes(x = group, y = value_minmax, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  labs(title = "Violin Plot (Trimmed & Min-Max Normalized)",
       x = "Group", y = "Value (min-max normalized, trimmed)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# ==== Pairwise Wilcoxon tests ====
df_trimmed$group <- droplevels(factor(df_trimmed$group))
comparisons <- combn(levels(df_trimmed$group), 2, simplify = FALSE)

p <- p + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",
  label = "p.format",
  size = 3,
  format.args = list(digits = 6)
)
print(p)

# ==== Export p-values table ====
pval_df <- do.call(rbind, lapply(comparisons, function(comp) {
  g1 <- comp[1]; g2 <- comp[2]
  v1 <- df_trimmed$value_minmax[df_trimmed$group == g1]
  v2 <- df_trimmed$value_minmax[df_trimmed$group == g2]
  if (length(v1) > 1 && length(v2) > 1) {
    pv <- wilcox.test(v1, v2)$p.value
    pf <- formatC(pv, format = "e", digits = 6)
  } else pf <- NA
  data.frame(group1 = g1, group2 = g2, p_value = pf)
}))
write.csv(pval_df, file.path(folder_path, out_csv), row.names = FALSE)
