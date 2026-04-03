# ============================ Density Plots Helper ============================
# This function reads a cluster summary table and produces two PPT files:
# (1) Overall density of Cc (log10 x-axis)
# (2) Grouped density of Cc by whether a cluster has any functional sites
#     (Fun_PTM or Fun_Mut non-empty) -> "Functional" vs "Non-functional"
# ============================================================================
plot_cluster_density <- function(input_file,
                                 out_density_all,
                                 out_density_group,
                                 vline_x = NULL) {
  # ---- dependencies ----
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(ggplot2)
    library(officer)
    library(rvg)
  })
  
  # ---- load data ----
  df <- readxl::read_excel(input_file)
  # Ensure Cc is numeric (it may come in as character from Excel)
  df <- df %>% mutate(Cc = suppressWarnings(as.numeric(Cc)))
  
  # ---- overall density ----
  p_all <- ggplot(df, aes(x = Cc)) +
    geom_density(alpha = 0.5, color = "black", fill = "#4D7EAE") +
    theme_bw() +
    scale_x_log10() +
    theme(
      axis.title   = element_text(size = 15, color = "black"),
      axis.text.y  = element_text(size = 15, color = "black"),
      axis.text.x  = element_text(size = 15, color = "black"),
      legend.text  = element_text(size = 15),
      legend.title = element_text(size = 15),
      strip.text   = element_text(size = 12)
    ) +
    ggtitle("Overall Density of Cc (log10 scale)")
  
  if (!is.null(vline_x)) {
    p_all <- p_all + geom_vline(xintercept = vline_x, linetype = "solid", color = "red", size = 0.5)
  }
  
  doc1 <- officer::read_pptx()
  doc1 <- officer::add_slide(doc1, layout = "Title and Content", master = "Office Theme")
  doc1 <- officer::ph_with(doc1, rvg::dml(ggobj = p_all), location = officer::ph_location_fullsize())
  print(doc1, target = out_density_all)
  
  # ---- grouped density by Functional vs Non-functional ----
  df2 <- df %>%
    mutate(
      Group = ifelse(
        (!is.na(Fun_PTM) & Fun_PTM != "") | (!is.na(Fun_Mut) & Fun_Mut != ""),
        "Functional", "Non-functional"
      )
    )
  
  p_grp <- ggplot(df2, aes(x = Cc, fill = Group)) +
    geom_density(alpha = 0.5, color = "black") +
    theme_bw() +
    scale_x_log10() +
    scale_fill_manual(
      values = c("Functional" = "#4D7EAE", "Non-functional" = "#5C996C"),
      name   = "Group"
    ) +
    theme(
      axis.title   = element_text(size = 15, color = "black"),
      axis.text.y  = element_text(size = 10, color = "black"),
      axis.text.x  = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 8),
      strip.text   = element_text(size = 10)
    ) +
    ggtitle("Density of Cc by Functional Group (log10 scale)")
  
  if (!is.null(vline_x)) {
    p_grp <- p_grp + geom_vline(xintercept = vline_x, linetype = "solid", color = "red", size = 0.5)
  }
  
  doc2 <- officer::read_pptx()
  doc2 <- officer::add_slide(doc2, layout = "Title and Content", master = "Office Theme")
  doc2 <- officer::ph_with(doc2, rvg::dml(ggobj = p_grp), location = officer::ph_location_fullsize())
  print(doc2, target = out_density_group)
  
  message("Done. Saved:\n  - ", out_density_all, "\n  - ", out_density_group)
}

# ============================ Example Run =====================================
# # Adjust paths if needed. These match your earlier directory layout.
# setwd("F:/软著/cluster软著示例/数据可视化")
# plot_cluster_density(
#   input_file        = "cluster_summary_cc30.xlsx",
#   out_density_all   = "density_all.pptx",
#   out_density_group = "density_group.pptx",
#   vline_x           = NULL   # or e.g., 1.21374 to draw a reference line
# )
