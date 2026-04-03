# ====================== Cluster Charts: Pie + Faceted Bars (+Venn) ======================
# This function:
#   1) Reads a cluster summary Excel (must contain: Pos_Type, Type, Fun_PTM, Fun_Mut,
#      Pocket_Site, Fpocket_Site, etc.)
#   2) Builds two pie charts:
#        - Distribution by Pos_Type
#        - Distribution by Type
#      and saves each as a PPT (editable vector graphics via rvg::dml).
#   3) Builds a faceted stacked bar chart (Pos_Type / Type / IF_Func vs IF_Pockets),
#      also saved as a PPT.
#   4) Optionally saves a Venn diagram (Pocket_Site vs Fpocket_Site) to PNG.
# All plots are self-contained and can be edited in PowerPoint.
# ========================================================================================

build_cluster_charts <- function(input_xlsx,
                                 out_pie_pos_ppt,
                                 out_pie_type_ppt,
                                 out_faceted_ppt,
                                 out_venn_png = NULL) {
  # ---- Dependencies ----
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(officer)
    library(rvg)
    library(ggh4x)
    library(scales)
  })
  
  # ---- Read data ----
  df <- readxl::read_excel(input_xlsx)
  
  # ---- Helper: safe string emptiness ----
  is_filled <- function(x) !is.na(x) & x != ""
  
  # ========================= 1) PIE: Pos_Type distribution =========================
  # Custom colors for Pos_Type
  pos_colors <- c(
    "Cross-Interface" = "#C54A68",
    "Interface"       = "#7BCB7A",
    "Non-Interface"   = "#4D7EAE",
    "NA"              = "#BDBDBD"  # fallback for missing
  )
  
  # Prepare data (keep NA as a visible category)
  df_pos <- df %>%
    mutate(Pos_Type = ifelse(is_filled(Pos_Type), Pos_Type, "NA")) %>%
    count(Pos_Type, .drop = FALSE) %>%
    mutate(pct = n / sum(n) * 100)
  
  p_pos <- ggplot(df_pos, aes(x = "", y = n, fill = Pos_Type)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = pos_colors) +
    geom_text(
      aes(label = paste0(Pos_Type, "\n", sprintf("%.1f%%", pct))),
      position = position_stack(vjust = 0.5), size = 4
    ) +
    labs(
      title    = "Distribution of Clusters by Pos_Type",
      subtitle = "Percentage of each Pos_Type class",
      fill     = "Pos_Type"
    ) +
    theme_void() +
    theme(legend.position = "right")
  
  # Save to PPT
  doc_pos <- officer::read_pptx()
  doc_pos <- officer::add_slide(doc_pos, layout = "Title and Content", master = "Office Theme")
  doc_pos <- officer::ph_with(doc_pos, rvg::dml(ggobj = p_pos), location = officer::ph_location_fullsize())
  print(doc_pos, target = out_pie_pos_ppt)
  
  # =========================== 2) PIE: Type distribution ===========================
  type_colors <- c(
    "Mut"    = "#C54A68",
    "PTM"    = "#7BCB7A",
    "Hybrid" = "#4D7EAE",
    "NA"     = "#BDBDBD"
  )
  
  df_type <- df %>%
    mutate(Type = ifelse(is_filled(Type), Type, "NA")) %>%
    count(Type, .drop = FALSE) %>%
    mutate(pct = n / sum(n) * 100)
  
  p_type <- ggplot(df_type, aes(x = "", y = n, fill = Type)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = type_colors) +
    geom_text(
      aes(label = paste0(Type, "\n", sprintf("%.1f%%", pct))),
      position = position_stack(vjust = 0.5), size = 4
    ) +
    labs(
      title    = "Distribution of Clusters by Type",
      subtitle = "Percentage of each composition type",
      fill     = "Type"
    ) +
    theme_void() +
    theme(legend.position = "right")
  
  # Save to PPT
  doc_type <- officer::read_pptx()
  doc_type <- officer::add_slide(doc_type, layout = "Title and Content", master = "Office Theme")
  doc_type <- officer::ph_with(doc_type, rvg::dml(ggobj = p_type), location = officer::ph_location_fullsize())
  print(doc_type, target = out_pie_type_ppt)
  
  # ===================== 3) Faceted stacked bars: IF_Pockets ======================
  # IF_Pockets: 1 if either pocket column is non-empty; else 0
  df_aug <- df %>%
    mutate(
      IF_Pockets = as.integer(
        (is_filled(Pocket_Site)) | (is_filled(Fpocket_Site))
      ),
      IF_Func = ifelse(
        (is_filled(Fun_PTM)) | (is_filled(Fun_Mut)),
        "Func_Cluster", "Non-func_Cluster"
      ),
      Pos_Type = ifelse(is_filled(Pos_Type), Pos_Type, "NA"),
      Type     = ifelse(is_filled(Type),     Type,     "NA")
    )
  
  # Long format for three categories
  df_long <- df_aug %>%
    mutate(IF_Pockets = factor(IF_Pockets, levels = c(0, 1), labels = c("0", "1"))) %>%
    pivot_longer(
      cols      = c(Pos_Type, Type, IF_Func),
      names_to  = "Category",
      values_to = "Value"
    ) %>%
    count(Category, Value, IF_Pockets) %>%
    mutate(
      label_count = paste0("n=", n),
      fill_color = dplyr::case_when(
        IF_Pockets == "0"                          ~ "lightgray",
        IF_Pockets == "1" & Category == "Pos_Type" ~ "#FF7F0E",
        IF_Pockets == "1" & Category == "Type"     ~ "#1F77B4",
        IF_Pockets == "1" & Category == "IF_Func"  ~ "#2CA02C",
        TRUE ~ "lightgray"
      ),
      Category = factor(Category, levels = c("Pos_Type", "Type", "IF_Func")),
      # make sure colored layer is below and gray on top (reverse stacking)
      IF_Pockets = factor(IF_Pockets, levels = c("1", "0"))
    )
  
  # Panel backgrounds per facet
  panel_bg <- tibble::tibble(
    Category = factor(c("Pos_Type", "Type", "IF_Func"),
                      levels = c("Pos_Type", "Type", "IF_Func")),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = c("#FF7F0E", "#1F77B4", "#2CA02C")
  )
  
  p_faceted <- ggplot() +
    # light panel tint per facet
    geom_rect(
      data = panel_bg,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      alpha = 0.08, inherit.aes = FALSE
    ) +
    # stacked bars; reverse so colored ("1") at bottom, gray ("0") on top
    geom_col(
      data     = df_long,
      aes(x = Value, y = n, fill = fill_color, order = IF_Pockets),
      width    = 0.6,
      position = position_stack(reverse = TRUE)
    ) +
    # labels on stacks
    geom_text(
      data     = df_long,
      aes(x = Value, y = n, label = label_count, order = IF_Pockets),
      position = position_stack(vjust = 0.7, reverse = TRUE),
      size     = 3,
      color    = "white"
    ) +
    facet_grid2(
      . ~ Category,
      scales = "free_x",
      space  = "free_x",
      switch = "x",
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = c(Pos_Type = "#FF7F0E",
                     Type     = "#1F77B4",
                     IF_Func  = "#2CA02C"),
          colour = c(Pos_Type = "#FF7F0E",
                     Type     = "#1F77B4",
                     IF_Func  = "#2CA02C"),
          size   = 0
        ),
        text_x = elem_list_text(
          colour = "white",
          size   = 12,
          face   = "bold"
        )
      )
    ) +
    scale_fill_identity() +
    guides(fill = "none") +
    labs(x = NULL, y = "Count") +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.text.x    = element_text(angle = 0),
      panel.spacing.x = unit(1, "cm"),
      axis.text.x     = element_text(size = 10)
    )
  
  doc_bar <- officer::read_pptx()
  doc_bar <- officer::add_slide(doc_bar, layout = "Title and Content", master = "Office Theme")
  doc_bar <- officer::ph_with(doc_bar, rvg::dml(ggobj = p_faceted), location = officer::ph_location_fullsize())
  print(doc_bar, target = out_faceted_ppt)
  
  # =========================== 4) Optional: Venn PNG ===========================
  if (!is.null(out_venn_png)) {
    # VennDiagram is grid-based; save as raster PNG then (optionally) place into PPT elsewhere.
    if (!requireNamespace("VennDiagram", quietly = TRUE) ||
        !requireNamespace("grid", quietly = TRUE)) {
      warning("VennDiagram/grid not available; skipping Venn export.")
    } else {
      # indices where columns are filled
      pocket_idx  <- which(is_filled(df$Pocket_Site))
      fpocket_idx <- which(is_filled(df$Fpocket_Site))
      png(out_venn_png, width = 1600, height = 1200, res = 200)
      grid::grid.newpage()
      VennDiagram::draw.pairwise.venn(
        area1        = length(pocket_idx),
        area2        = length(fpocket_idx),
        cross.area   = length(intersect(pocket_idx, fpocket_idx)),
        category     = c("Pocket_Site", "Fpocket_Site"),
        fill         = c("skyblue", "lightgreen"),
        alpha        = rep(0.5, 2),
        cat.pos      = c(-20, 20),
        cat.dist     = rep(0.05, 2),
        scaled       = TRUE
      )
      dev.off()
    }
  }
  
  message("Saved outputs:",
          "\n  - Pos_Type pie PPT: ", out_pie_pos_ppt,
          "\n  - Type pie PPT:     ", out_pie_type_ppt,
          "\n  - Faceted bars PPT: ", out_faceted_ppt,
          if (!is.null(out_venn_png)) paste0("\n  - Venn PNG:          ", out_venn_png) else "")
  
  invisible(list(
    pie_pos_data   = df_pos,
    pie_type_data  = df_type,
    faceted_data   = df_long
  ))
}

# =============================== Example run ==================================
# Adjust paths as needed. These match your earlier structure.
# setwd("F:/软著/cluster软著示例/数据可视化")
# build_cluster_charts(
#   input_xlsx      = "cluster_summary_cc30.xlsx",
#   out_pie_pos_ppt = "PosType_PieChart_Edit.pptx",
#   out_pie_type_ppt= "ComType_PieChart_Edit.pptx",
#   out_faceted_ppt = "Pocket_Plot.pptx",
#   out_venn_png    = "Pocket_Venn.png"
# )
