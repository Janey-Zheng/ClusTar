# =========================== Heatmap Pipeline (All-in-One) ===========================
# This script:
#   1) Loads a cluster summary table.
#   2) Counts DrugBank-mapped genes per PPI (with homodimer fix).
#   3) Splits PPIs into gene_count groups (0/1/2), and within each, into Type (Mut/PTM/Hybrid).
#   4) Deduplicates per PPI by keeping rows with max Cc, takes top-N by Cc, and writes top files.
#   5) For each top file, parses Site_New into a PPI x Cancer wide table of sites and writes it.
#   6) Builds frequency matrices by joining with pan-cancer rate files (Mut or PTM).
#   7) Draws min-max normalized heatmaps and saves as editable PPT (DrawingML) with clear names.
# Notes:
#   - No setwd() calls in the pipeline; paths are explicit.
#   - All comments are in English.
#   - Dependencies: readxl, dplyr, tidyr, stringr, purrr, writexl, openxlsx,
#                   pheatmap, RColorBrewer, rvg, officer, grid
# =====================================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(writexl)
  library(openxlsx)
  library(pheatmap)
  library(RColorBrewer)
  library(rvg)
  library(officer)
  library(grid)
})

# -------- Helpers: basic I/O and checks --------------------------------------
.read_xlsx <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  readxl::read_excel(path)
}

.write_xlsx <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writexl::write_xlsx(df, path)
}

.write_xlsx_openxlsx <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  openxlsx::write.xlsx(df, path, overwrite = TRUE)
}

# -------- Step 1–3: counts and splits ----------------------------------------
add_drugbank_counts <- function(df) {
  req <- c("Gene1","Gene2","Drugbank_gene_mapped","Cc","Type","Site_New")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing columns in cluster summary: ", paste(miss, collapse = ", "))
  df %>%
    rowwise() %>%
    mutate(
      drugbank_count = case_when(
        is.na(Drugbank_gene_mapped) | Drugbank_gene_mapped == "" ~ 0L,
        TRUE ~ length(strsplit(Drugbank_gene_mapped, ",")[[1]])
      ),
      # If homodimer (Gene1 == Gene2) and has mapping, force gene_count = 2
      gene_count = if_else(Gene1 == Gene2 & drugbank_count > 0, 2L, drugbank_count),
      PPI = paste0(pmin(Gene1, Gene2), "-", pmax(Gene1, Gene2))
    ) %>%
    ungroup()
}

split_by_gene_and_type <- function(df) {
  # Returns nested list: [[ "0"|"1"|"2" ]][["Mut"|"PTM"|"Hybrid"]] -> data.frame
  out <- list()
  for (g in c(0L, 1L, 2L)) {
    gdf <- df %>% filter(gene_count == g)
    out[[as.character(g)]] <- list(
      Mut    = gdf %>% filter(Type == "Mut"),
      PTM    = gdf %>% filter(Type == "PTM"),
      Hybrid = gdf %>% filter(Type == "Hybrid")
    )
  }
  out
}

# -------- Step 4: deduplicate per-PPI by max Cc and take top-N ----------------
dedup_max_cc <- function(df) {
  if (!nrow(df)) return(df)
  df[with(df, ave(Cc, PPI, FUN = max) == Cc), ]
}

top_n_by_cc <- function(df, n = 30) {
  if (!nrow(df)) return(df)
  df$Cc <- suppressWarnings(as.numeric(df$Cc))
  df <- df[order(-df$Cc), , drop = FALSE]
  df[seq_len(min(n, nrow(df))), , drop = FALSE]
}

# -------- Step 5: parse Site_New into PPI x cancer wide of sites ----------
parse_sites_to_wide <- function(in_top_xlsx, out_wide_xlsx) {
  dat <- .read_xlsx(in_top_xlsx)
  need <- c("PPI","Site_New")
  if (!all(need %in% names(dat))) {
    stop("Columns PPI and Site_New are required in: ", in_top_xlsx)
  }
  df_long <- dat %>%
    filter(!is.na(Site_New) & Site_New != "") %>%
    mutate(pairs = str_extract_all(Site_New, "[^,]+?\\([^)]*\\)")) %>%
    tidyr::unnest(pairs) %>%
    mutate(
      Site    = str_remove(pairs, "\\(.*$"),
      Cancers = str_extract(pairs, "(?<=\\()[^)]*")
    ) %>%
    separate_rows(Cancers, sep = ",") %>%
    mutate(Cancer = str_trim(Cancers)) %>%
    select(PPI, Cancer, Site)
  df_wide <- df_long %>%
    distinct() %>%
    pivot_wider(
      id_cols     = PPI,
      names_from  = Cancer,
      values_from = Site,
      values_fn   = list(Site = ~paste(.x, collapse = ", ")),
      values_fill = ""
    )
  .write_xlsx(df_wide, out_wide_xlsx)
  invisible(df_wide)
}

# -------- Step 6: build frequency matrices -----------------------------------
# 6A) Mutation matrices from "Mut_Freq_Pan_cancer.xlsx" (needs columns Mut_Inf, mutation_rate)
build_mut_matrix <- function(in_wide_xlsx, mut_rate_xlsx, out_matrix_xlsx) {
  wide <- .read_xlsx(in_wide_xlsx)
  cancer_sheets <- readxl::excel_sheets(mut_rate_xlsx)
  rates_long <- purrr::map_dfr(cancer_sheets, function(ca) {
    readxl::read_excel(mut_rate_xlsx, sheet = ca) %>%
      select(Mut_Inf, mutation_rate) %>%
      mutate(Cancer = ca)
  })
  mat <- wide %>%
    pivot_longer(-PPI, names_to = "Cancer", values_to = "SitesString") %>%
    filter(SitesString != "") %>%
    separate_rows(SitesString, sep = ",") %>%
    mutate(Mut_Inf = str_trim(SitesString)) %>%
    left_join(rates_long, by = c("Mut_Inf", "Cancer")) %>%
    group_by(PPI, Cancer) %>%
    summarise(Cluster_Mut_Freq = sum(mutation_rate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols     = PPI,
      names_from  = Cancer,
      values_from = Cluster_Mut_Freq,
      values_fill = 0
    )
  .write_xlsx(mat, out_matrix_xlsx)
  invisible(mat)
}

# 6B) PTM matrices from "PTM_Freq_Pan_9_cancer.xlsx" (needs columns Site, Prop_sig)
build_ptm_matrix <- function(in_wide_xlsx, ptm_rate_xlsx, out_matrix_xlsx) {
  wide <- .read_xlsx(in_wide_xlsx)
  cancer_sheets <- readxl::excel_sheets(ptm_rate_xlsx)
  rates_long <- purrr::map_dfr(cancer_sheets, function(ca) {
    readxl::read_excel(ptm_rate_xlsx, sheet = ca) %>%
      select(Site, Prop_sig) %>%
      mutate(Cancer = ca)
  }) %>% distinct(Site, Cancer, .keep_all = TRUE)
  mat <- wide %>%
    pivot_longer(-PPI, names_to = "Cancer", values_to = "SitesString") %>%
    filter(SitesString != "") %>%
    separate_rows(SitesString, sep = ",") %>%
    mutate(PTM_Inf = str_trim(SitesString)) %>%
    left_join(rates_long, by = c("PTM_Inf" = "Site", "Cancer" = "Cancer")) %>%
    group_by(PPI, Cancer) %>%
    summarise(Cluster_PTM_Freq = sum(Prop_sig, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols     = PPI,
      names_from  = Cancer,
      values_from = Cluster_PTM_Freq,
      values_fill = 0
    )
  .write_xlsx(mat, out_matrix_xlsx)
  invisible(mat)
}

# -------- Step 7: heatmap -> PPT ---------------------------------------------
# Standard column order for cancers (first col is PPI)
STD_ORDER <- c("PPI","BRCA","CCRCC","COAD","GBM","HNSCC","LSCC","LUAD","OV","PDAC","UCEC")

add_missing_cols <- function(df, col_order = STD_ORDER) {
  miss <- setdiff(col_order, colnames(df))
  if (length(miss) > 0) {
    miss_df <- as.data.frame(matrix(0, nrow = nrow(df), ncol = length(miss)))
    colnames(miss_df) <- miss
    df <- bind_cols(df, miss_df)
  }
  df %>% select(all_of(col_order))
}

save_heatmap_to_ppt <- function(in_matrix_xlsx, out_pptx, palette = "Blues") {
  df <- .read_xlsx(in_matrix_xlsx)
  df <- add_missing_cols(df, STD_ORDER)
  mat <- as.matrix(df[ , -1, drop = FALSE])
  rownames(mat) <- df$PPI
  
  # Min-max normalize per row; keep zeros if constant rows
  mat_norm <- t(apply(mat, 1, function(x){
    if (all(x == 0) || max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }))
  rownames(mat_norm) <- rownames(mat); colnames(mat_norm) <- colnames(mat)
  
  cols <- colorRampPalette(brewer.pal(9, palette))(50)
  ph <- pheatmap(
    mat_norm,
    silent        = TRUE,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    border_color  = NA,
    color         = cols,
    na_col        = "#2e2e2e",
    cellwidth     = 15,
    cellheight    = 15
  )
  
  dml_plot <- rvg::dml(width = 12, height = 9, code = {
    grid.newpage()
    grid.draw(rectGrob(gp = gpar(fill = "white", col = NA)))
    grid.draw(ph$gtable)
  })
  doc <- read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = dml_plot, location = ph_location_fullsize())
  print(doc, target = out_pptx)
  invisible(out_pptx)
}

# -------- File-naming helpers (clear, consistent) ----------------------------
# These names are used throughout the pipeline.
fname_top      <- function(out_dir_top, g, type, top_n)    file.path(out_dir_top, sprintf("Top%d_geneCount%s_%s.xlsx", top_n, g, type))
fname_sites    <- function(out_dir_top, g, type)            file.path(out_dir_top, sprintf("Sites_geneCount%s_%s.xlsx", g, type))
fname_matrix_m <- function(out_dir_top, g)                  file.path(out_dir_top, sprintf("Matrix_Mutation_geneCount%s.xlsx", g))
fname_matrix_p <- function(out_dir_top, g)                  file.path(out_dir_top, sprintf("Matrix_PTM_geneCount%s.xlsx", g))
fname_matrix_hm<- function(out_dir_top, g)                  file.path(out_dir_top, sprintf("Matrix_Hybrid-Mutation_geneCount%s.xlsx", g))
fname_matrix_hp<- function(out_dir_top, g)                  file.path(out_dir_top, sprintf("Matrix_Hybrid-PTM_geneCount%s.xlsx", g))
fname_ppt_m    <- function(out_dir_ppt, g)                  file.path(out_dir_ppt, sprintf("Heatmap_Mutation_geneCount%s.pptx", g))
fname_ppt_p    <- function(out_dir_ppt, g)                  file.path(out_dir_ppt, sprintf("Heatmap_PTM_geneCount%s.pptx", g))
fname_ppt_hm   <- function(out_dir_ppt, g)                  file.path(out_dir_ppt, sprintf("Heatmap_Hybrid-Mutation_geneCount%s.pptx", g))
fname_ppt_hp   <- function(out_dir_ppt, g)                  file.path(out_dir_ppt, sprintf("Heatmap_Hybrid-PTM_geneCount%s.pptx", g))

# -------- Orchestrator: end-to-end run ---------------------------------------
run_heatmap_pipeline <- function(
    # input summary table (must contain: Gene1,Gene2,Drugbank_gene_mapped,Cc,Type,Site_New)
  cluster_summary_xlsx,
  # where to write intermediate top tables and wide tables
  out_dir_top      = "F:/Data_0822/pictures/Heatmap",
  # pan-cancer reference files
  mut_rate_xlsx    = "F:/Data_0822/pictures/Heatmap/Mut_Freq_Pan_cancer.xlsx",
  ptm_rate_xlsx    = "F:/Data_0822/pictures/Heatmap/PTM_Freq_Pan_9_cancer.xlsx",
  # output PPT directory
  out_dir_ppt      = "F:/Data_0822/pictures/Heatmap",
  # top-N
  top_n            = 30
) {
  dir.create(out_dir_top, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_ppt, showWarnings = FALSE, recursive = TRUE)
  
  # Load and pre-compute counts
  cc <- .read_xlsx(cluster_summary_xlsx)
  cc$Cc <- suppressWarnings(as.numeric(cc$Cc))   # ensure numeric
  cc_ex <- add_drugbank_counts(cc)
  
  # Split by gene_count and Type
  parts <- split_by_gene_and_type(cc_ex)
  
  # Track generated matrices (for the final plotting pass)
  generated <- list(Mut = character(0), PTM = character(0),
                    Hybrid_Mut = character(0), Hybrid_PTM = character(0))
  
  for (g in names(parts)) {             # "0", "1", "2"
    for (tp in names(parts[[g]])) {     # "Mut", "PTM", "Hybrid"
      df0 <- parts[[g]][[tp]]
      if (!nrow(df0)) next
      df1 <- dedup_max_cc(df0)
      df_top <- top_n_by_cc(df1, n = top_n)
      
      # Write Top list with clear name
      top_file <- fname_top(out_dir_top, g, tp, top_n)
      .write_xlsx_openxlsx(df_top, top_file)
      
      # Parse site strings to PPI x Cancer wide table
      wide_file <- fname_sites(out_dir_top, g, tp)
      parse_sites_to_wide(top_file, wide_file)
      
      # Build frequency matrices and remember paths
      if (tp == "Mut") {
        mat_file <- fname_matrix_m(out_dir_top, g)
        build_mut_matrix(wide_file, mut_rate_xlsx, mat_file)
        generated$Mut <- c(generated$Mut, mat_file)
      } else if (tp == "PTM") {
        mat_file <- fname_matrix_p(out_dir_top, g)
        build_ptm_matrix(wide_file, ptm_rate_xlsx, mat_file)
        generated$PTM <- c(generated$PTM, mat_file)
      } else if (tp == "Hybrid") {
        # For Hybrid clusters, build both mutation-based and ptm-based matrices
        mat_file_mut <- fname_matrix_hm(out_dir_top, g)
        mat_file_ptm <- fname_matrix_hp(out_dir_top, g)
        build_mut_matrix(wide_file, mut_rate_xlsx, mat_file_mut)
        build_ptm_matrix(wide_file, ptm_rate_xlsx, mat_file_ptm)
        generated$Hybrid_Mut <- c(generated$Hybrid_Mut, mat_file_mut)
        generated$Hybrid_PTM <- c(generated$Hybrid_PTM, mat_file_ptm)
      }
    }
  }
  
  # Draw heatmaps -> PPT with clear names
  # Mutation (Reds)
  for (m in generated$Mut) {
    g <- sub("^.*geneCount(\\d+)\\.xlsx$", "\\1", m)
    save_heatmap_to_ppt(m, fname_ppt_m(out_dir_ppt, g), palette = "Reds")
  }
  # PTM (Blues)
  for (m in generated$PTM) {
    g <- sub("^.*geneCount(\\d+)\\.xlsx$", "\\1", m)
    save_heatmap_to_ppt(m, fname_ppt_p(out_dir_ppt, g), palette = "Blues")
  }
  # Hybrid-Mutation (Reds)
  for (m in generated$Hybrid_Mut) {
    g <- sub("^.*geneCount(\\d+)\\.xlsx$", "\\1", m)
    save_heatmap_to_ppt(m, fname_ppt_hm(out_dir_ppt, g), palette = "Reds")
  }
  # Hybrid-PTM (Blues)
  for (m in generated$Hybrid_PTM) {
    g <- sub("^.*geneCount(\\d+)\\.xlsx$", "\\1", m)
    save_heatmap_to_ppt(m, fname_ppt_hp(out_dir_ppt, g), palette = "Blues")
  }
  
  invisible(generated)
}

# ================================ Example run =================================
# Adjust the paths below to your environment, then run.
# (You can uncomment the setwd() if you prefer, but it’s not required.)

# setwd("F:/软著/cluster软著示例/数据可视化/Heatmap")

run_heatmap_pipeline(
  cluster_summary_xlsx = "F:/软著/cluster软著示例/数据可视化/cluster_summary_cc30.xlsx",
  out_dir_top          = "F:/软著/cluster软著示例/数据可视化/Heatmap/output_tables",
  mut_rate_xlsx        = "F:/软著/cluster软著示例/数据可视化/Heatmap/Mut_Freq_Pan_cancer.xlsx",
  ptm_rate_xlsx        = "F:/软著/cluster软著示例/数据可视化/Heatmap/PTM_Freq_Pan_9_cancer.xlsx",
  out_dir_ppt          = "F:/软著/cluster软著示例/数据可视化/Heatmap/ppt",
  top_n                = 30
)
