# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(readxl)
library(writexl)

# ------------------------------------------------------------------------------
# 1) Merge all per-structure TXT tables in a folder
# ------------------------------------------------------------------------------
merge_cluster_txt <- function(input_dir,
                              pattern = "\\.txt$",
                              out_csv = NULL) {
  # List files
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    warning("No TXT files found in: ", input_dir)
    return(invisible(NULL))
  }
  
  # Read non-empty tables
  dt_list <- lapply(files, function(f) {
    x <- tryCatch(fread(f, header = TRUE), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x
  })
  dt_list <- Filter(Negate(is.null), dt_list)
  if (length(dt_list) == 0) {
    warning("All TXT files are empty or unreadable in: ", input_dir)
    return(invisible(NULL))
  }
  
  merged <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  
  # Optional: write merged CSV
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    fwrite(merged, out_csv)
  }
  merged
}

# ------------------------------------------------------------------------------
# 2) Assign stable cluster IDs within each PDB
#    Cluster_New format: "Cluster.<PDB_ID>.<Cluster_ID>"
# ------------------------------------------------------------------------------
assign_cluster_ids <- function(df_merged) {
  stopifnot(all(c("PDB", "Cluster") %in% names(df_merged)))
  
  df <- as_tibble(df_merged)
  
  # PDB_ID: dense integer per unique PDB
  # Cluster_ID: dense integer per unique Cluster within each PDB, zero-based
  df2 <- df %>%
    mutate(PDB_ID = match(PDB, unique(PDB))) %>%
    group_by(PDB) %>%
    mutate(Cluster_ID = match(Cluster, unique(Cluster)) - 1L) %>%
    ungroup() %>%
    mutate(Cluster_New = paste0("Cluster.", PDB_ID, ".", Cluster_ID)) %>%
    select(-PDB_ID, -Cluster_ID)
  
  df2
}

# ------------------------------------------------------------------------------
# 3) Load PDB meta info (protein/gene etc.) and reduce to useful distinct columns
# ------------------------------------------------------------------------------
load_pdb_info <- function(pdb_inf_xlsx) {
  stopifnot(file.exists(pdb_inf_xlsx))
  PDB_Inf <- read_excel(pdb_inf_xlsx)
  useful_Inf <- PDB_Inf %>%
    select(PDB_ID, Uniprot1, Uniprot2, Gene1, Gene2, Chain1, Chain2, length1, length2, Resolution) %>%
    distinct()
  useful_Inf
}

# ------------------------------------------------------------------------------
# 4) Merge clusters with PDB info and site annotations (cancer/type)
#    Produces:
#      - Single_Site_Data_Raw.xlsx
#      - Single_Site_Data.xlsx (filtered columns)
#      - Cluster_Raw.xlsx (collapsed by Cluster_New)
# ------------------------------------------------------------------------------
build_outputs <- function(cluster_with_ids,
                          pdb_useful_inf,
                          site_all_txt,
                          out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Drop original Cluster (we use Cluster_New)
  cluster_data <- cluster_with_ids %>%
    select(-Cluster)
  
  # Merge PDB meta by PDB vs PDB_ID
  merged_Inf <- merge(cluster_data, pdb_useful_inf,
                      by.x = "PDB", by.y = "PDB_ID", all.x = TRUE)
  
  # Load site annotations (Gene/Type/Cancer_Name by Site_Inf)
  stopifnot(file.exists(site_all_txt))
  site_cancer <- fread(site_all_txt, sep = "\t", header = TRUE)
  site_cancer_selected <- site_cancer %>%
    select(Gene, Site_Inf, Type, Cancer_Name)
  
  # Merge
  all_sites_Inf <- merge(merged_Inf, site_cancer_selected,
                         by = "Site_Inf", all.x = TRUE)
  
  # Order/select columns (keep names consistent with your originals)
  # Note: your original used "C.vi." and "Cc" column names; keep them if present.
  cols_existing <- intersect(
    c("PDB","Uniprot1","Uniprot2","Gene1","Gene2","Chain1","Chain2",
      "length1","length2","Resolution","Gene","Chain","Site_Num",
      "Site_Inf","Cluster_New","C.vi.","Cc","Type","Cancer_Name"),
    names(all_sites_Inf)
  )
  all_sites_Inf <- all_sites_Inf %>% select(all_of(cols_existing))
  
  # Create Site_New for display
  if (all(c("Site_Inf","Cancer_Name") %in% names(all_sites_Inf))) {
    all_sites_Inf$Site_New <- paste0(all_sites_Inf$Site_Inf, "(", all_sites_Inf$Cancer_Name, ")")
  } else {
    all_sites_Inf$Site_New <- all_sites_Inf$Site_Inf
  }
  
  # Write single-site raw
  write_xlsx(all_sites_Inf, file.path(out_dir, "Single_Site_Data_Raw.xlsx"))
  
  # # Filtered single-site view (drop some columns if present)
  # filtered_data <- all_sites_Inf %>%
  #   select(-any_of(c("Site_Inf", "Site_Num", "Cancer_Name", "C.vi.", "Chain")))
  # write_xlsx(filtered_data, file.path(out_dir, "Single_Site_Data.xlsx"))
  
  # Collapse into per-cluster records:
  # Keep first values for stable metadata; merge unique Type/Site_New with commas
  final_cluster_data <- filtered_data %>%
    group_by(Cluster_New) %>%
    summarise(
      Type     = paste(sort(unique(na.omit(Type))), collapse = ","),
      Site_New = paste(sort(unique(na.omit(Site_New))), collapse = ","),
      across(.cols = -c(Type, Site_New, Cluster_New), .fns = ~ dplyr::first(., default = NA)),
      .groups = "drop"
    )
  write_xlsx(final_cluster_data, file.path(out_dir, "Cluster_Raw.xlsx"))
  
  list(
    single_site_raw = all_sites_Inf,
    single_site     = filtered_data,
    cluster_raw     = final_cluster_data
  )
}

# ------------------------------------------------------------------------------
# 5) Summary stats: count unique PPIs and unique structures
#    - PPI counted as unordered gene pair (Gene1, Gene2)
# ------------------------------------------------------------------------------
summarize_ppi_and_structs <- function(cluster_raw_xlsx) {
  stopifnot(file.exists(cluster_raw_xlsx))
  df <- read_excel(cluster_raw_xlsx)
  
  # Normalize unordered gene pairs
  df2 <- df %>%
    mutate(
      gA = pmin(Gene1, Gene2, na.rm = TRUE),
      gB = pmax(Gene1, Gene2, na.rm = TRUE)
    ) %>%
    distinct(gA, gB)
  
  n_ppi <- nrow(df2)
  n_structs <- length(unique(df$PDB))
  
  list(n_ppi = n_ppi, n_structs = n_structs)
}

# ------------------------------------------------------------------------------
# Orchestrator: one-call end-to-end
# ------------------------------------------------------------------------------
run_all <- function(
    cluster_txt_dir,                       # e.g. "F:/Data_0822/clusters/cluster_result_new"
    merged_csv_out,                        # e.g. "F:/Data_0822/clusters/cluster_results_merged_tables/Cluster_Raw_Merged.csv"
    cluster_raw_id_csv_out = "Cluster_Raw_Merged_ID.csv",
    pdb_info_xlsx,                         # e.g. "F:/Data_0822/clusters/cluster_results_merged_tables/All_Cluster_Sites_Inf_PPI_Database.xlsx"
    site_all_txt,                          # e.g. "F:/Data_0822/clusters/cluster_results_merged_tables/All_Cluster_Sites.txt"
    outputs_dir                            # e.g. "F:/Data_0822/clusters/cluster_results_merged_tables"
) {
  # 1) Merge
  merged <- merge_cluster_txt(cluster_txt_dir, out_csv = merged_csv_out)
  if (is.null(merged)) stop("No merged data produced. Aborting.")
  
  # 2) Assign cluster IDs
  merged_id <- assign_cluster_ids(merged)
  fwrite(merged_id, file.path(outputs_dir, cluster_raw_id_csv_out))
  
  # 3) Load PDB info
  pdb_useful <- load_pdb_info(pdb_info_xlsx)
  
  # 4) Build outputs
  outs <- build_outputs(
    cluster_with_ids = merged_id,
    pdb_useful_inf   = pdb_useful,
    site_all_txt     = site_all_txt,
    out_dir          = outputs_dir
  )
  
  # 5) Summary stats
  stats <- summarize_ppi_and_structs(file.path(outputs_dir, "Cluster_Raw.xlsx"))
  message(sprintf("Unique PPIs: %d | Unique PDB structures: %d", stats$n_ppi, stats$n_structs))
  
  invisible(list(data = outs, stats = stats))
}

# ------------------------------------------------------------------------------
# Example usage (uncomment and edit paths)
# ------------------------------------------------------------------------------
# run_all(
#   cluster_txt_dir = "F:/Data_0822/clusters/cluster_result_new",
#   merged_csv_out  = "F:/Data_0822/clusters/cluster_results_merged_tables/Cluster_Raw_Merged.csv",
#   cluster_raw_id_csv_out = "Cluster_Raw_Merged_ID.csv",
#   pdb_info_xlsx   = "F:/Data_0822/clusters/cluster_results_merged_tables/All_Cluster_Sites_Inf_PPI_Database.xlsx",
#   site_all_txt    = "F:/Data_0822/clusters/cluster_results_merged_tables/All_Cluster_Sites.txt",
#   outputs_dir     = "F:/Data_0822/clusters/cluster_results_merged_tables"
# )
