# ==============================================================================
# Cluster pipeline
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(writexl)
  library(stringr)
  library(purrr)
})

# ----------------------- logging & small helpers ------------------------------

options(dplyr.summarise.inform = FALSE)
options(cluster.pipeline.verbose = TRUE)

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && !(is.character(a) && a == "")) a else b

log_msg <- function(...){
  if (isTRUE(getOption("cluster.pipeline.verbose", TRUE))) {
    message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " ")))
    flush.console()
  }
}

time_it <- function(label, expr) {
  log_msg("START:", label)
  t0 <- proc.time()[["elapsed"]]
  on.exit({
    t1 <- proc.time()[["elapsed"]]
    log_msg("DONE :", label, "| elapsed:", sprintf("%.2fs", t1 - t0))
  }, add = TRUE)
  force(expr)
}

# ----------------------- IO helper -------------------------------------------

# Read xlsx/csv/tsv/txt into data.frame with logging and basic checks.
.read_tab <- function(path, sep = NULL) {
  log_msg("Reading file:", path)
  path_norm <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
  if (!file.exists(path_norm)) {
    stop(
      "File not found: ", path, "\n",
      "Normalized:     ", path_norm, "\n",
      "Working dir:    ", getwd()
    )
  }
  finfo <- tryCatch(file.info(path_norm), error = function(e) NULL)
  if (!is.null(finfo)) log_msg("Size:", sprintf("%.1f MB", finfo$size / (1024^2)))
  ext <- tolower(tools::file_ext(path_norm))
  out <- switch(
    ext,
    "xls"  = as.data.frame(readxl::read_excel(path_norm)),
    "xlsx" = as.data.frame(readxl::read_excel(path_norm)),
    "csv"  = as.data.frame(data.table::fread(path_norm, fill = TRUE)),
    "tsv"  = as.data.frame(data.table::fread(path_norm, sep = "\t", fill = TRUE)),
    "txt"  = as.data.frame(data.table::fread(path_norm, sep = sep %||% "\t", fill = TRUE)),
    {
      if (!is.null(sep)) as.data.frame(data.table::fread(path_norm, sep = sep, fill = TRUE))
      else stop("Unsupported file type: ", ext)
    }
  )
  log_msg("Loaded:", nrow(out), "rows x", ncol(out), "cols")
  out
}

# ----------------------- A1: interface site annotation -----------------------

# Mark per-site IF_Interface_Site by comparing Site_Num with adj_Interface1/2.
annotate_interface_sites <- function(single_site_raw_path, pdb_structure_path) {
  cluster_raw  <- .read_tab(single_site_raw_path)
  pdb_struct   <- .read_tab(pdb_structure_path)
  
  pdb_inf <- pdb_struct |>
    dplyr::select(PDB_ID, adj_Interface1, adj_Interface2, structure1, structure2) |>
    mutate(
      adj_Interface1 = as.character(adj_Interface1),
      adj_Interface2 = as.character(adj_Interface2)
    )
  
  merged <- merge(cluster_raw, pdb_inf, by.x = "PDB", by.y = "PDB_ID", all.x = TRUE)
  
  out <- merged |>
    mutate(
      adj1 = purrr::map(stringr::str_extract_all(dplyr::coalesce(adj_Interface1, ""), "\\d+"), as.integer),
      adj2 = purrr::map(stringr::str_extract_all(dplyr::coalesce(adj_Interface2, ""), "\\d+"), as.integer),
      site = purrr::map(stringr::str_extract_all(dplyr::coalesce(as.character(Site_Num), ""), "\\d+"), as.integer)
    ) |>
    rowwise() |>
    mutate(
      IF_Interface_Site = as.integer(
        if (Gene == Gene1) any(site %in% adj1)
        else if (Gene == Gene2) any(site %in% adj2)
        else FALSE
      )
    ) |>
    ungroup() |>
    select(-adj1, -adj2, -site)
  
  out
}

# ----------------------- A2: BioLiP pockets (fast) ---------------------------

# Expand Corrected_Residues to (Uniprot, SiteNum) and join per site.
annotate_biolip_pockets <- function(df, biolip_path) {
  log_msg("Biolip fast: explode residues and join")
  pocket_data <- .read_tab(biolip_path)
  stopifnot(all(c("Uniprot", "Corrected_Residues", "Pocket_ID") %in% names(pocket_data)))
  
  df$Uniprot  <- with(df, dplyr::case_when(Gene == Gene1 ~ Uniprot1, Gene == Gene2 ~ Uniprot2, TRUE ~ NA_character_))
  df$Site_Num <- suppressWarnings(as.integer(df$Site_Num))
  wanted_uniprots <- unique(df$Uniprot)
  
  pocket_dt <- data.table::as.data.table(pocket_data)[Uniprot %in% wanted_uniprots]
  pocket_long <- pocket_dt[
    , .(SiteNum = as.integer(unlist(stringr::str_extract_all(Corrected_Residues, "\\d+")))),
    by = .(Uniprot, Pocket_ID)
  ][!is.na(SiteNum)]
  
  df_dt <- data.table::as.data.table(df)
  df_dt[, rid := .I]
  data.table::setkey(df_dt, Uniprot, Site_Num)
  data.table::setkey(pocket_long, Uniprot, SiteNum)
  
  joined <- pocket_long[df_dt, on = .(Uniprot, SiteNum = Site_Num)]
  agg <- joined[
    , .(Pocket_IDs = if (all(is.na(Pocket_ID))) NA_character_ else paste(unique(Pocket_ID[!is.na(Pocket_ID)]), collapse = ",")),
    by = rid
  ]
  
  out <- merge(df_dt, agg, by = "rid", all.x = TRUE, sort = FALSE)[, rid := NULL]
  out <- out |>
    dplyr::mutate(Site_Name = stringr::str_sub(stringr::str_split_fixed(dplyr::coalesce(Site_New, ""), "-", 2)[, 2], 1, 1))
  
  as.data.frame(out)
}

# ----------------------- A3: fpocket IDs (fast) ------------------------------

# Tokenize fpocket Site_one and join with per-site token "Chain:ResidueNameResidueNum".
annotate_fpocket <- function(df, fpocket_summary_path) {
  log_msg("fpocket fast: explode tokens and join")
  fpocket <- .read_tab(fpocket_summary_path)
  stopifnot(all(c("pdb_id", "Fpocket_Name", "Site_one") %in% names(fpocket)))
  
  needed_pdbs <- unique(tolower(df$PDB))
  fpk <- fpocket[tolower(fpocket$pdb_id) %in% needed_pdbs, c("pdb_id","Fpocket_Name","Site_one")]
  
  split_tokens <- function(x) {
    if (is.na(x) || x == "") return(character(0))
    toks <- unlist(strsplit(x, "[,;\\s]+"))
    toks <- trimws(toks)
    toks[toks != ""]
  }
  
  fpk_dt <- data.table::as.data.table(fpk)
  fpk_long <- fpk_dt[, .(token = split_tokens(Site_one)), by = .(pdb_id, Fpocket_Name)]
  fpk_long <- fpk_long[grepl("^[A-Za-z]:[A-Z][0-9]+$", token)]
  data.table::setkey(fpk_long, pdb_id, token)
  
  df_dt <- data.table::as.data.table(df)
  df_dt[, rid := .I]
  df_dt[, token := paste0(Chain, ":", Site_Name, Site_Num)]
  df_dt[, pdb_id := tolower(PDB)]
  data.table::setkey(df_dt, pdb_id, token)
  
  joined <- fpk_long[df_dt, on = .(pdb_id, token)]
  agg <- joined[
    , .(Fpocket_ID = if (all(is.na(Fpocket_Name))) NA_character_ else paste(unique(Fpocket_Name[!is.na(Fpocket_Name)]), collapse = ",")),
    by = rid
  ]
  
  out <- merge(df_dt, agg, by = "rid", all.x = TRUE, sort = FALSE)[, c("rid","token","pdb_id") := NULL]
  as.data.frame(out)
}

# ----------------------- A4: functional site lists ---------------------------

# Tag per-site functional lists by intersecting tokens with PTM/mutation tables.
annotate_functional_sites <- function(df, ptm_sites_path, mut_sites_path) {
  ptm_sites <- read_excel(ptm_sites_path, skip = 1)
  mut_sites <- read_excel(mut_sites_path, skip = 1)
  
  df2 <- df %>%
    rowwise() %>%
    mutate(
      Fun_PTM = {
        tokens <- Site_Inf
        hits   <- intersect(tokens, ptm_sites$PTM_Site)
        if (length(hits) > 0) paste(hits, collapse = ",") else NA_character_
      },
      Fun_Mut = {
        tokens <- Site_Inf
        hits   <- intersect(tokens, mut_sites$gene_protein)
        if (length(hits) > 0) paste(hits, collapse = ",") else NA_character_
      }
    ) %>%
    ungroup()
  
  df2
}

# ----------------------- A5: pocket strings ----------------------------------

# Attach human-readable Pocket_Site / Fpocket_Site per site.
attach_pocket_strings <- function(df) {
  out <- df |>
    rowwise() |>
    mutate(
      Pocket_Site = if (is.na(Pocket_IDs) || is.na(Site_New)) NA_character_ else {
        paste(paste0(unique(str_split(Pocket_IDs, ",")[[1]]), "(", Site_New, ")"), collapse = ",")
      },
      Fpocket_Site = if (is.na(Fpocket_ID) || is.na(Site_New)) NA_character_ else {
        paste(paste0(unique(str_split(Fpocket_ID, ",")[[1]]), "(", Site_New, ")"), collapse = ",")
      }
    ) |>
    ungroup()
  out
}

# ----------------------- A6: single-site minimal export ----------------------

# Keep essential columns for cluster-level aggregation.
make_cluster_inf_raw <- function(df) {
  keep <- c("PDB","Cluster_New","Uniprot1","Uniprot2","Gene1","Gene2","Chain1","Chain2",
            "length1","length2","Uniprot","Gene","Chain","Type","Cc","IF_Interface_Site",
            "Fun_PTM","Fun_Mut","Site_New","Pocket_Site","Fpocket_Site")
  keep <- intersect(keep, names(df))
  dplyr::select(df, all_of(keep))
}

# ----------------------- B1: per-cluster summary -----------------------------

# Append chain tag by gene; used to form gene-token site strings per cluster.
.append_chain_suffix_gene_tokens <- function(tokens_chr, g1, g2, ch1, ch2) {
  if (is.null(tokens_chr) || is.na(tokens_chr) || !nzchar(trimws(tokens_chr))) return(tokens_chr)
  toks <- trimws(unlist(strsplit(tokens_chr, ",", fixed = TRUE)))
  toks <- toks[nzchar(toks)]
  if (!length(toks)) return(tokens_chr)
  suffixed <- vapply(toks, function(tk){
    gene <- sub("-.*$", "", tk)
    ch <- if (gene == g1) ch1 else if (gene == g2) ch2 else ""
    if (nzchar(ch)) paste0(tk, "-", ch) else tk
  }, character(1))
  paste(unique(suffixed), collapse = ",")
}

# Aggregate per cluster and derive Type, Pos_Type, and site strings.
summarize_clusters <- function(single_df) {
  update_type <- function(site_inf_new) {
    if (is.na(site_inf_new) || site_inf_new == "") return("Mut")
    has_ptm <- grepl("-[a-z]+\\(", site_inf_new)
    has_mut <- grepl("[A-Z][0-9]+[A-Z]\\(", site_inf_new)
    if (has_ptm && has_mut) "Hybrid" else if (has_ptm) "PTM" else "Mut"
  }
  
  out <- single_df %>%
    mutate(Type = ifelse(is.na(Type) | Type == "", sapply(Site_New, update_type), Type)) %>%
    group_by(PDB, Cluster_New) %>%
    summarise(
      Uniprot1          = first(Uniprot1),
      Uniprot2          = first(Uniprot2),
      Gene1             = first(Gene1),
      Gene2             = first(Gene2),
      Chain1            = first(Chain1),
      Chain2            = first(Chain2),
      length1           = first(length1),
      length2           = first(length2),
      Cc                = suppressWarnings(as.numeric(first(Cc))),
      IF_Interface_Site = as.integer(any(IF_Interface_Site == 1, na.rm = TRUE)),
      Type = {
        has_ptm <- any(Type %in% c("PTM", "Hybrid"), na.rm = TRUE)
        has_mut <- any(Type %in% c("Mut", "Hybrid"), na.rm = TRUE)
        dplyr::case_when(
          has_ptm && has_mut ~ "Hybrid",
          has_ptm            ~ "PTM",
          has_mut            ~ "Mut",
          TRUE               ~ NA_character_
        )
      },
      Fun_PTM      = paste(unique(na.omit(Fun_PTM)),      collapse = ","),
      Fun_Mut      = paste(unique(na.omit(Fun_Mut)),      collapse = ","),
      Site_Inf_New = paste(unique(na.omit(Site_New)),     collapse = ","),
      Pocket_Site  = paste(unique(na.omit(Pocket_Site)),  collapse = ","),
      Fpocket_Site = paste(unique(na.omit(Fpocket_Site)), collapse = ","),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(Site_New = .append_chain_suffix_gene_tokens(Site_Inf_New, Gene1, Gene2, Chain1, Chain2)) %>%
    ungroup()
  
  out
}

# ----------------------- B1.5: MD formatting string --------------------------

# Build MD site string "UniProt_Chain_X<pos>" and attach as Site_New_for_MD.
format_sites_for_md <- function(cluster_summary_df) {
  df <- cluster_summary_df |>
    mutate(Site_New_tmp = gsub("\\([^)]*\\)", "", Site_Inf_New))
  
  df1 <- df |>
    group_by(Cluster_New) |>
    summarise(
      Site_New_tmp = paste(unique(na.omit(Site_New_tmp)), collapse = ","),
      across(.cols = -c(Site_New_tmp), .fns = dplyr::first),
      .groups = "drop"
    )
  
  df1$Site_New_tmp <- mapply(function(s, g1, u1, g2, u2) {
    if (is.na(s) || s == "") return(s)
    sites <- trimws(unlist(strsplit(s, ",")))
    mapped <- vapply(sites, function(site) {
      gene <- sub("-.*$", "", site)
      gene_new <- if (gene == g1) u1 else if (gene == g2) u2 else gene
      sub(paste0("^", gene, "(?=-)"), gene_new, site, perl = TRUE)
    }, character(1))
    paste(mapped, collapse = ",")
  }, df1$Site_New_tmp, df1$Gene1, df1$Uniprot1, df1$Gene2, df1$Uniprot2, USE.NAMES = FALSE)
  
  convert_site <- function(site) {
    parts <- unlist(strsplit(site, "-"))
    if (length(parts) == 3) paste0(parts[1], "_", parts[2], "_X", gsub("[A-Za-z]", "", parts[3])) else site
  }
  df1$Site_New_tmp <- sapply(df1$Site_New_tmp, function(x) {
    if (is.na(x) || x == "") return(x)
    sites <- trimws(unlist(strsplit(x, ",")))
    paste(vapply(sites, convert_site, character(1)), collapse = ",")
  })
  
  df_out <- cluster_summary_df
  df_out$Site_New_for_MD <- df1$Site_New_tmp[match(df_out$Cluster_New, df1$Cluster_New)]
  df_out
}

# ----------------------- B2: COSMIC ------------------------------------------

# Add Cancer_gene by intersecting Gene1/Gene2 with COSMIC gene list.
add_cosmic <- function(cluster_summary_df, cosmic_csv_path) {
  cosmic <- read.csv(cosmic_csv_path)
  
  # Use fixed column name "Gene.Symbol" to extract gene list
  cosmic_genes <- unique(trimws(as.character(cosmic$Gene.Symbol)))
  
  df <- cluster_summary_df %>%
    rowwise() %>%
    mutate(
      Cancer_gene = {
        hits <- intersect(c(Gene1, Gene2), cosmic_genes)
        if (length(hits) == 0) NA_character_ else paste(hits, collapse = ",")
      }
    ) %>%
    ungroup()
  
  df
}


# ----------------------- B3: DrugBank ----------------------------------------

# Map DrugBank targets to Uniprot1/2 and provide gene names mapping.
add_drugbank <- function(cluster_df, drug_items_xlsx, proteins_csv, out_drugbank_inf_path = NULL) {
  data_item    <- readxl::read_excel(drug_items_xlsx)
  data_protein <- read.csv(proteins_csv, header = FALSE, stringsAsFactors = FALSE)
  
  pro_approved <- data_protein[data_protein$V1 %in% data_item$drugbank_id, ]
  merge_data <- merge(pro_approved, data_item, by.x = "V1", by.y = "drugbank_id", all.x = TRUE)
  
  if (!is.null(out_drugbank_inf_path)) writexl::write_xlsx(merge_data, out_drugbank_inf_path)
  
  drugbank_uniprots <- merge_data$V4
  
  cluster_data <- cluster_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Drugbank_gene = {
        hits <- intersect(c(Uniprot1, Uniprot2), drugbank_uniprots)
        if (length(hits) == 0) NA_character_ else paste(hits, collapse = ",")
      }
    ) %>%
    dplyr::ungroup()
  
  up2gene <- c(cluster_data$Gene1, cluster_data$Gene2)
  names(up2gene) <- c(cluster_data$Uniprot1, cluster_data$Uniprot2)
  
  map_ids <- function(id_string, dict) {
    if (is.na(id_string) || id_string == "") return("")
    tokens <- stringr::str_split(id_string, ",")[[1]]
    tokens <- trimws(tokens)
    mapped <- dict[tokens]
    mapped[is.na(mapped)] <- ""
    keep   <- mapped[mapped != ""]
    if (length(keep) == 0) "" else paste(unique(keep), collapse = ",")
  }
  
  cluster_out <- cluster_data %>%
    dplyr::mutate(Drugbank_gene_mapped = vapply(Drugbank_gene, map_ids, character(1), dict = up2gene))
  
  attr(cluster_out, "drugbank_inf") <- merge_data
  cluster_out
}

# ----------------------- B4: ASD mapping -------------------------------------

# Attach ASD_Target_Protein by PDB; also map to UniProt when possible.
add_asd_mapping <- function(df, asd_pairs_xlsx) {
  asd <- .read_tab(asd_pairs_xlsx)
  norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
  nms  <- names(asd); nml <- norm(nms)
  
  ap_idx <- which(nml == "allosteric_pdb")
  tg_idx <- which(nml == "target_gene")
  if (length(ap_idx) == 0 || length(tg_idx) == 0) {
    stop("ASD_HUMAN.xlsx must contain columns allosteric_pdb and target_gene.")
  }
  ap_col <- nms[ap_idx[1]]
  tg_col <- nms[tg_idx[1]]
  
  asd2 <- asd |>
    dplyr::mutate(PDB = tolower(.data[[ap_col]]), Target_Protein = .data[[tg_col]]) |>
    dplyr::group_by(PDB) |>
    dplyr::summarise(ASD_Target_Protein = paste(unique(Target_Protein), collapse = ","), .groups = "drop")
  
  out <- df |>
    dplyr::mutate(PDB_l = tolower(PDB)) |>
    dplyr::left_join(asd2, by = c("PDB_l" = "PDB")) |>
    dplyr::select(-PDB_l)
  
  # Keep also mapped UniProt names if needed downstream (not exported now).
  gene2up <- c(setNames(out$Uniprot1, out$Gene1), setNames(out$Uniprot2, out$Gene2))
  out$ASD_Target_Protein_New <- vapply(out$ASD_Target_Protein %||% NA_character_, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    toks <- trimws(unlist(strsplit(s, ",", fixed = TRUE)))
    m <- gene2up[toks]; m[is.na(m)] <- ""
    m <- unique(m[nzchar(m)])
    if (!length(m)) NA_character_ else paste(m, collapse = ",")
  }, character(1))
  
  out
}

# ----------------------- B5: build Site_New from InfRaw ----------------------

# Build per-cluster Site_New by concatenating "Site_New-Chain" from single-site table.
.compute_site_map_from_infraw_df <- function(single_site_df) {
  stopifnot(all(c("Cluster_New","Site_New","Chain") %in% names(single_site_df)))
  tmp <- single_site_df %>%
    mutate(
      Site_Inf_all = ifelse(
        is.na(Site_New) | is.na(Chain) | Site_New == "" | Chain == "",
        NA_character_,
        paste0(Site_New, "-", Chain)
      )
    ) %>%
    group_by(Cluster_New) %>%
    summarise(Site_New = paste(unique(na.omit(Site_Inf_all)), collapse = ","), .groups = "drop")
  tmp
}

# ----------------------- B6: apply Site_New + recompute Pos_Type  ------------

# Recompute Pos_Type using the requested rule.
# - If Site_New contains ≥2 distinct chain tags (e.g., "-A" and "-B"): "Cross-Interface"
# - Else (0/1 chain only): if any IF_Interface_Site==1 within that cluster: "Interface"; else "Non-Interface"
.apply_site_new_and_recompute_pos <- function(cluster_df, site_map_df, single_site_df) {
  # join Site_New from site_map_df
  out <- dplyr::left_join(
    cluster_df,
    dplyr::select(site_map_df, Cluster_New, Site_New),
    by = "Cluster_New",
    suffix = c("", ".inf")
  )
  if ("Site_New.inf" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(Site_New = dplyr::coalesce(.data$Site_New.inf, .data$Site_New)) %>%
      dplyr::select(-dplyr::all_of("Site_New.inf"))
  }
  
  # per-cluster interface flag aggregated from single-site table
  iface_map <- single_site_df %>%
    dplyr::group_by(Cluster_New) %>%
    dplyr::summarise(has_interface = any(IF_Interface_Site == 1, na.rm = TRUE), .groups = "drop")
  
  out <- dplyr::left_join(out, iface_map, by = "Cluster_New")
  
  # count distinct chains present in Site_New
  count_chains <- function(s) {
    if (is.null(s) || is.na(s) || !nzchar(s)) return(0L)
    toks <- trimws(unlist(strsplit(s, ",", fixed = TRUE)))
    toks <- toks[nzchar(toks)]
    if (!length(toks)) return(0L)
    tags <- stringr::str_extract_all(paste(toks, collapse=","), "-[A-Za-z](?=,|$)")[[1]]
    chains <- sub("^-", "", tags)
    length(unique(chains[!is.na(chains) & nzchar(chains)]))
  }
  
  out$.chains_n <- vapply(out$Site_New, count_chains, integer(1))
  
  out$Pos_Type <- ifelse(
    out$.chains_n >= 2,                         "Cross-Interface",
    ifelse(out$has_interface %||% FALSE,        "Interface", "Non-Interface")
  )
  
  out$.chains_n <- NULL
  out$has_interface[is.na(out$has_interface)] <- FALSE
  out
}

# ----------------------- C: MD merge with interface sites --------------------

# Merge MD site string with interface site collection; do not recompute Pos_Type here.
postprocess_md_sites <- function(md_sites_df,
                                 cluster_interface_sites_path,
                                 out_md_sites_path,
                                 out_pos_summary_path = NULL) {
  site_interface <- .read_tab(cluster_interface_sites_path)
  need <- c("PDB_ID","Interface_Site")
  if (!all(need %in% names(site_interface))) {
    stop("cluster_interface_sites.xlsx must contain columns: PDB_ID, Interface_Site.")
  }
  merged <- merge(md_sites_df, site_interface[, need], by.x = "PDB", by.y = "PDB_ID", all.x = TRUE)
  if (!is.null(out_md_sites_path)) writexl::write_xlsx(merged, out_md_sites_path)
  if (!is.null(out_pos_summary_path)) {
    # disabled per request (kept for optional debugging)
    writexl::write_xlsx(merged, out_pos_summary_path)
  }
  merged
}

# ----------------------- export column order ---------------------------------

# Define final export columns, removing Drugbank_gene and ASD_Target_Protein_New.
reorder_columns_as_example <- function(df) {
  wanted <- c(
    "PDB","Cluster_New","Uniprot1","Uniprot2","Gene1","Gene2","Chain1","Chain2",
    "length1","length2","Cc","Type","Pos_Type","Fun_PTM","Fun_Mut",
    "Pocket_Site","Fpocket_Site","Site_New",
    "Cancer_gene","Drugbank_gene_mapped","ASD_Target_Protein"
  )
  for (nm in wanted) if (!nm %in% names(df)) df[[nm]] <- NA
  df[, wanted]
}

# ----------------------- main pipeline ---------------------------------------

# End-to-end pipeline: builds single-site, cluster summary, MD sites; no Pos.xlsx.
run_pipeline <- function(
    single_site_raw_path,
    pdb_structure_path,
    biolip_path,
    fpocket_summary_path,
    ptm_sites_path,
    mut_sites_path,
    out_single_site_path,
    out_cluster_summary_path,
    out_cluster_sites_for_md_path,
    cosmic_csv_path,
    drug_items_xlsx,
    proteins_csv,
    asd_pairs_xlsx,
    cluster_interface_sites_path
) {
  # A) per-site enrichment
  df0 <- time_it("A1 annotate_interface_sites", {
    annotate_interface_sites(single_site_raw_path, pdb_structure_path)
  })
  df1 <- time_it("A2 annotate_biolip_pockets", {
    annotate_biolip_pockets(df0, biolip_path)
  })
  df2 <- time_it("A3 annotate_fpocket", {
    annotate_fpocket(df1, fpocket_summary_path)
  })
  df3 <- time_it("A4 annotate_functional_sites", {
    annotate_functional_sites(df2, ptm_sites_path, mut_sites_path)
  })
  df4 <- time_it("A5 attach_pocket_strings", {
    attach_pocket_strings(df3)
  })
  single_site_min <- time_it("A6 make_cluster_inf_raw", {
    make_cluster_inf_raw(df4)
  })
  if (!is.null(out_single_site_path)) {
    log_msg("Writing:", out_single_site_path)
    writexl::write_xlsx(single_site_min, out_single_site_path)
  }
  
  # B) cluster-level aggregation
  cluster_sum <- time_it("B1 summarize_clusters", {
    summarize_clusters(single_site_min)
  })
  cluster_sum <- time_it("B1.5 format_sites_for_md", {
    format_sites_for_md(cluster_sum)
  })
  if (!is.null(cosmic_csv_path)) {
    cluster_sum <- time_it("B2 add_cosmic", {
      add_cosmic(cluster_sum, cosmic_csv_path)
    })
  }
  if (!is.null(drug_items_xlsx) && !is.null(proteins_csv)) {
    cluster_sum <- time_it("B3 add_drugbank", {
      add_drugbank(cluster_sum, drug_items_xlsx, proteins_csv)
    })
  }
  if (!is.null(asd_pairs_xlsx)) {
    cluster_sum <- time_it("B4 add_asd_mapping", {
      add_asd_mapping(cluster_sum, asd_pairs_xlsx)
    })
  }
  
  # Build Site_New from per-site table, then recompute Pos_Type with the simple rule
  site_map_from_infraw <- time_it("B5 build Site_New from InfRaw", {
    .compute_site_map_from_infraw_df(single_site_min)
  })
  cluster_sum <- time_it("B6 apply Site_New & recompute Pos_Type", {
    .apply_site_new_and_recompute_pos(cluster_sum, site_map_from_infraw, single_site_min)
  })
  
  # Export
  cluster_sum_export <- reorder_columns_as_example(cluster_sum)
  if (!is.null(out_cluster_summary_path)) {
    log_msg("Writing:", out_cluster_summary_path)
    writexl::write_xlsx(cluster_sum_export, out_cluster_summary_path)
  }
  
  # MD helper file (no Pos_Type recompute here)
  md_sites <- cluster_sum[, c("PDB","Cluster_New","Uniprot1","Uniprot2","Gene1","Gene2","Chain1","Chain2","Site_New_for_MD")]
  names(md_sites)[names(md_sites)=="Site_New_for_MD"] <- "Site_New"
  invisible(postprocess_md_sites(
    md_sites_df                   = md_sites,
    cluster_interface_sites_path  = cluster_interface_sites_path,
    out_md_sites_path             = out_cluster_sites_for_md_path,
    out_pos_summary_path          = NULL
  ))
  
  log_msg("Pipeline finished.")
  invisible(list(
    single_site       = single_site_min,
    cluster_summary   = cluster_sum_export,
    md_sites          = md_sites
  ))
}

# ----------------------- example call ----------------------------------------
options(cluster.pipeline.verbose = TRUE)
setwd("F:/软著/cluster软著示例/")
res <- run_pipeline(
  single_site_raw_path          = "Single_Site_Data_Raw.xlsx",
  pdb_structure_path            = "F:/软著/cluster软著示例/supplement_data/PDB_Structure_Final.xlsx",
  biolip_path                   = "F:/软著/cluster软著示例/supplement_data/biolip_9606_with_PocketID.xlsx",
  fpocket_summary_path          = "F:/软著/cluster软著示例/supplement_data/fpocket_summary.xlsx",
  ptm_sites_path                = "F:/软著/cluster软著示例/supplement_data/merge_PTM_Inf_Final.xlsx",
  mut_sites_path                = "F:/软著/cluster软著示例/supplement_data/merge_mut_Inf_Final.xlsx",
  out_single_site_path          = "cluster_Inf_Raw.xlsx",
  out_cluster_summary_path      = "cluster_summary_final.xlsx",
  out_cluster_sites_for_md_path = "cluster_sites_for_MD.xlsx",
  cosmic_csv_path               = "F:/软著/cluster软著示例/supplement_data/Census_allFri May 30 06_55_30 2025.csv",
  drug_items_xlsx               = "F:/软著/cluster软著示例/supplement_data/Drug_data_approved.xlsx",
  proteins_csv                  = "F:/软著/cluster软著示例/supplement_data/items.csv",
  asd_pairs_xlsx                = "F:/软著/cluster软著示例/supplement_data/ASD_HUMAN.xlsx",
  cluster_interface_sites_path  = "F:/软著/cluster软著示例/supplement_data/cluster_interface_sites.xlsx"
)
