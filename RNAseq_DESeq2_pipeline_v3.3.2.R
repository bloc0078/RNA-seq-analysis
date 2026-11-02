# ==============================================================================
# SCRIPT: RNAseq_DESeq2_pipeline_v3.3.2.R
# VERSION: v3.3.2
# DESCRIPTION: Audit-safe DESeq2 pipeline with durable provenance and complete
#              output capture.
# AUTHOR: Gemini AI Collaboration, refactored by Codex
# ORIGINAL DATE: 2026-05-15
# VERSION DATE: 2026-07-24
#
# v3.3.2 AUDIT CHANGES:
# - Aligns DESeq2 independent filtering with the configured padj_cutoff by passing
#   alpha = padj_cutoff to results().
# - Clusters the sample-distance heatmap directly from its distance object rather
#   than calculating distances between distance profiles.
# - Rejects count values above R's integer range before integer conversion.
# - Corrects the script basename used for self-archiving and provenance.
# - Uses the USER INPUTS block as the sole configuration surface; input and
#   output roots are derived automatically from the R working directory.
#
# v3.3.0-v3.3.1 ROBUSTNESS CHANGES:
# - Organizes future runs by reference, project, contrast, and unique run ID.
# - Uses numbered, purpose-specific subfolders for inputs, results, QC, visuals,
#   plot data, serialized objects, logs, provenance, and checksums.
# - Copies the exact executed script and writes the complete manual run settings.
# - Records input and output SHA-256 checksums plus a final file inventory.
# - Saves every figure as both PNG and vector PDF.
# - Saves PCA scores/loadings/variance, volcano classifications and label choices,
#   MA plot data, histogram bins, heatmap ordering, and model diagnostic tables.
# - Saves named R objects explicitly instead of using save.image(), which archives
#   the global workspace rather than the local objects created inside main().
# - Records structured warnings and explicit RUN_COMPLETE/RUN_FAILED status files.
# - Improves input/configuration validation and captures dependency failures in
#   the run log after the run folder has been created.
# - Uses interpretable heatmaps with scales/condition annotations and accurate
#   adjusted-p-value volcano labeling; statistical calculations are unchanged.
#
# ==============================================================================

script_version <- "v3.3.2"
script_basename <- "RNAseq_DESeq2_pipeline_v3.3.2.R"

.deseq2_detect_script_path <- function(expected_basename) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  candidates <- character(0)

  if (length(file_arg) > 0) {
    candidates <- c(candidates, sub("^--file=", "", file_arg[1]))
  }

  for (frame_index in rev(seq_len(sys.nframe()))) {
    source_path <- tryCatch(sys.frame(frame_index)$ofile, error = function(e) NULL)
    if (!is.null(source_path) && length(source_path) == 1) {
      candidates <- c(candidates, source_path)
    }
  }

  candidates <- c(candidates, file.path(getwd(), expected_basename))
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    return(NA_character_)
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

.deseq2_checksum <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }

  unname(tools::md5sum(path))
}

.deseq2_checksum_algorithm <- function() {
  if (requireNamespace("digest", quietly = TRUE)) "SHA-256" else "MD5"
}

.deseq2_relative_path <- function(path, root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  normalized_root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", normalized_root), "/?"), "", normalized_path)
}

.deseq2_finalize_run <- function(context, warnings_log, run_error = NULL) {
  if (is.null(context$output_root) || !dir.exists(context$output_root)) {
    return(invisible(NULL))
  }

  end_time <- Sys.time()
  status <- if (is.null(run_error)) "COMPLETE" else "FAILED"
  error_message <- if (is.null(run_error)) "" else conditionMessage(run_error)

  warning_file <- file.path(context$paths$logs, paste0("Warnings_", context$run_id, ".csv"))
  if (length(warnings_log) == 0) {
    warning_table <- data.frame(
      timestamp = character(0),
      message = character(0),
      call = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    warning_table <- do.call(rbind, warnings_log)
  }
  utils::write.csv(warning_table, warning_file, row.names = FALSE)

  status_table <- data.frame(
    field = c(
      "status",
      "script_version",
      "run_id",
      "start_time_iso8601",
      "end_time_iso8601",
      "duration_seconds",
      "error_message"
    ),
    value = c(
      status,
      context$script_version,
      context$run_id,
      format(context$start_time, "%Y-%m-%dT%H:%M:%S%z"),
      format(end_time, "%Y-%m-%dT%H:%M:%S%z"),
      as.character(round(as.numeric(difftime(end_time, context$start_time, units = "secs")), 3)),
      error_message
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    status_table,
    file.path(context$paths$provenance, "Run_Status.csv"),
    row.names = FALSE
  )

  marker_name <- if (is.null(run_error)) "RUN_COMPLETE.txt" else "RUN_FAILED.txt"
  writeLines(
    c(
      paste0("status=", status),
      paste0("run_id=", context$run_id),
      paste0("script_version=", context$script_version),
      paste0("end_time=", format(end_time, "%Y-%m-%dT%H:%M:%S%z")),
      paste0("error_message=", error_message)
    ),
    file.path(context$paths$provenance, marker_name)
  )

  if (!is.null(context$log_filename)) {
    cat(
      "\n========================================================\n",
      "FINAL RUN STATUS: ", status, "\n",
      if (nzchar(error_message)) paste0("ERROR: ", error_message, "\n") else "",
      "End Time: ", format(end_time, "%Y-%m-%dT%H:%M:%S%z"), "\n",
      "========================================================\n",
      sep = "",
      file = context$log_filename,
      append = TRUE
    )
  }

  inventory_path <- file.path(context$paths$checksums, "File_Inventory.csv")
  checksum_name <- if (.deseq2_checksum_algorithm() == "SHA-256") "SHA256SUMS.txt" else "MD5SUMS.txt"
  checksum_path <- file.path(context$paths$checksums, checksum_name)

  all_files <- list.files(
    context$output_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
  all_files <- setdiff(normalizePath(all_files, winslash = "/", mustWork = TRUE),
                       normalizePath(c(inventory_path, checksum_path),
                                     winslash = "/", mustWork = FALSE))
  all_files <- sort(all_files)

  if (length(all_files) > 0) {
    file_info <- file.info(all_files)
    file_inventory <- data.frame(
      relative_path = vapply(
        all_files,
        .deseq2_relative_path,
        character(1),
        root = context$output_root
      ),
      folder_role = dirname(vapply(
        all_files,
        .deseq2_relative_path,
        character(1),
        root = context$output_root
      )),
      size_bytes = file_info$size,
      modified_time_iso8601 = format(file_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
      checksum_algorithm = .deseq2_checksum_algorithm(),
      checksum = vapply(all_files, .deseq2_checksum, character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    file_inventory <- data.frame(
      relative_path = character(0),
      folder_role = character(0),
      size_bytes = numeric(0),
      modified_time_iso8601 = character(0),
      checksum_algorithm = character(0),
      checksum = character(0),
      stringsAsFactors = FALSE
    )
  }

  utils::write.csv(file_inventory, inventory_path, row.names = FALSE)
  checksum_lines <- paste(file_inventory$checksum, file_inventory$relative_path, sep = "  ")
  writeLines(checksum_lines, checksum_path)

  invisible(file_inventory)
}

main <- function(context) {
  # ============================================================================
  # SECTION 1: LOAD LIBRARIES & DEFINE INPUTS
  # ============================================================================
  cat("\n========================================================\n")
  cat("SECTION 1: VALIDATING MANUAL INPUTS\n")
  cat("========================================================\n")

  safe_name <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(nchar(x) == 0, "unnamed", x)
  }

  require_exactly_one <- function(files, label) {
    if (length(files) != 1) {
      listing <- if (length(files) == 0) "  <none>" else paste0("  ", files, collapse = "\n")
      stop(
        paste0(
          "ERROR: Expected exactly one ", label, " file, found ", length(files), ".\n",
          "Matching files:\n", listing
        ),
        call. = FALSE
      )
    }
    files
  }

  annotate_results <- function(res_obj, gene_map) {
    res_df <- as.data.frame(res_obj)
    res_df$gene_id <- rownames(res_df)
    match_idx <- match(res_df$gene_id, gene_map$gene_id)
    res_df$gene <- gene_map$gene[match_idx]
    res_df$gene_label <- ifelse(!is.na(res_df$gene) & nzchar(res_df$gene), res_df$gene, res_df$gene_id)
    res_df[, c("gene_id", "gene", "gene_label", setdiff(names(res_df), c("gene_id", "gene", "gene_label")))]
  }

  annotate_vst <- function(vst_matrix, gene_map) {
    vst_df <- data.frame(gene_id = rownames(vst_matrix), vst_matrix, check.names = FALSE)
    match_idx <- match(vst_df$gene_id, gene_map$gene_id)
    vst_df$gene <- gene_map$gene[match_idx]
    vst_df$gene_label <- ifelse(!is.na(vst_df$gene) & nzchar(vst_df$gene), vst_df$gene, vst_df$gene_id)
    vst_df[, c("gene_id", "gene", "gene_label", setdiff(names(vst_df), c("gene_id", "gene", "gene_label")))]
  }

  annotate_gene_ids <- function(gene_ids, gene_map) {
    gene_df <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    match_idx <- match(gene_df$gene_id, gene_map$gene_id)
    gene_df$gene <- gene_map$gene[match_idx]
    gene_df$gene_label <- ifelse(!is.na(gene_df$gene) & nzchar(gene_df$gene), gene_df$gene, gene_df$gene_id)
    gene_df[, c("gene_id", "gene", "gene_label")]
  }

  save_png <- function(filename, width, height, res, expr) {
    expression_to_draw <- substitute(expr)
    drawing_environment <- parent.frame()
    pdf_filename <- file.path(
      paths$visuals_pdf,
      sub("\\.png$", ".pdf", basename(filename), ignore.case = TRUE)
    )

    png(filename, width = width, height = height, units = "in", res = res)
    png_result <- tryCatch(
      eval(expression_to_draw, envir = drawing_environment),
      finally = grDevices::dev.off()
    )

    grDevices::pdf(
      pdf_filename,
      width = width,
      height = height,
      useDingbats = FALSE,
      onefile = TRUE
    )
    tryCatch(
      eval(expression_to_draw, envir = drawing_environment),
      finally = grDevices::dev.off()
    )

    invisible(png_result)
  }

  # --- USER INPUTS ---
  reference_tag   <- "Erdman"
  ctrl_name      <- "ribA2_d4"
  test_name      <- "ribA2_atc_d4"
  ctrl_indices   <- c(33, 34, 35)
  test_indices   <- c(21, 22, 23)
  project_tag    <- "ribA2 Tet-OFF RNAseq"

  # --- THRESHOLDS ---
  padj_cutoff <- 0.05
  sig_lfc_cutoff <- 1.0
  volcano_lfc_cutoff <- 1.0
  top_n_volcano_labels <- 10

  # --- AUTO-DETECTED SETTINGS: DO NOT EDIT FOR A NORMAL RUN ---
  input_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  output_parent <- file.path(input_dir, "DESeq2_results")

  threshold_label <- function(prefix, value) {
    paste0(prefix, gsub("\\.", "p", as.character(value)))
  }

  padj_label <- threshold_label("FDR", padj_cutoff)
  sig_lfc_label <- threshold_label("LFC", sig_lfc_cutoff)
  volcano_lfc_label <- threshold_label("LFC", volcano_lfc_cutoff)

  validate_nonempty_scalar <- function(value, label) {
    if (length(value) != 1 || is.na(value) || !nzchar(trimws(as.character(value)))) {
      stop("ERROR: ", label, " must be one non-empty value.", call. = FALSE)
    }
  }

  validate_nonempty_scalar(reference_tag, "reference_tag")
  validate_nonempty_scalar(project_tag, "project_tag")
  validate_nonempty_scalar(ctrl_name, "ctrl_name")
  validate_nonempty_scalar(test_name, "test_name")
  validate_nonempty_scalar(input_dir, "input_dir")
  validate_nonempty_scalar(output_parent, "output_parent")

  if (identical(ctrl_name, test_name)) {
    stop("ERROR: ctrl_name and test_name must be different.", call. = FALSE)
  }
  if (length(ctrl_indices) == 0 || length(test_indices) == 0) {
    stop("ERROR: ctrl_indices and test_indices must each contain at least one column.", call. = FALSE)
  }
  if (!is.numeric(ctrl_indices) || !is.numeric(test_indices) ||
      any(!is.finite(c(ctrl_indices, test_indices))) ||
      any(c(ctrl_indices, test_indices) != as.integer(c(ctrl_indices, test_indices)))) {
    stop("ERROR: Sample indices must be finite whole numbers.", call. = FALSE)
  }
  if (!is.numeric(padj_cutoff) || length(padj_cutoff) != 1 ||
      is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
    stop("ERROR: padj_cutoff must be greater than 0 and no greater than 1.", call. = FALSE)
  }
  if (!is.numeric(sig_lfc_cutoff) || length(sig_lfc_cutoff) != 1 ||
      is.na(sig_lfc_cutoff) || sig_lfc_cutoff < 0) {
    stop("ERROR: sig_lfc_cutoff must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(volcano_lfc_cutoff) || length(volcano_lfc_cutoff) != 1 ||
      is.na(volcano_lfc_cutoff) || volcano_lfc_cutoff < 0) {
    stop("ERROR: volcano_lfc_cutoff must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(top_n_volcano_labels) || length(top_n_volcano_labels) != 1 ||
      is.na(top_n_volcano_labels) || top_n_volcano_labels < 0 ||
      top_n_volcano_labels != as.integer(top_n_volcano_labels)) {
    stop("ERROR: top_n_volcano_labels must be a non-negative whole number.", call. = FALSE)
  }

  # --- SAFE FILE NAMING ---
  reference_tag_safe <- safe_name(reference_tag)
  project_tag_safe <- safe_name(project_tag)
  ctrl_name_safe <- safe_name(ctrl_name)
  test_name_safe <- safe_name(test_name)
  contrast_tag <- paste0(test_name_safe, "_vs_", ctrl_name_safe)
  file_tag <- paste0(project_tag_safe, "_", contrast_tag)

  cat("[SUCCESS] Manual inputs and naming tags validated.\n")

  # ============================================================================
  # SECTION 2: INPUT DETECTION, DIRECTORY ARCHITECTURE & LOGGING
  # ============================================================================
  cat("\n========================================================\n")
  cat("SECTION 2: INITIALIZING SYSTEM LOGS & DIRECTORIES\n")
  cat("========================================================\n")

  start_time <- Sys.time()
  run_stamp <- format(start_time, "%Y%m%dT%H%M%S%z")
  run_id <- paste0("run_", run_stamp, "_pid", Sys.getpid())
  output_root <- file.path(
    output_parent,
    reference_tag_safe,
    project_tag_safe,
    contrast_tag,
    run_id
  )

  paths <- list(
    provenance = file.path(output_root, "00_provenance"),
    inputs = file.path(output_root, "01_inputs"),
    raw_counts = file.path(output_root, "02_results", "01_raw_counts"),
    normalized_counts = file.path(output_root, "02_results", "02_normalized_counts"),
    differential_expression = file.path(output_root, "02_results", "03_differential_expression"),
    transformed_expression = file.path(output_root, "02_results", "04_transformed_expression"),
    downstream = file.path(output_root, "02_results", "05_downstream_ready"),
    qc_tables = file.path(output_root, "03_qc", "tables"),
    visuals_png = file.path(output_root, "04_visuals", "png"),
    visuals_pdf = file.path(output_root, "04_visuals", "pdf"),
    plot_data = file.path(output_root, "05_plot_data"),
    objects = file.path(output_root, "06_serialized_objects"),
    logs = file.path(output_root, "07_logs"),
    checksums = file.path(output_root, "08_checksums")
  )

  dir_creation <- vapply(
    unique(unname(unlist(paths))),
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  )
  if (!all(dir.exists(unique(unname(unlist(paths)))))) {
    stop("ERROR: Failed to create the complete output folder structure.", call. = FALSE)
  }

  context$output_root <- output_root
  context$paths <- paths
  context$run_id <- run_id
  context$run_stamp <- run_stamp
  context$start_time <- start_time
  context$script_version <- script_version

  log_filename <- file.path(paths$logs, paste0("Execution_Log_", run_id, ".txt"))
  context$log_filename <- log_filename
  log_con <- file(log_filename, open = "wt")
  output_sink_start <- sink.number(type = "output")
  message_sink_start <- sink.number(type = "message")

  on.exit({
    while (sink.number(type = "message") > message_sink_start) {
      sink(type = "message")
    }
    while (sink.number(type = "output") > output_sink_start) {
      sink(type = "output")
    }
    if (!is.null(log_con) && isOpen(log_con)) {
      close(log_con)
    }
  }, add = TRUE)

  sink(log_con, type = "output", split = TRUE)
  sink(log_con, type = "message")

  cat("Creating Output Directory: ", output_root, "\n")

  suppressPackageStartupMessages({
    library("DESeq2")
    library("PCAtools")
    library("EnhancedVolcano")
    library("apeglm")
    library("rtracklayer")
    library("ggplot2")
    library("pheatmap")
    library("utils")
  })
  cat("[SUCCESS] Required libraries loaded.\n")

  script_path <- .deseq2_detect_script_path(script_basename)
  if (is.na(script_path)) {
    warning("The executing script path could not be determined; no script copy was archived.")
    archived_script_path <- NA_character_
    script_checksum <- NA_character_
  } else {
    archived_script_path <- file.path(paths$provenance, basename(script_path))
    copied_script <- file.copy(script_path, archived_script_path, overwrite = FALSE)
    if (!copied_script) {
      stop("ERROR: Failed to archive the exact executed script.", call. = FALSE)
    }
    script_checksum <- .deseq2_checksum(archived_script_path)
  }

  manual_inputs <- data.frame(
    field = c(
      "script_version",
      "reference_tag",
      "project_tag",
      "ctrl_name",
      "test_name",
      "ctrl_indices",
      "test_indices",
      "input_dir",
      "output_parent",
      "padj_cutoff",
      "sig_lfc_cutoff",
      "volcano_lfc_cutoff",
      "top_n_volcano_labels"
    ),
    value = c(
      script_version,
      reference_tag,
      project_tag,
      ctrl_name,
      test_name,
      paste(ctrl_indices, collapse = ","),
      paste(test_indices, collapse = ","),
      normalizePath(input_dir, winslash = "/", mustWork = FALSE),
      normalizePath(output_parent, winslash = "/", mustWork = FALSE),
      as.character(padj_cutoff),
      as.character(sig_lfc_cutoff),
      as.character(volcano_lfc_cutoff),
      as.character(top_n_volcano_labels)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    manual_inputs,
    file.path(paths$provenance, "Manual_Run_Settings.csv"),
    row.names = FALSE
  )

  run_readme <- c(
    "DESeq2 reproducible run bundle",
    paste0("Run ID: ", run_id),
    paste0("Script version: ", script_version),
    paste0("Reference: ", reference_tag),
    paste0("Project: ", project_tag),
    paste0("Contrast: ", test_name, " vs ", ctrl_name),
    "",
    "Folder map:",
    "00_provenance          exact script, manual settings, manifest, and run status",
    "01_inputs              immutable copies of the count table and GTF used",
    "02_results             raw/normalized counts, DE results, VST, and downstream exports",
    "03_qc/tables           numerical quality-control summaries",
    "04_visuals/png         high-resolution raster figures",
    "04_visuals/pdf         vector figures",
    "05_plot_data           numerical data and ordering behind every figure",
    "06_serialized_objects  named R objects needed to recreate results and figures",
    "07_logs                execution log, warnings, and session information",
    "08_checksums           final file inventory and integrity checksums",
    "",
    "Treat this run directory as immutable after completion."
  )
  writeLines(run_readme, file.path(paths$provenance, "README_RUN.txt"))

  found_counts <- list.files(path = input_dir, pattern = "^Master_.*\\.txt$", full.names = TRUE)
  found_gtf <- list.files(path = input_dir, pattern = "\\.gtf$", full.names = TRUE)
  input_counts <- require_exactly_one(found_counts, "Master_*.txt count")
  input_gtf <- require_exactly_one(found_gtf, ".gtf annotation")

  copied_counts_path <- file.path(paths$inputs, basename(input_counts))
  copied_gtf_path <- file.path(paths$inputs, basename(input_gtf))
  copied_counts <- file.copy(input_counts, copied_counts_path, overwrite = FALSE)
  copied_gtf <- file.copy(input_gtf, copied_gtf_path, overwrite = FALSE)
  if (!copied_counts || !copied_gtf) {
    stop("ERROR: Failed to copy selected input files into the audit input folder.", call. = FALSE)
  }

  input_provenance <- data.frame(
    role = c("count_table", "gtf_annotation"),
    original_path = normalizePath(c(input_counts, input_gtf), winslash = "/", mustWork = TRUE),
    archived_path = c(
      .deseq2_relative_path(copied_counts_path, output_root),
      .deseq2_relative_path(copied_gtf_path, output_root)
    ),
    size_bytes = file.info(c(input_counts, input_gtf))$size,
    checksum_algorithm = .deseq2_checksum_algorithm(),
    checksum = vapply(c(copied_counts_path, copied_gtf_path), .deseq2_checksum, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(
    input_provenance,
    file.path(paths$provenance, "Input_Provenance.csv"),
    row.names = FALSE
  )

  cat("========================================================\n")
  cat("DESeq2 PIPELINE REPORT\n")
  cat("========================================================\n")
  cat("Project Tag:      ", project_tag, "\n")
  cat("Reference Tag:    ", reference_tag, "\n")
  cat("Control Label:    ", ctrl_name, "\n")
  cat("Test Label:       ", test_name, "\n")
  cat("Run ID:           ", run_id, "\n")
  cat("Run Timestamp:    ", format(start_time, "%Y-%m-%dT%H:%M:%S%z"), "\n")
  cat("Script Version:   ", script_version, "\n")
  cat("Script SHA:       ", script_checksum, "\n")
  cat("File Tag Used:    ", file_tag, "\n")
  cat("Output Folder:    ", output_root, "\n")
  cat("Count File:       ", input_counts, "\n")
  cat("GTF File:         ", input_gtf, "\n")
  cat("--------------------------------------------------------\n")

  # ============================================================================
  # SECTION 3: DATA LOADING, GTF ANNOTATION & SAMPLE AUDIT
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 3: IMPORTING MATRICES, GTF MAPPING & SAMPLE AUDIT\n")
  cat("--------------------------------------------------------\n")
  
  gtf_obj <- rtracklayer::import(copied_gtf_path)
  gtf_df <- as.data.frame(gtf_obj)
  
  # featureCounts was run with:
  #   -t CDS -g locus_tag
  # Therefore the DESeq2 count row IDs should be matched to GTF locus_tag
  # from CDS rows, not only to type == "gene" rows.
  
  if ("type" %in% names(gtf_df) && any(gtf_df$type == "CDS", na.rm = TRUE)) {
    gtf_df_map <- gtf_df[gtf_df$type == "CDS", , drop = FALSE]
    cat("Using GTF feature type for annotation map: CDS\n")
  } else {
    gtf_df_map <- gtf_df
    cat("No CDS feature rows found; using all GTF rows for annotation map.\n")
  }
  
  id_candidates <- c("locus_tag", "gene_id", "gene", "Name")
  id_col <- id_candidates[id_candidates %in% names(gtf_df_map)][1]
  
  if (is.na(id_col)) {
    stop(
      "ERROR: Could not find a usable GTF ID column. Tried: ",
      paste(id_candidates, collapse = ", "),
      call. = FALSE
    )
  }
  
  gene_name_candidates <- c("gene", "gene_name", "Name", "locus_tag", "gene_id")
  gene_col <- gene_name_candidates[gene_name_candidates %in% names(gtf_df_map)][1]
  
  if (is.na(gene_col)) {
    cat("No gene-name-like GTF column found; using ", id_col, " as gene label.\n", sep = "")
    gene_values <- gtf_df_map[[id_col]]
  } else {
    cat("Using GTF column for gene IDs: ", id_col, "\n", sep = "")
    cat("Using GTF column for gene names: ", gene_col, "\n", sep = "")
    gene_values <- gtf_df_map[[gene_col]]
  }
  
  gene_map <- unique(data.frame(
    gene_id = gtf_df_map[[id_col]],
    gene = gene_values,
    stringsAsFactors = FALSE
  ))
  
  gene_map <- gene_map[!is.na(gene_map$gene_id) & nzchar(gene_map$gene_id), , drop = FALSE]
  
  duplicated_gene_ids <- gene_map$gene_id[duplicated(gene_map$gene_id)]
  
  if (length(duplicated_gene_ids) > 0) {
    duplicated_gene_ids <- unique(duplicated_gene_ids)
    
    cat("[WARNING] Duplicate gene_id entries found in gene_map: ",
        length(duplicated_gene_ids),
        " duplicated gene IDs. Keeping first occurrence.\n",
        sep = "")
    
    write.csv(
      gene_map[gene_map$gene_id %in% duplicated_gene_ids, ],
      file.path(paths$qc_tables, paste0("Duplicate_GTF_Gene_IDs_", file_tag, ".csv")),
      row.names = FALSE
    )
    
    gene_map <- gene_map[!duplicated(gene_map$gene_id), , drop = FALSE]
  }
  
  cat("[REPORT] Gene map rows: ", nrow(gene_map), "\n", sep = "") 

# Manual label fixes specific to this ribA2 Tet-OFF dataset.
# Do not include older coaX/panG relabeling from prior inherited scripts;
# those were project-specific carryovers and are not appropriate here.
manual_gene_fixes <- c()

applied_fixes <- character(0)
for (fix_id in names(manual_gene_fixes)) {
  fix_idx <- which(gene_map$gene_id == fix_id)
  if (length(fix_idx) > 0) {
    gene_map$gene[fix_idx] <- manual_gene_fixes[[fix_id]]
    applied_fixes <- c(applied_fixes, paste0(fix_id, "=", manual_gene_fixes[[fix_id]]))
  }
}
  if (length(applied_fixes) > 0) {
    cat("Manual gene label fixes applied: ", paste(applied_fixes, collapse = ", "), "\n")
  } else {
    cat("Manual gene label fixes applied: none matched gene_map.\n")
  }

  raw_counts <- read.table(copied_counts_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  sample_indices <- c(ctrl_indices, test_indices)

  if (any(sample_indices < 1)) {
    stop("ERROR: Sample indices must be >= 1.", call. = FALSE)
  }
  if (any(duplicated(sample_indices))) {
    stop("ERROR: Duplicate sample indices detected: ", paste(sample_indices[duplicated(sample_indices)], collapse = ", "), call. = FALSE)
  }
  if (any(sample_indices > ncol(raw_counts))) {
    stop(
      "ERROR: Sample index greater than available count columns. Max requested: ",
      max(sample_indices),
      "; available columns: ",
      ncol(raw_counts),
      call. = FALSE
    )
  }

  renamed_sample_ids <- c(
    paste0(ctrl_name, "_", seq_along(ctrl_indices)),
    paste0(test_name, "_", seq_along(test_indices))
  )
  sample_conditions <- c(
    rep(ctrl_name, length(ctrl_indices)),
    rep(test_name, length(test_indices))
  )
  sample_mapping <- data.frame(
    original_column_index = sample_indices,
    original_column_name = colnames(raw_counts)[sample_indices],
    renamed_sample_id = renamed_sample_ids,
    condition = sample_conditions,
    stringsAsFactors = FALSE
  )
  sample_mapping_filename <- paste0("Sample_Column_Mapping_", file_tag, ".csv")
  write.csv(sample_mapping, file.path(paths$provenance, sample_mapping_filename), row.names = FALSE)
  cat("Sample column mapping:\n")
  print(sample_mapping)

  count_matrix <- as.matrix(raw_counts[, sample_indices, drop = FALSE])
  colnames(count_matrix) <- renamed_sample_ids
  numeric_count_matrix <- suppressWarnings(matrix(as.numeric(count_matrix), nrow = nrow(count_matrix), dimnames = dimnames(count_matrix)))
  if (anyNA(numeric_count_matrix)) {
    stop("ERROR: Count matrix contains NA or non-numeric values.", call. = FALSE)
  }
  if (any(numeric_count_matrix < 0)) {
    stop("ERROR: Count matrix contains negative values.", call. = FALSE)
  }
  if (any(numeric_count_matrix != round(numeric_count_matrix))) {
    stop("ERROR: Count matrix contains non-integer values.", call. = FALSE)
  }
  if (any(numeric_count_matrix > .Machine$integer.max)) {
    stop(
      "ERROR: Count matrix contains values above R's maximum integer value (",
      .Machine$integer.max,
      ").",
      call. = FALSE
    )
  }
  count_matrix <- numeric_count_matrix
  storage.mode(count_matrix) <- "integer"
  if (anyNA(count_matrix)) {
    stop("ERROR: Integer conversion introduced NA values.", call. = FALSE)
  }
  
  raw_selected_counts_filename <- paste0("Raw_Selected_Counts_", file_tag, ".csv")
  raw_selected_counts_export <- data.frame(
    gene_id = rownames(count_matrix),
    count_matrix,
    check.names = FALSE
  )
  
  write.csv(
    raw_selected_counts_export,
    file.path(paths$raw_counts, raw_selected_counts_filename),
    row.names = FALSE,
    quote = FALSE
  )

  raw_count_qc <- data.frame(
    sample_id = colnames(count_matrix),
    condition = sample_conditions,
    raw_library_size = colSums(count_matrix),
    detected_genes_count_gt0 = colSums(count_matrix > 0),
    detected_genes_count_ge10 = colSums(count_matrix >= 10),
    percent_zero_counts = 100 * colSums(count_matrix == 0) / nrow(count_matrix),
    row.names = NULL
  )
  
  write.csv(raw_count_qc, file.path(paths$qc_tables, paste0("Raw_Count_QC_", file_tag, ".csv")), row.names = FALSE)
  cat("[QC] Raw library size min/max: ",
      min(raw_count_qc$raw_library_size), " / ", max(raw_count_qc$raw_library_size), "\n", sep = "")
  cat("[QC] Detected genes >0 min/max: ",
      min(raw_count_qc$detected_genes_count_gt0), " / ", max(raw_count_qc$detected_genes_count_gt0), "\n", sep = "")
  cat("[QC] Detected genes >=10 min/max: ",
      min(raw_count_qc$detected_genes_count_ge10), " / ", max(raw_count_qc$detected_genes_count_ge10), "\n", sep = "")

  # ============================================================================
  # SECTION 4: RUNNING DESEQ2 STATISTICS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 4: RUNNING DESEQ2 STATISTICS\n")
  cat("--------------------------------------------------------\n")

  condition_factor <- factor(sample_conditions, levels = c(ctrl_name, test_name))
  sampleTable <- data.frame(
    sample_id = colnames(count_matrix),
    condition = condition_factor,
    row.names = colnames(count_matrix)
  )
  
  cat("[REPORT] Sample counts by condition:\n")
  print(table(sampleTable$condition))
  
  if (any(table(sampleTable$condition) < 2)) {
    warning(
      "One or more conditions has fewer than 2 samples. DESeq2 can run, but replicate-based inference will be weak."
    )
  }

  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sampleTable, design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = ctrl_name)
  all_gene_ids_before_filter <- rownames(dds)
  genes_before_filter <- nrow(dds)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  kept_gene_ids <- rownames(dds)
  removed_gene_ids <- all_gene_ids_before_filter[!keep]
  genes_after_filter <- nrow(dds)
  genes_removed_filter <- length(removed_gene_ids)

  genes_kept_annotated <- annotate_gene_ids(kept_gene_ids, gene_map)
  genes_removed_annotated <- annotate_gene_ids(removed_gene_ids, gene_map)
  write.csv(
    genes_kept_annotated,
    file.path(paths$qc_tables, paste0("Genes_Kept_After_LowCount_Filter_", file_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    genes_removed_annotated,
    file.path(paths$qc_tables, paste0("Genes_Removed_By_LowCount_Filter_", file_tag, ".csv")),
    row.names = FALSE
  )

  cat("[REPORT] Genes before low-count filter: ", genes_before_filter, "\n", sep = "")
  cat("[REPORT] Genes after low-count filter:  ", genes_after_filter, "\n", sep = "")
  cat("[REPORT] Genes removed by low-count filter: ", genes_removed_filter, "\n", sep = "")

  dds <- DESeq(dds)

  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_export <- data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE)
  write.csv(norm_counts_export, file.path(paths$normalized_counts, paste0("Normalized_Counts_", file_tag, ".csv")), row.names = FALSE, quote = FALSE)

  normalized_count_qc <- data.frame(
    sample_id = colnames(norm_counts),
    condition = as.character(sampleTable[colnames(norm_counts), "condition"]),
    normalized_library_size = colSums(norm_counts),
    normalized_detected_genes_gt0 = colSums(norm_counts > 0),
    normalized_detected_genes_ge10 = colSums(norm_counts >= 10),
    row.names = NULL
  )
  write.csv(normalized_count_qc, file.path(paths$qc_tables, paste0("Normalized_Count_QC_", file_tag, ".csv")), row.names = FALSE)

  available_coefs <- resultsNames(dds)
  coef_pattern <- paste0("condition_", test_name, "_vs_", ctrl_name)
  coef_matches <- grep(coef_pattern, available_coefs, fixed = TRUE, value = TRUE)
  if (length(coef_matches) != 1) {
    stop(
      paste0(
        "ERROR: Expected exactly one apeglm coefficient matching '", coef_pattern, "', found ",
        length(coef_matches), ".\nAvailable coefficients:\n",
        paste0("  ", available_coefs, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  target_coef <- coef_matches
  res_unshrunken <- results(
    dds,
    contrast = c("condition", test_name, ctrl_name),
    alpha = padj_cutoff
  )
  resLFC <- lfcShrink(dds, coef = target_coef, type = "apeglm")

  # ============================================================================
  # SECTION 5: EXPORTING DATA, METADATA & QUALITY METRICS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 5: EXPORTING DATA, METADATA & QUALITY METRICS\n")
  cat("--------------------------------------------------------\n")

  write.csv(sampleTable, file.path(paths$provenance, paste0("Sample_Metadata_", file_tag, ".csv")), row.names = FALSE)
  write.csv(gene_map, file.path(paths$provenance, paste0("Gene_Annotation_Map_", project_tag_safe, ".csv")), row.names = FALSE)

  res_unshrunken_annotated <- annotate_results(res_unshrunken, gene_map)
  res_shrunken_annotated <- annotate_results(resLFC, gene_map)

  id_matches <- res_shrunken_annotated$gene_id %in% gene_map$gene_id
  id_match_count <- sum(id_matches)
  annotation_total <- nrow(res_shrunken_annotated)
  id_match_rate <- 100 * id_match_count / annotation_total
  
  gene_symbol_count <- sum(!is.na(res_shrunken_annotated$gene) & nzchar(res_shrunken_annotated$gene))
  gene_symbol_rate <- 100 * gene_symbol_count / annotation_total
  
  cat("[REPORT] Gene ID matches: ",
      id_match_count, " / ", annotation_total,
      " (", round(id_match_rate, 1), "%)\n", sep = "")
  
  cat("[REPORT] Non-empty gene symbols: ",
      gene_symbol_count, " / ", annotation_total,
      " (", round(gene_symbol_rate, 1), "%)\n", sep = "")
  
  if (id_match_rate < 90) {
    warning(
      "Gene ID match rate is below 90%. Check whether count IDs and GTF locus_tag/gene_id values use the same ID system."
    )
  }
  
  unmatched_count_ids <- res_shrunken_annotated$gene_id[!id_matches]
  write.csv(
    data.frame(gene_id = unmatched_count_ids),
    file.path(paths$qc_tables, paste0("Unmatched_Count_IDs_Not_In_GTF_", file_tag, ".csv")),
    row.names = FALSE
  )

  missing_gene_symbol_ids <- res_shrunken_annotated$gene_id[
    is.na(res_shrunken_annotated$gene) | !nzchar(res_shrunken_annotated$gene)
  ]

  write.csv(
    data.frame(gene_id = missing_gene_symbol_ids),
    file.path(paths$qc_tables, paste0("Count_IDs_Missing_Gene_Symbol_", file_tag, ".csv")),
    row.names = FALSE
  )

  unshrunken_result_filename <- paste0("DESeq2_Unshrunken_Results_", file_tag, ".csv")
  shrunken_result_filename <- paste0("DESeq2_apeglm_Shrunken_Results_", file_tag, ".csv")
  write.csv(res_unshrunken_annotated, file.path(paths$differential_expression, unshrunken_result_filename), row.names = FALSE)
  write.csv(res_shrunken_annotated, file.path(paths$differential_expression, shrunken_result_filename), row.names = FALSE)

  crosswalk_gene_ids <- unique(c(gene_map$gene_id, res_shrunken_annotated$gene_id))
  crosswalk <- annotate_gene_ids(crosswalk_gene_ids, gene_map)
  crosswalk$original_count_id_present <- crosswalk$gene_id %in% rownames(raw_counts)
  write.csv(
    crosswalk,
    file.path(paths$provenance, paste0("Gene_ID_Crosswalk_", file_tag, ".csv")),
    row.names = FALSE
  )

  ma_unshrunken_data <- DESeq2::plotMA(
    res_unshrunken,
    alpha = padj_cutoff,
    returnData = TRUE
  )
  ma_unshrunken_data$gene_id <- rownames(res_unshrunken)
  ma_unshrunken_data <- ma_unshrunken_data[, c("gene_id", "mean", "lfc", "isDE")]
  write.csv(
    ma_unshrunken_data,
    file.path(paths$plot_data, paste0("MA_Unshrunken_Data_", file_tag, ".csv")),
    row.names = FALSE
  )

  ma_shrunken_data <- DESeq2::plotMA(
    resLFC,
    alpha = padj_cutoff,
    returnData = TRUE
  )
  ma_shrunken_data$gene_id <- rownames(resLFC)
  ma_shrunken_data <- ma_shrunken_data[, c("gene_id", "mean", "lfc", "isDE")]
  write.csv(
    ma_shrunken_data,
    file.path(paths$plot_data, paste0("MA_apeglm_Shrunken_Data_", file_tag, ".csv")),
    row.names = FALSE
  )

  save_png(file.path(paths$visuals_png, paste0("MA_Unshrunken_", file_tag, ".png")), 6, 6, 300, {
    plotMA(
      res_unshrunken,
      alpha = padj_cutoff,
      ylim = c(-5, 5),
      main = paste("MA Unshrunken:", project_tag)
    )
  })
  save_png(file.path(paths$visuals_png, paste0("MA_apeglm_Shrunken_", file_tag, ".png")), 6, 6, 300, {
    plotMA(
      resLFC,
      alpha = padj_cutoff,
      ylim = c(-5, 5),
      main = paste("MA apeglm Shrunken:", project_tag)
    )
  })

  pvalue_non_na <- res_unshrunken_annotated$pvalue[!is.na(res_unshrunken_annotated$pvalue)]
  pvalue_histogram <- if (length(pvalue_non_na) > 0) {
    hist(pvalue_non_na, breaks = 50, plot = FALSE)
  } else {
    NULL
  }
  if (!is.null(pvalue_histogram)) {
    pvalue_histogram_table <- data.frame(
      bin_lower = head(pvalue_histogram$breaks, -1),
      bin_upper = tail(pvalue_histogram$breaks, -1),
      midpoint = pvalue_histogram$mids,
      count = pvalue_histogram$counts,
      density = pvalue_histogram$density
    )
  } else {
    pvalue_histogram_table <- data.frame(
      bin_lower = numeric(0),
      bin_upper = numeric(0),
      midpoint = numeric(0),
      count = integer(0),
      density = numeric(0)
    )
  }
  write.csv(
    pvalue_histogram_table,
    file.path(paths$plot_data, paste0("PValue_Histogram_Bins_", file_tag, ".csv")),
    row.names = FALSE
  )

  save_png(file.path(paths$visuals_png, paste0("PValue_Histogram_Unshrunken_", file_tag, ".png")), 6, 5, 300, {
    if (length(pvalue_non_na) > 0) {
      plot(
        pvalue_histogram,
        main = paste("Unshrunken p-value Distribution:", project_tag),
        xlab = "p-value",
        col = "grey80",
        border = "white"
      )
    } else {
      plot.new()
      title(main = paste("Unshrunken p-value Distribution:", project_tag))
      text(0.5, 0.5, "No non-NA p-values available")
    }
  })
  pvalue_summary <- data.frame(
    n_total_genes = nrow(res_unshrunken_annotated),
    n_non_na_pvalue = sum(!is.na(res_unshrunken_annotated$pvalue)),
    n_na_pvalue = sum(is.na(res_unshrunken_annotated$pvalue)),
    n_non_na_padj = sum(!is.na(res_unshrunken_annotated$padj)),
    n_na_padj = sum(is.na(res_unshrunken_annotated$padj)),
    min_pvalue = if (length(pvalue_non_na) > 0) min(pvalue_non_na) else NA_real_,
    median_pvalue = if (length(pvalue_non_na) > 0) median(pvalue_non_na) else NA_real_,
    min_padj = if (any(!is.na(res_unshrunken_annotated$padj))) min(res_unshrunken_annotated$padj, na.rm = TRUE) else NA_real_,
    median_padj = if (any(!is.na(res_unshrunken_annotated$padj))) median(res_unshrunken_annotated$padj, na.rm = TRUE) else NA_real_
  )
  write.csv(pvalue_summary, file.path(paths$qc_tables, paste0("PValue_Summary_Unshrunken_", file_tag, ".csv")), row.names = FALSE)

  sf_df <- data.frame(Sample = names(sizeFactors(dds)), SizeFactor = sizeFactors(dds), row.names = NULL)
  write.csv(sf_df, file.path(paths$qc_tables, paste0("Size_Factors_", run_id, ".csv")), row.names = FALSE)

  model_diagnostics <- as.data.frame(S4Vectors::mcols(dds))
  model_diagnostics$gene_id <- rownames(dds)
  model_diagnostics <- model_diagnostics[, c("gene_id", setdiff(names(model_diagnostics), "gene_id"))]
  atomic_diagnostic_columns <- vapply(
    model_diagnostics,
    function(column) is.atomic(column) && is.null(dim(column)),
    logical(1)
  )
  model_diagnostics_export <- model_diagnostics[, atomic_diagnostic_columns, drop = FALSE]
  write.csv(
    model_diagnostics_export,
    file.path(paths$plot_data, paste0("DESeq2_Model_Diagnostics_", file_tag, ".csv")),
    row.names = FALSE
  )

  disp_plot_filename <- paste0("Dispersion_Plot_", file_tag, ".png")
  save_png(file.path(paths$visuals_png, disp_plot_filename), 6, 5, 300, {
    plotDispEsts(dds, main = paste("Dispersion Estimates:", project_tag))
  })

  # ============================================================================
  # SECTION 6: VST NORMALIZED COUNTS & GROUP MEANS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 6: GENERATING VST NORMALIZED COUNTS\n")
  cat("--------------------------------------------------------\n")

  vsd <- DESeq2::vst(dds, blind = FALSE)
  vst_matrix <- SummarizedExperiment::assay(vsd)

  sample_cor <- cor(vst_matrix, method = "pearson")
  write.csv(sample_cor, file.path(paths$qc_tables, paste0("Sample_Correlation_VST_", file_tag, ".csv")), quote = FALSE)
  sample_annotation <- data.frame(
    condition = as.character(sampleTable[colnames(vst_matrix), "condition"]),
    row.names = colnames(vst_matrix),
    stringsAsFactors = FALSE
  )
  condition_colors <- setNames(c("#1f1f1f", "#7B2CBF"), c(ctrl_name, test_name))
  annotation_colors <- list(condition = condition_colors)
  correlation_hclust <- hclust(dist(sample_cor))
  correlation_heatmap <- pheatmap::pheatmap(
    sample_cor,
    cluster_rows = correlation_hclust,
    cluster_cols = correlation_hclust,
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    border_color = NA,
    main = "VST Pearson Sample Correlation",
    silent = TRUE
  )
  save_png(file.path(paths$visuals_png, paste0("Sample_Correlation_Heatmap_VST_", file_tag, ".png")), 7, 6, 300, {
    grid::grid.newpage()
    grid::grid.draw(correlation_heatmap$gtable)
  })
  correlation_order <- data.frame(
    position = seq_along(correlation_hclust$order),
    sample_id = labels(correlation_hclust)[correlation_hclust$order],
    stringsAsFactors = FALSE
  )
  write.csv(
    correlation_order,
    file.path(paths$plot_data, paste0("Sample_Correlation_Heatmap_Order_", file_tag, ".csv")),
    row.names = FALSE
  )

  sample_dist <- as.matrix(dist(t(vst_matrix)))
  write.csv(sample_dist, file.path(paths$qc_tables, paste0("Sample_Distance_VST_", file_tag, ".csv")), quote = FALSE)
  distance_hclust <- hclust(stats::as.dist(sample_dist))
  distance_heatmap <- pheatmap::pheatmap(
    sample_dist,
    cluster_rows = distance_hclust,
    cluster_cols = distance_hclust,
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    border_color = NA,
    main = "VST Euclidean Sample Distance",
    silent = TRUE
  )
  save_png(file.path(paths$visuals_png, paste0("Sample_Distance_Heatmap_VST_", file_tag, ".png")), 7, 6, 300, {
    grid::grid.newpage()
    grid::grid.draw(distance_heatmap$gtable)
  })
  distance_order <- data.frame(
    position = seq_along(distance_hclust$order),
    sample_id = labels(distance_hclust)[distance_hclust$order],
    stringsAsFactors = FALSE
  )
  write.csv(
    distance_order,
    file.path(paths$plot_data, paste0("Sample_Distance_Heatmap_Order_", file_tag, ".csv")),
    row.names = FALSE
  )

  vst_gene_variance <- apply(vst_matrix, 1, var)
  top_variable_genes <- annotate_gene_ids(names(vst_gene_variance), gene_map)
  top_variable_genes$vst_variance <- as.numeric(vst_gene_variance[top_variable_genes$gene_id])
  top_variable_genes <- top_variable_genes[order(top_variable_genes$vst_variance, decreasing = TRUE), ]
  write.csv(top_variable_genes, file.path(paths$qc_tables, paste0("Top_Variable_Genes_VST_", file_tag, ".csv")), row.names = FALSE)
  write.csv(
    head(top_variable_genes, min(500, nrow(top_variable_genes))),
    file.path(paths$qc_tables, paste0("Top_500_Variable_Genes_VST_", file_tag, ".csv")),
    row.names = FALSE
  )

  vst_matrix_filename <- paste0("VST_Matrix_", file_tag, ".csv")
  vst_matrix_export <- data.frame(gene_id = rownames(vst_matrix), vst_matrix, check.names = FALSE)
  write.csv(vst_matrix_export, file.path(paths$transformed_expression, vst_matrix_filename), row.names = FALSE, quote = FALSE)

  vst_annotated <- annotate_vst(vst_matrix, gene_map)
  vst_annotated_filename <- paste0("VST_Table_Annotated_", file_tag, ".csv")
  write.csv(vst_annotated, file.path(paths$transformed_expression, vst_annotated_filename), row.names = FALSE)

  group_levels <- levels(sampleTable$condition)
  vst_group_means <- data.frame(gene_id = rownames(vst_matrix), stringsAsFactors = FALSE)
  for (grp in group_levels) {
    group_samples <- rownames(sampleTable)[sampleTable$condition == grp]
    vst_group_means[[paste0("mean_", safe_name(grp))]] <- rowMeans(vst_matrix[, group_samples, drop = FALSE])
  }
  group_match_idx <- match(vst_group_means$gene_id, gene_map$gene_id)
  vst_group_means$gene <- gene_map$gene[group_match_idx]
  vst_group_means$gene_label <- ifelse(!is.na(vst_group_means$gene) & nzchar(vst_group_means$gene), vst_group_means$gene, vst_group_means$gene_id)
  vst_group_means <- vst_group_means[, c("gene_id", "gene", "gene_label", setdiff(names(vst_group_means), c("gene_id", "gene", "gene_label")))]
  group_means_filename <- paste0("Group_Means_VST_", file_tag, ".csv")
  write.csv(vst_group_means, file.path(paths$transformed_expression, group_means_filename), row.names = FALSE)

  # ============================================================================
  # SECTION 7: SIGNIFICANT GENE EXPORTS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 7: EXPORTING SIGNIFICANT GENE TABLES\n")
  cat("--------------------------------------------------------\n")

  sig_fdr <- !is.na(res_shrunken_annotated$padj) & res_shrunken_annotated$padj < padj_cutoff
  sig_fdr_lfc <- sig_fdr & !is.na(res_shrunken_annotated$log2FoldChange) & abs(res_shrunken_annotated$log2FoldChange) >= sig_lfc_cutoff
  sig_fdr_lfc_up <- sig_fdr & !is.na(res_shrunken_annotated$log2FoldChange) &
    res_shrunken_annotated$log2FoldChange >= sig_lfc_cutoff
  sig_fdr_lfc_down <- sig_fdr & !is.na(res_shrunken_annotated$log2FoldChange) &
    res_shrunken_annotated$log2FoldChange <= -sig_lfc_cutoff
  padj_na <- is.na(res_shrunken_annotated$padj)

  sig_fdr_filename <- paste0("Significant_Genes_", padj_label, "_", file_tag, ".csv")
  sig_fdr_lfc_filename <- paste0("Significant_Genes_", padj_label, "_", sig_lfc_label, "_", file_tag, ".csv")
  write.csv(res_shrunken_annotated[sig_fdr, ], file.path(paths$differential_expression, sig_fdr_filename), row.names = FALSE)
  write.csv(res_shrunken_annotated[sig_fdr_lfc, ], file.path(paths$differential_expression, sig_fdr_lfc_filename), row.names = FALSE)

  sig_gene_ids <- res_shrunken_annotated$gene_id[sig_fdr]
  sig_lfc_gene_ids <- res_shrunken_annotated$gene_id[sig_fdr_lfc]
  sig_vst_filename <- paste0("Significant_Genes_VST_", padj_label, "_Annotated_", file_tag, ".csv")
  sig_vst_lfc_filename <- paste0("Significant_Genes_VST_", padj_label, "_", sig_lfc_label, "_Annotated_", file_tag, ".csv")
  write.csv(vst_annotated[vst_annotated$gene_id %in% sig_gene_ids, ], file.path(paths$transformed_expression, sig_vst_filename), row.names = FALSE)
  write.csv(vst_annotated[vst_annotated$gene_id %in% sig_lfc_gene_ids, ], file.path(paths$transformed_expression, sig_vst_lfc_filename), row.names = FALSE)

  deg_count_summary <- data.frame(
    category = c(
      "FDR_only_total",
      "FDR_LFC_total",
      "FDR_LFC_up",
      "FDR_LFC_down",
      "non_significant",
      "padj_NA"
    ),
    count = c(
      sum(sig_fdr),
      sum(sig_fdr_lfc),
      sum(sig_fdr_lfc_up),
      sum(sig_fdr_lfc_down),
      sum(!sig_fdr & !padj_na),
      sum(padj_na)
    )
  )
  
  write.csv(
    deg_count_summary,
    file.path(paths$qc_tables, paste0("DEG_Count_Summary_", file_tag, ".csv")),
    row.names = FALSE
  )
  
  cat("[REPORT] DEG count summary:\n")
  print(deg_count_summary)
  
  filtering_summary <- data.frame(
    genes_before_low_count_filter = genes_before_filter,
    genes_after_low_count_filter = genes_after_filter,
    genes_removed_low_count_filter = genes_removed_filter,
    genes_with_non_na_pvalue = sum(!is.na(res_shrunken_annotated$pvalue)),
    genes_with_na_pvalue = sum(is.na(res_shrunken_annotated$pvalue)),
    genes_with_non_na_padj = sum(!is.na(res_shrunken_annotated$padj)),
    genes_with_na_padj = sum(is.na(res_shrunken_annotated$padj)),
    padj_cutoff = padj_cutoff,
    sig_lfc_cutoff = sig_lfc_cutoff,
    significant_FDR_count = sum(sig_fdr),
    significant_FDR_LFC_count = sum(sig_fdr_lfc),
    upregulated_FDR_LFC_count = sum(sig_fdr_lfc_up),
    downregulated_FDR_LFC_count = sum(sig_fdr_lfc_down)
  )
  write.csv(filtering_summary, file.path(paths$qc_tables, paste0("DESeq2_Filtering_Summary_", file_tag, ".csv")), row.names = FALSE)

  top_50_by_padj <- res_shrunken_annotated[!is.na(res_shrunken_annotated$padj), ]
  top_50_by_padj <- top_50_by_padj[order(top_50_by_padj$padj), ]
  top_50_up <- res_shrunken_annotated[sig_fdr_lfc_up, ]
  top_50_up <- top_50_up[order(top_50_up$padj, -top_50_up$log2FoldChange), ]
  top_50_down <- res_shrunken_annotated[sig_fdr_lfc_down, ]
  top_50_down <- top_50_down[order(top_50_down$padj, top_50_down$log2FoldChange), ]
  write.csv(head(top_50_by_padj, 50), file.path(paths$differential_expression, paste0("Top_50_DEGs_By_padj_", file_tag, ".csv")), row.names = FALSE)
  write.csv(head(top_50_up, 50), file.path(paths$differential_expression, paste0("Top_50_Upregulated_DEGs_", file_tag, ".csv")), row.names = FALSE)
  write.csv(head(top_50_down, 50), file.path(paths$differential_expression, paste0("Top_50_Downregulated_DEGs_", file_tag, ".csv")), row.names = FALSE)

  # ============================================================================
  # SECTION 8: FUNCTIONAL ANALYSIS / DAVID EXPORTS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 8: EXPORTING FUNCTIONAL ANALYSIS LISTS\n")
  cat("--------------------------------------------------------\n")

  all_filtered_ids <- rownames(dds)
  non_na_padj <- !is.na(res_shrunken_annotated$padj)
  non_na_padj_ids <- res_shrunken_annotated$gene_id[non_na_padj]

  total_deg_ids <- res_shrunken_annotated$gene_id[sig_fdr_lfc]
  up_deg_ids <- res_shrunken_annotated$gene_id[sig_fdr_lfc_up]
  down_deg_ids <- res_shrunken_annotated$gene_id[sig_fdr_lfc_down]

  david_background_all_filename <- paste0("DAVID_Background_All_Filtered_Genes_", file_tag, ".txt")
  david_background_non_na_filename <- paste0("DAVID_Background_NonNA_padj_Genes_", file_tag, ".txt")
  david_total_filename <- paste0("DAVID_Target_TOTAL_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")
  david_up_filename <- paste0("DAVID_Target_UP_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")
  david_down_filename <- paste0("DAVID_Target_DOWN_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")
  david_total_symbols_filename <- paste0("DAVID_Target_TOTAL_gene_symbols_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")
  david_up_symbols_filename <- paste0("DAVID_Target_UP_gene_symbols_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")
  david_down_symbols_filename <- paste0("DAVID_Target_DOWN_gene_symbols_", padj_label, "_", sig_lfc_label, "_", file_tag, ".txt")

  write.table(all_filtered_ids, file.path(paths$downstream, david_background_all_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(non_na_padj_ids, file.path(paths$downstream, david_background_non_na_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(total_deg_ids, file.path(paths$downstream, david_total_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(up_deg_ids, file.path(paths$downstream, david_up_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down_deg_ids, file.path(paths$downstream, david_down_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)

  total_deg_symbols <- unique(res_shrunken_annotated$gene[sig_fdr_lfc])
  up_deg_symbols <- unique(res_shrunken_annotated$gene[sig_fdr_lfc_up])
  down_deg_symbols <- unique(res_shrunken_annotated$gene[sig_fdr_lfc_down])
  total_deg_symbols <- total_deg_symbols[!is.na(total_deg_symbols) & nzchar(total_deg_symbols)]
  up_deg_symbols <- up_deg_symbols[!is.na(up_deg_symbols) & nzchar(up_deg_symbols)]
  down_deg_symbols <- down_deg_symbols[!is.na(down_deg_symbols) & nzchar(down_deg_symbols)]

  write.table(total_deg_symbols, file.path(paths$downstream, david_total_symbols_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(up_deg_symbols, file.path(paths$downstream, david_up_symbols_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down_deg_symbols, file.path(paths$downstream, david_down_symbols_filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (length(total_deg_symbols) == 0 && length(up_deg_symbols) == 0 && length(down_deg_symbols) == 0) {
    cat("[QC] No gene symbols available for DAVID gene-symbol target exports; empty files written.\n")
  }

  david_total_count <- length(total_deg_ids)
  david_up_count <- length(up_deg_ids)
  david_down_count <- length(down_deg_ids)
  cat("DAVID export counts:\n")
  cat("  ", david_background_all_filename, ": ", length(all_filtered_ids), " genes\n", sep = "")
  cat("  ", david_background_non_na_filename, ": ", length(non_na_padj_ids), " genes\n", sep = "")
  cat("  ", david_total_filename, ": ", david_total_count, " genes\n", sep = "")
  cat("  ", david_up_filename, ": ", david_up_count, " genes\n", sep = "")
  cat("  ", david_down_filename, ": ", david_down_count, " genes\n", sep = "")

  # ============================================================================
  # SECTION 9: ENRICHMENT-READY & CLUSTER-READY EXPORTS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 9: EXPORTING ENRICHMENT-READY AND CLUSTER-READY TABLES\n")
  cat("--------------------------------------------------------\n")

  rank_metric <- sign(res_shrunken_annotated$log2FoldChange) * -log10(res_shrunken_annotated$pvalue)
  rank_metric[is.infinite(rank_metric)] <- NA
  direction <- ifelse(
    is.na(res_shrunken_annotated$log2FoldChange),
    NA,
    ifelse(
      res_shrunken_annotated$log2FoldChange > 0,
      "UP",
      ifelse(res_shrunken_annotated$log2FoldChange < 0, "DOWN", "NO_CHANGE")
    )
  )
  
  direction_significant <- ifelse(
    sig_fdr_lfc_up,
    "UP",
    ifelse(sig_fdr_lfc_down, "DOWN", "NS")
  )
  enrichment_ready_cols <- list(
    gene_id = res_shrunken_annotated$gene_id,
    gene = res_shrunken_annotated$gene,
    gene_label = res_shrunken_annotated$gene_label,
    baseMean = res_shrunken_annotated$baseMean,
    log2FoldChange = res_shrunken_annotated$log2FoldChange,
    lfcSE = res_shrunken_annotated$lfcSE
  )
  if ("stat" %in% names(res_shrunken_annotated)) {
    enrichment_ready_cols$stat <- res_shrunken_annotated$stat
  }
  enrichment_ready_cols <- c(
    enrichment_ready_cols,
    list(
      pvalue = res_shrunken_annotated$pvalue,
      padj = res_shrunken_annotated$padj,
      rank_metric = rank_metric,
      direction = direction,
      direction_significant = direction_significant,
      significant_FDR = sig_fdr,
      significant_FDR_LFC = sig_fdr_lfc
    )
  )
  enrichment_ready <- data.frame(enrichment_ready_cols, stringsAsFactors = FALSE)
  enrichment_ready[[paste0("significant_", padj_label)]] <- sig_fdr
  enrichment_ready[[paste0("significant_", padj_label, "_", sig_lfc_label)]] <- sig_fdr_lfc
  enrichment_ready_filename <- paste0("Enrichment_Ready_", file_tag, ".csv")
  write.csv(enrichment_ready, file.path(paths$downstream, enrichment_ready_filename), row.names = FALSE)

  cluster_ready <- data.frame(
    gene_id = res_shrunken_annotated$gene_id,
    gene = res_shrunken_annotated$gene,
    gene_label = res_shrunken_annotated$gene_label,
    padj = res_shrunken_annotated$padj,
    stringsAsFactors = FALSE
  )
  cluster_ready[[file_tag]] <- res_shrunken_annotated$log2FoldChange
  cluster_ready <- cluster_ready[, c("gene_id", "gene", "gene_label", file_tag, "padj")]
  cluster_ready_filename <- paste0("Cluster_Ready_Log2FC_", file_tag, ".csv")
  write.csv(cluster_ready, file.path(paths$downstream, cluster_ready_filename), row.names = FALSE)

  # ============================================================================
  # SECTION 10: PCA & VOLCANO PLOTS
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 10: GENERATING PCA AND VOLCANO PLOTS\n")
  cat("--------------------------------------------------------\n")

  p_pca <- PCAtools::pca(vst_matrix, metadata = sampleTable, removeVar = 0.1)
  pc1_var <- p_pca$variance[1]
  pc2_var <- p_pca$variance[2]
  cat("PCA variance: PC1 = ", pc1_var, "%; PC2 = ", pc2_var, "%\n", sep = "")

  pca_scores <- as.data.frame(p_pca$rotated, check.names = FALSE)
  pca_scores$sample_id <- rownames(pca_scores)
  pca_scores$condition <- as.character(sampleTable[pca_scores$sample_id, "condition"])
  pca_scores <- pca_scores[, c("sample_id", "condition", setdiff(names(pca_scores), c("sample_id", "condition")))]
  write.csv(
    pca_scores,
    file.path(paths$plot_data, paste0("PCA_Sample_Scores_", file_tag, ".csv")),
    row.names = FALSE
  )

  pca_loadings <- as.data.frame(p_pca$loadings, check.names = FALSE)
  pca_loadings$gene_id <- rownames(pca_loadings)
  pca_loadings <- pca_loadings[, c("gene_id", setdiff(names(pca_loadings), "gene_id"))]
  write.csv(
    pca_loadings,
    file.path(paths$plot_data, paste0("PCA_Gene_Loadings_", file_tag, ".csv")),
    row.names = FALSE
  )

  pca_variance <- data.frame(
    component = names(p_pca$variance),
    variance_percent = as.numeric(p_pca$variance),
    stringsAsFactors = FALSE
  )
  write.csv(
    pca_variance,
    file.path(paths$plot_data, paste0("PCA_Variance_Explained_", file_tag, ".csv")),
    row.names = FALSE
  )

  pca_plot <- PCAtools::biplot(
    p_pca,
    colby = "condition",
    colkey = condition_colors,
    legendPosition = "right",
    encircle = TRUE,
    encircleFill = TRUE,
    title = paste("PCA:", project_tag)
  )
  pca_plot_filename <- paste0("PCA_", file_tag, ".png")
  save_png(file.path(paths$visuals_png, pca_plot_filename), 6, 5, 300, {
    print(pca_plot)
  })

  volcano_candidates <- res_shrunken_annotated[
    !is.na(res_shrunken_annotated$padj) &
      !is.na(res_shrunken_annotated$log2FoldChange) &
      res_shrunken_annotated$padj < padj_cutoff &
      abs(res_shrunken_annotated$log2FoldChange) >= volcano_lfc_cutoff,
    ,
    drop = FALSE
  ]
  volcano_candidates <- volcano_candidates[order(volcano_candidates$padj), , drop = FALSE]
  select_labels <- head(volcano_candidates$gene_label, top_n_volcano_labels)
  riba2_idx <- (!is.na(res_shrunken_annotated$gene) &
                  tolower(res_shrunken_annotated$gene) == "riba2") |
    res_shrunken_annotated$gene_id == "ACP3PV_RS07545"
  riba2_labels <- res_shrunken_annotated$gene_label[riba2_idx]
  select_labels <- unique(c(select_labels, riba2_labels))

  volcano_plot_data <- res_shrunken_annotated
  volcano_plot_data$passes_padj <- !is.na(volcano_plot_data$padj) &
    volcano_plot_data$padj < padj_cutoff
  volcano_plot_data$passes_abs_lfc <- !is.na(volcano_plot_data$log2FoldChange) &
    abs(volcano_plot_data$log2FoldChange) >= volcano_lfc_cutoff
  volcano_plot_data$volcano_category <- ifelse(
    volcano_plot_data$passes_padj & volcano_plot_data$passes_abs_lfc,
    "adjusted_p_and_abs_log2FC",
    ifelse(
      volcano_plot_data$passes_padj,
      "adjusted_p_only",
      ifelse(volcano_plot_data$passes_abs_lfc, "abs_log2FC_only", "not_significant")
    )
  )
  volcano_plot_data$selected_for_label <- volcano_plot_data$gene_label %in% select_labels
  volcano_plot_data$negative_log10_padj <- ifelse(
    is.na(volcano_plot_data$padj),
    NA_real_,
    -log10(volcano_plot_data$padj)
  )
  write.csv(
    volcano_plot_data,
    file.path(paths$plot_data, paste0("Volcano_Plot_Data_", file_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    volcano_plot_data[volcano_plot_data$selected_for_label, c(
      "gene_id", "gene", "gene_label", "log2FoldChange", "padj", "volcano_category"
    )],
    file.path(paths$plot_data, paste0("Volcano_Selected_Labels_", file_tag, ".csv")),
    row.names = FALSE
  )

  v_plot <- EnhancedVolcano(
    res_shrunken_annotated,
    lab = res_shrunken_annotated$gene_label,
    selectLab = select_labels,
    x = "log2FoldChange",
    y = "padj",
    title = project_tag,
    subtitle = paste(test_name, "vs", ctrl_name),
    xlab = bquote(~Log[2]~"fold change"),
    ylab = bquote(~-Log[10]~"adjusted p-value"),
    pCutoff = padj_cutoff,
    FCcutoff = volcano_lfc_cutoff,
    legendLabels = c(
      "NS",
      "|log2FC|",
      "adjusted p-value",
      "adjusted p-value and |log2FC|"
    ),
    legendPosition = "bottom",
    drawConnectors = TRUE
  ) + ggplot2::theme_bw()

  volcano_plot_filename <- paste0("Volcano_", file_tag, ".png")
  save_png(file.path(paths$visuals_png, volcano_plot_filename), 10, 8.5, 600, {
    print(v_plot)
  })

  plot_specifications <- data.frame(
    figure = c(
      "MA_Unshrunken",
      "MA_apeglm_Shrunken",
      "PValue_Histogram_Unshrunken",
      "Dispersion_Plot",
      "Sample_Correlation_Heatmap_VST",
      "Sample_Distance_Heatmap_VST",
      "PCA",
      "Volcano"
    ),
    width_inches = c(6, 6, 6, 6, 7, 7, 6, 10),
    height_inches = c(6, 6, 5, 5, 6, 6, 5, 8.5),
    png_dpi = c(300, 300, 300, 300, 300, 300, 300, 600),
    formats = "PNG;PDF",
    padj_cutoff = padj_cutoff,
    sig_lfc_cutoff = sig_lfc_cutoff,
    volcano_lfc_cutoff = volcano_lfc_cutoff,
    stringsAsFactors = FALSE
  )
  write.csv(
    plot_specifications,
    file.path(paths$plot_data, "Plot_Specifications.csv"),
    row.names = FALSE
  )

  # ============================================================================
  # SECTION 11: SESSION INFO & FINAL ARCHIVE
  # ============================================================================
  cat("\n--------------------------------------------------------\n")
  cat("SECTION 11: SAVING SESSION INFO AND NAMED OBJECT ARCHIVES\n")
  cat("--------------------------------------------------------\n")

  session_info_filename <- paste0("Session_Info_", run_id, ".txt")
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(paths$logs, session_info_filename)
  )

  de_results_alpha <- tryCatch(
    as.character(S4Vectors::metadata(res_unshrunken)$alpha),
    error = function(e) NA_character_
  )
  dispersion_fit_type <- tryCatch(
    as.character(attr(DESeq2::dispersionFunction(dds), "fitType")),
    error = function(e) NA_character_
  )
  if (length(dispersion_fit_type) == 0) {
    dispersion_fit_type <- NA_character_
  }

  run_manifest <- data.frame(
    field = c(
      "script_name",
      "script_version",
      "script_source_path",
      "script_archived_path",
      "script_checksum_algorithm",
      "script_checksum",
      "reference_tag",
      "project_tag",
      "ctrl_name",
      "test_name",
      "contrast_tag",
      "file_tag",
      "run_id",
      "run_stamp",
      "start_time_iso8601",
      "timezone",
      "output_root",
      "input_counts_original",
      "input_counts_archived",
      "input_counts_checksum",
      "input_gtf_original",
      "input_gtf_archived",
      "input_gtf_checksum",
      "count_file_copied",
      "gtf_file_copied",
      "ctrl_indices",
      "test_indices",
      "padj_cutoff",
      "sig_lfc_cutoff",
      "volcano_lfc_cutoff",
      "top_n_volcano_labels",
      "low_count_filter",
      "design_formula",
      "shrinkage_type",
      "DESeq2_results_alpha",
      "dispersion_fit_type",
      "VST_blind",
      "PCA_removeVar",
      "visual_formats",
      "folder_layout_version",
      "genes_before_filter",
      "genes_after_filter",
      "gene_id_match_rate_percent",
      "R_version",
      "platform",
      "DESeq2_version",
      "PCAtools_version",
      "EnhancedVolcano_version",
      "apeglm_version",
      "rtracklayer_version",
      "ggplot2_version",
      "pheatmap_version"
    ),
    value = c(
      if (is.na(script_path)) script_basename else basename(script_path),
      script_version,
      if (is.na(script_path)) "" else script_path,
      if (is.na(archived_script_path)) "" else .deseq2_relative_path(archived_script_path, output_root),
      .deseq2_checksum_algorithm(),
      script_checksum,
      reference_tag,
      project_tag,
      ctrl_name,
      test_name,
      contrast_tag,
      file_tag,
      run_id,
      run_stamp,
      format(start_time, "%Y-%m-%dT%H:%M:%S%z"),
      Sys.timezone(),
      output_root,
      normalizePath(input_counts, winslash = "/", mustWork = TRUE),
      .deseq2_relative_path(copied_counts_path, output_root),
      input_provenance$checksum[input_provenance$role == "count_table"],
      normalizePath(input_gtf, winslash = "/", mustWork = TRUE),
      .deseq2_relative_path(copied_gtf_path, output_root),
      input_provenance$checksum[input_provenance$role == "gtf_annotation"],
      as.character(copied_counts),
      as.character(copied_gtf),
      paste(ctrl_indices, collapse = ","),
      paste(test_indices, collapse = ","),
      as.character(padj_cutoff),
      as.character(sig_lfc_cutoff),
      as.character(volcano_lfc_cutoff),
      as.character(top_n_volcano_labels),
      "rowSums(counts(dds)) >= 10",
      "~ condition",
      "apeglm",
      de_results_alpha,
      dispersion_fit_type,
      "FALSE",
      "0.1",
      "PNG;PDF",
      "v1_numbered_reproducible_bundle",
      as.character(genes_before_filter),
      as.character(genes_after_filter),
      as.character(round(id_match_rate, 1)),
      R.version.string,
      R.version$platform,
      as.character(packageVersion("DESeq2")),
      as.character(packageVersion("PCAtools")),
      as.character(packageVersion("EnhancedVolcano")),
      as.character(packageVersion("apeglm")),
      as.character(packageVersion("rtracklayer")),
      as.character(packageVersion("ggplot2")),
      as.character(packageVersion("pheatmap"))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    run_manifest,
    file.path(paths$provenance, paste0("Run_Manifest_", file_tag, ".csv")),
    row.names = FALSE
  )

  saveRDS(dds, file.path(paths$objects, paste0("DESeqDataSet_", file_tag, ".rds")))
  saveRDS(vsd, file.path(paths$objects, paste0("VST_Object_", file_tag, ".rds")))
  saveRDS(res_unshrunken, file.path(paths$objects, paste0("DESeq2_Unshrunken_Results_", file_tag, ".rds")))
  saveRDS(resLFC, file.path(paths$objects, paste0("DESeq2_apeglm_Results_", file_tag, ".rds")))
  saveRDS(p_pca, file.path(paths$objects, paste0("PCA_Object_", file_tag, ".rds")))
  saveRDS(pca_plot, file.path(paths$objects, paste0("PCA_Plot_", file_tag, ".rds")))
  saveRDS(v_plot, file.path(paths$objects, paste0("Volcano_Plot_", file_tag, ".rds")))
  saveRDS(
    correlation_heatmap,
    file.path(paths$objects, paste0("Correlation_Heatmap_", file_tag, ".rds"))
  )
  saveRDS(
    distance_heatmap,
    file.path(paths$objects, paste0("Distance_Heatmap_", file_tag, ".rds"))
  )
  saveRDS(
    utils::sessionInfo(),
    file.path(paths$objects, paste0("Session_Info_", run_id, ".rds"))
  )

  named_archive_filename <- paste0("Named_Analysis_Objects_", run_id, ".RData")
  save(
    dds,
    vsd,
    res_unshrunken,
    resLFC,
    p_pca,
    pca_plot,
    v_plot,
    correlation_heatmap,
    distance_heatmap,
    sampleTable,
    sample_mapping,
    gene_map,
    run_manifest,
    file = file.path(paths$objects, named_archive_filename),
    envir = environment()
  )

  cat("\n========================================================\n")
  cat("PIPELINE COMPLETE\n")
  cat("========================================================\n")
  cat("output_root: ", output_root, "\n", sep = "")
  cat("unshrunken result filename: ", unshrunken_result_filename, "\n", sep = "")
  cat("shrunken result filename: ", shrunken_result_filename, "\n", sep = "")
  cat("VST annotated filename: ", vst_annotated_filename, "\n", sep = "")
  cat("execution log filename: ", basename(log_filename), "\n", sep = "")
  cat("DAVID total/up/down counts: ", david_total_count, "/", david_up_count, "/", david_down_count, "\n", sep = "")
  cat("========================================================\n")

  invisible(output_root)
}

.deseq2_run <- function() {
  context <- new.env(parent = emptyenv())
  warnings_log <- list()
  run_error <- NULL

  withCallingHandlers(
    tryCatch(
      main(context),
      error = function(error_condition) {
        run_error <<- error_condition
        invisible(NULL)
      }
    ),
    warning = function(warning_condition) {
      warning_call <- conditionCall(warning_condition)
      warnings_log[[length(warnings_log) + 1]] <<- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        message = conditionMessage(warning_condition),
        call = if (is.null(warning_call)) "" else paste(deparse(warning_call), collapse = " "),
        stringsAsFactors = FALSE
      )
    }
  )

  .deseq2_finalize_run(context, warnings_log, run_error)

  if (!is.null(run_error)) {
    stop(conditionMessage(run_error), call. = FALSE)
  }

  invisible(context$output_root)
}

.deseq2_run()
