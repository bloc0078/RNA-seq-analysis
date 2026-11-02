#!/usr/bin/env Rscript

# Flavoprotein enrichment analysis for DESeq2 result files
#
# Run this script from a directory containing:
#   1. One or more DESeq2 CSV files beginning with "DESeq2"
#   2. One UniProt flavoprotein CSV file named "UniProt_flavoproteins.csv"
#
# Example:
#   Rscript flavoprotein_enrichment.R


# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

config <- list(
  input_directory = ".",
  deseq2_pattern = "^DESeq2.*\\.csv$",
  flavoprotein_pattern = "^UniProt_flavoproteins\\.csv$",
  output_directory = "flavoprotein_enrichment_output",
  padj_cutoff = 0.05,
  log2fc_cutoff = 1,
  enrichment_fdr_cutoff = 0.05,
  minimum_class_size_for_plot = 3,
  maximum_genes_per_ranked_plot = 40,
  save_pdf_plots = TRUE
)


# -----------------------------------------------------------------------------
# 2. General helpers and validation
# -----------------------------------------------------------------------------

run_started <- Sys.time()
run_log <- character()

log_message <- function(...) {
  text <- paste0(...)
  timestamped <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    text
  )
  run_log <<- c(run_log, timestamped)
  message(text)
}

stop_with_log <- function(...) {
  text <- paste0(...)
  log_message("ERROR: ", text)
  stop(text, call. = FALSE)
}

collapse_values <- function(values) {
  values <- as.character(values)
  values <- unique(values[
    !is.na(values) & trimws(values) != ""
  ])
  paste(values, collapse = "; ")
}

validate_single_number <- function(
  value,
  setting_name,
  lower,
  upper = Inf,
  lower_inclusive = TRUE,
  upper_inclusive = TRUE
) {
  valid <- is.numeric(value) &&
    length(value) == 1 &&
    !is.na(value) &&
    is.finite(value)

  if (valid) {
    lower_ok <- if (lower_inclusive) value >= lower else value > lower
    upper_ok <- if (upper_inclusive) value <= upper else value < upper
    valid <- lower_ok && upper_ok
  }

  if (!valid) {
    stop_with_log("Invalid setting: ", setting_name, " = ", value)
  }
}

validate_single_number(
  config$padj_cutoff,
  "padj_cutoff",
  lower = 0,
  upper = 1,
  lower_inclusive = FALSE
)

validate_single_number(
  config$log2fc_cutoff,
  "log2fc_cutoff",
  lower = 0
)

validate_single_number(
  config$enrichment_fdr_cutoff,
  "enrichment_fdr_cutoff",
  lower = 0,
  upper = 1,
  lower_inclusive = FALSE
)

validate_single_number(
  config$minimum_class_size_for_plot,
  "minimum_class_size_for_plot",
  lower = 1
)

validate_single_number(
  config$maximum_genes_per_ranked_plot,
  "maximum_genes_per_ranked_plot",
  lower = 1
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop_with_log(
    "The ggplot2 package is required. Install it with ",
    "install.packages('ggplot2') and run the script again."
  )
}

input_directory <- normalizePath(
  config$input_directory,
  mustWork = TRUE
)

output_directory <- file.path(
  input_directory,
  config$output_directory
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

table_directory <- file.path(output_directory, "tables")
plot_directory <- file.path(output_directory, "plots")
classified_directory <- file.path(output_directory, "classified_deseq2")

for (directory in c(
  table_directory,
  plot_directory,
  classified_directory
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

log_message("Analysis started in: ", input_directory)


# -----------------------------------------------------------------------------
# 3. Locate and import input files
# -----------------------------------------------------------------------------

deseq2_files <- list.files(
  path = input_directory,
  pattern = config$deseq2_pattern,
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = FALSE
)

if (length(deseq2_files) == 0) {
  stop_with_log(
    "No DESeq2 CSV files matched '",
    config$deseq2_pattern,
    "' in ",
    input_directory
  )
}

flavoprotein_files <- list.files(
  path = input_directory,
  pattern = config$flavoprotein_pattern,
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = FALSE
)

if (length(flavoprotein_files) != 1) {
  stop_with_log(
    "Expected exactly one flavoprotein CSV matching '",
    config$flavoprotein_pattern,
    "', but found ",
    length(flavoprotein_files),
    "."
  )
}

flavoprotein_file <- flavoprotein_files[[1]]

log_message(
  "Found ",
  length(deseq2_files),
  " DESeq2 file(s): ",
  paste(basename(deseq2_files), collapse = ", ")
)
log_message("Using flavoprotein file: ", basename(flavoprotein_file))

read_csv_preserving_names <- function(file_path) {
  tryCatch(
    read.csv(
      file_path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(error) {
      stop_with_log(
        "Could not read ",
        basename(file_path),
        ": ",
        conditionMessage(error)
      )
    }
  )
}

coerce_numeric_strict <- function(values, column_name, file_name) {
  if (is.numeric(values)) {
    return(values)
  }

  converted <- suppressWarnings(as.numeric(values))
  bad_value <- !is.na(values) &
    trimws(as.character(values)) != "" &
    is.na(converted)

  if (any(bad_value)) {
    stop_with_log(
      "Column '",
      column_name,
      "' contains non-numeric values in ",
      file_name,
      "."
    )
  }

  converted
}

required_deseq2_columns <- c(
  "gene_id",
  "log2FoldChange",
  "padj"
)

read_deseq2_file <- function(file_path) {
  results <- read_csv_preserving_names(file_path)
  missing_columns <- setdiff(required_deseq2_columns, names(results))

  if (length(missing_columns) > 0) {
    stop_with_log(
      basename(file_path),
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  results$gene_id <- trimws(as.character(results$gene_id))
  results$log2FoldChange <- coerce_numeric_strict(
    results$log2FoldChange,
    "log2FoldChange",
    basename(file_path)
  )
  results$padj <- coerce_numeric_strict(
    results$padj,
    "padj",
    basename(file_path)
  )

  if ("pvalue" %in% names(results)) {
    results$pvalue <- coerce_numeric_strict(
      results$pvalue,
      "pvalue",
      basename(file_path)
    )
  }

  nonempty_ids <- results$gene_id[
    !is.na(results$gene_id) & results$gene_id != ""
  ]
  duplicate_ids <- unique(nonempty_ids[duplicated(nonempty_ids)])

  if (length(duplicate_ids) > 0) {
    stop_with_log(
      basename(file_path),
      " contains duplicate gene_id values: ",
      paste(head(duplicate_ids, 10), collapse = ", "),
      if (length(duplicate_ids) > 10) " ..." else ""
    )
  }

  results$source_file <- basename(file_path)
  results
}

dataset_names <- tools::file_path_sans_ext(basename(deseq2_files))

if (anyDuplicated(dataset_names)) {
  stop_with_log("DESeq2 input filenames do not produce unique dataset names.")
}

deseq2_results <- setNames(
  lapply(deseq2_files, read_deseq2_file),
  dataset_names
)

make_dataset_label <- function(dataset_name) {
  label <- sub(
    "^DESeq2_apeglm_Shrunken_Results_",
    "",
    dataset_name,
    ignore.case = TRUE
  )
  label <- sub("^DESeq2_", "", label, ignore.case = TRUE)
  label <- gsub("_+", " ", label)
  trimws(label)
}

dataset_labels <- setNames(
  vapply(dataset_names, make_dataset_label, character(1)),
  dataset_names
)

flavoproteins <- read_csv_preserving_names(flavoprotein_file)
required_flavoprotein_columns <- c(
  "Rv_ID",
  "Gene_name",
  "flavin_class"
)

missing_flavoprotein_columns <- setdiff(
  required_flavoprotein_columns,
  names(flavoproteins)
)

if (length(missing_flavoprotein_columns) > 0) {
  stop_with_log(
    basename(flavoprotein_file),
    " is missing required columns: ",
    paste(missing_flavoprotein_columns, collapse = ", ")
  )
}

flavoproteins$Rv_ID <- trimws(as.character(flavoproteins$Rv_ID))
flavoproteins$Gene_name <- trimws(as.character(flavoproteins$Gene_name))
flavoproteins$flavin_class <- trimws(
  as.character(flavoproteins$flavin_class)
)

empty_flavoprotein_id <- is.na(flavoproteins$Rv_ID) |
  flavoproteins$Rv_ID == ""

if (any(empty_flavoprotein_id)) {
  stop_with_log(
    "The flavoprotein file contains ",
    sum(empty_flavoprotein_id),
    " row(s) without an Rv_ID."
  )
}

duplicate_flavoprotein_ids <- unique(
  flavoproteins$Rv_ID[duplicated(flavoproteins$Rv_ID)]
)

if (length(duplicate_flavoprotein_ids) > 0) {
  stop_with_log(
    "The flavoprotein file contains duplicate Rv_ID values: ",
    paste(head(duplicate_flavoprotein_ids, 10), collapse = ", "),
    if (length(duplicate_flavoprotein_ids) > 10) " ..." else ""
  )
}

flavoproteins$flavin_class[
  is.na(flavoproteins$flavin_class) |
    flavoproteins$flavin_class == ""
] <- "Unclassified"

flavoprotein_ids <- flavoproteins$Rv_ID
flavin_classes <- sort(unique(flavoproteins$flavin_class))

log_message(
  "Loaded ",
  nrow(flavoproteins),
  " uniquely identified flavoproteins across ",
  length(flavin_classes),
  " flavin classes."
)


# -----------------------------------------------------------------------------
# 4. Classify differential expression
# -----------------------------------------------------------------------------

expression_class_levels <- c(
  "upregulated",
  "downregulated",
  "not_significant",
  "not_testable",
  "missing_gene_id"
)

classify_deseq2_results <- function(results) {
  missing_gene_id <- is.na(results$gene_id) | results$gene_id == ""
  missing_statistics <- is.na(results$padj) |
    is.na(results$log2FoldChange)

  passes_padj <- !is.na(results$padj) &
    results$padj <= config$padj_cutoff
  passes_log2fc <- !is.na(results$log2FoldChange) &
    abs(results$log2FoldChange) >= config$log2fc_cutoff

  results$is_significant <- !missing_gene_id &
    !missing_statistics &
    passes_padj &
    passes_log2fc

  expression_class <- rep("not_significant", nrow(results))
  expression_class[missing_statistics] <- "not_testable"
  expression_class[missing_gene_id] <- "missing_gene_id"
  expression_class[
    results$is_significant &
      results$log2FoldChange >= config$log2fc_cutoff
  ] <- "upregulated"
  expression_class[
    results$is_significant &
      results$log2FoldChange <= -config$log2fc_cutoff
  ] <- "downregulated"

  results$expression_class <- factor(
    expression_class,
    levels = expression_class_levels
  )

  flavin_match <- match(results$gene_id, flavoproteins$Rv_ID)
  results$is_flavoprotein <- !is.na(flavin_match)
  results$flavin_class <- flavoproteins$flavin_class[flavin_match]
  results
}

classified_results <- lapply(
  deseq2_results,
  classify_deseq2_results
)

classification_summary <- do.call(
  rbind,
  lapply(names(classified_results), function(dataset_name) {
    counts <- table(
      factor(
        classified_results[[dataset_name]]$expression_class,
        levels = expression_class_levels
      )
    )

    data.frame(
      dataset = dataset_name,
      dataset_label = dataset_labels[[dataset_name]],
      expression_class = names(counts),
      gene_count = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
)
row.names(classification_summary) <- NULL


# -----------------------------------------------------------------------------
# 5. Input matching diagnostics
# -----------------------------------------------------------------------------

matching_summary <- do.call(
  rbind,
  lapply(names(classified_results), function(dataset_name) {
    results <- classified_results[[dataset_name]]
    valid_id <- !is.na(results$gene_id) & results$gene_id != ""
    result_ids <- unique(results$gene_id[valid_id])
    matched_ids <- intersect(flavoprotein_ids, result_ids)
    testable_ids <- unique(results$gene_id[
      valid_id &
        !is.na(results$padj) &
        !is.na(results$log2FoldChange)
    ])

    data.frame(
      dataset = dataset_name,
      dataset_label = dataset_labels[[dataset_name]],
      deseq2_unique_gene_ids = length(result_ids),
      flavoproteins_in_deseq2 = length(matched_ids),
      flavoproteins_absent_from_deseq2 = length(
        setdiff(flavoprotein_ids, result_ids)
      ),
      testable_deseq2_genes = length(testable_ids),
      testable_flavoproteins = length(
        intersect(flavoprotein_ids, testable_ids)
      ),
      stringsAsFactors = FALSE
    )
  })
)

unmatched_flavoproteins <- do.call(
  rbind,
  lapply(names(classified_results), function(dataset_name) {
    result_ids <- unique(classified_results[[dataset_name]]$gene_id)
    unmatched_ids <- setdiff(flavoprotein_ids, result_ids)

    if (length(unmatched_ids) == 0) {
      return(NULL)
    }

    rows <- match(unmatched_ids, flavoproteins$Rv_ID)
    data.frame(
      dataset = dataset_name,
      dataset_label = dataset_labels[[dataset_name]],
      flavoproteins[rows, , drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
)

if (is.null(unmatched_flavoproteins)) {
  unmatched_flavoproteins <- data.frame(
    dataset = character(),
    dataset_label = character(),
    Rv_ID = character(),
    stringsAsFactors = FALSE
  )
}

for (row_index in seq_len(nrow(matching_summary))) {
  log_message(
    matching_summary$dataset_label[[row_index]],
    ": matched ",
    matching_summary$flavoproteins_in_deseq2[[row_index]],
    " flavoproteins; ",
    matching_summary$flavoproteins_absent_from_deseq2[[row_index]],
    " absent from DESeq2; ",
    matching_summary$testable_flavoproteins[[row_index]],
    " testable."
  )
}


# -----------------------------------------------------------------------------
# 6. Fisher enrichment helpers
# -----------------------------------------------------------------------------

directions_to_test <- c("upregulated", "downregulated")

run_fisher_test <- function(contingency_table) {
  test_is_possible <- all(rowSums(contingency_table) > 0) &&
    all(colSums(contingency_table) > 0)

  if (!test_is_possible) {
    return(list(
      status = "insufficient_categories",
      p_value = NA_real_,
      odds_ratio = NA_real_,
      odds_ratio_ci_low = NA_real_,
      odds_ratio_ci_high = NA_real_
    ))
  }

  result <- stats::fisher.test(
    contingency_table,
    alternative = "greater",
    conf.level = 0.95
  )

  list(
    status = "completed",
    p_value = result$p.value,
    odds_ratio = unname(result$estimate),
    odds_ratio_ci_low = result$conf.int[[1]],
    odds_ratio_ci_high = result$conf.int[[2]]
  )
}

get_testable_and_selected_ids <- function(results, direction) {
  valid_universe_row <- !is.na(results$gene_id) &
    results$gene_id != "" &
    !is.na(results$padj) &
    !is.na(results$log2FoldChange)

  selected_row <- valid_universe_row &
    !is.na(results$expression_class) &
    results$expression_class == direction

  list(
    universe_ids = unique(results$gene_id[valid_universe_row]),
    selected_ids = unique(results$gene_id[selected_row])
  )
}

make_long_contingency_table <- function(
  contingency_table,
  dataset_name,
  direction,
  test_type,
  flavin_class = NA_character_,
  reference_population = NA_character_
) {
  expanded <- as.data.frame(as.table(contingency_table))
  names(expanded) <- c(
    "selection_status",
    "feature_status",
    "gene_count"
  )

  data.frame(
    dataset = dataset_name,
    dataset_label = dataset_labels[[dataset_name]],
    direction = direction,
    test_type = test_type,
    flavin_class = flavin_class,
    reference_population = reference_population,
    expanded,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# 7. Overall flavoprotein enrichment
# -----------------------------------------------------------------------------

run_overall_enrichment <- function(results, dataset_name, direction) {
  id_sets <- get_testable_and_selected_ids(results, direction)
  universe_ids <- id_sets$universe_ids
  selected_ids <- id_sets$selected_ids
  background_ids <- setdiff(universe_ids, selected_ids)

  universe_flavoproteins <- intersect(universe_ids, flavoprotein_ids)
  selected_flavoproteins <- intersect(
    selected_ids,
    universe_flavoproteins
  )
  background_flavoproteins <- intersect(
    background_ids,
    universe_flavoproteins
  )

  selected_nonflavoproteins <- setdiff(
    selected_ids,
    universe_flavoproteins
  )
  background_nonflavoproteins <- setdiff(
    background_ids,
    universe_flavoproteins
  )

  contingency_table <- matrix(
    c(
      length(selected_flavoproteins),
      length(selected_nonflavoproteins),
      length(background_flavoproteins),
      length(background_nonflavoproteins)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      selection_status = c("selected", "not_selected"),
      feature_status = c("flavoprotein", "not_flavoprotein")
    )
  )

  fisher_result <- run_fisher_test(contingency_table)

  universe_size <- length(universe_ids)
  selected_size <- length(selected_ids)
  universe_flavoprotein_count <- length(universe_flavoproteins)
  selected_flavoprotein_count <- length(selected_flavoproteins)

  expected_flavoproteins <- if (universe_size > 0) {
    selected_size * universe_flavoprotein_count / universe_size
  } else {
    NA_real_
  }

  fold_enrichment <- if (
    !is.na(expected_flavoproteins) && expected_flavoproteins > 0
  ) {
    selected_flavoprotein_count / expected_flavoproteins
  } else {
    NA_real_
  }

  summary <- data.frame(
    dataset = dataset_name,
    dataset_label = dataset_labels[[dataset_name]],
    direction = direction,
    test_status = fisher_result$status,
    universe_genes = universe_size,
    universe_flavoproteins = universe_flavoprotein_count,
    selected_genes = selected_size,
    selected_flavoproteins = selected_flavoprotein_count,
    expected_flavoproteins = expected_flavoproteins,
    fold_enrichment = fold_enrichment,
    odds_ratio = fisher_result$odds_ratio,
    odds_ratio_ci_low = fisher_result$odds_ratio_ci_low,
    odds_ratio_ci_high = fisher_result$odds_ratio_ci_high,
    p_value = fisher_result$p_value,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    contingency_table = contingency_table,
    contingency_long = make_long_contingency_table(
      contingency_table,
      dataset_name,
      direction,
      test_type = "overall_flavoprotein"
    ),
    flavoprotein_hits = selected_flavoproteins
  )
}

overall_enrichment_tests <- list()

for (dataset_name in names(classified_results)) {
  for (direction in directions_to_test) {
    test_key <- paste(dataset_name, direction, sep = "__")
    overall_enrichment_tests[[test_key]] <- run_overall_enrichment(
      classified_results[[dataset_name]],
      dataset_name,
      direction
    )
  }
}

overall_enrichment_summary <- do.call(
  rbind,
  lapply(overall_enrichment_tests, function(test) test$summary)
)
row.names(overall_enrichment_summary) <- NULL

overall_enrichment_summary$enrichment_FDR <- NA_real_
valid_overall_pvalue <- !is.na(overall_enrichment_summary$p_value)
overall_enrichment_summary$enrichment_FDR[valid_overall_pvalue] <- p.adjust(
  overall_enrichment_summary$p_value[valid_overall_pvalue],
  method = "BH"
)

get_hit_annotation <- function(hit_ids, column_name) {
  hit_rows <- match(hit_ids, flavoproteins$Rv_ID)
  collapse_values(flavoproteins[[column_name]][hit_rows])
}

overall_enrichment_summary$flavoprotein_hit_ids <- vapply(
  overall_enrichment_tests,
  function(test) collapse_values(test$flavoprotein_hits),
  character(1)
)
overall_enrichment_summary$flavoprotein_hit_genes <- vapply(
  overall_enrichment_tests,
  function(test) {
    get_hit_annotation(test$flavoprotein_hits, "Gene_name")
  },
  character(1)
)

overall_contingency_tables <- do.call(
  rbind,
  lapply(
    overall_enrichment_tests,
    function(test) test$contingency_long
  )
)
row.names(overall_contingency_tables) <- NULL


# -----------------------------------------------------------------------------
# 8. Enrichment by flavin class
# -----------------------------------------------------------------------------

# Two reference populations are evaluated:
#   all_testable_genes  - Is this class enriched among all DE genes?
#   flavoproteins_only  - Is this class enriched within DE flavoproteins?

class_reference_populations <- c(
  "all_testable_genes",
  "flavoproteins_only"
)

run_class_enrichment <- function(
  results,
  dataset_name,
  direction,
  class_name,
  class_ids,
  reference_population
) {
  id_sets <- get_testable_and_selected_ids(results, direction)
  all_universe_ids <- id_sets$universe_ids
  all_selected_ids <- id_sets$selected_ids

  if (reference_population == "all_testable_genes") {
    universe_ids <- all_universe_ids
    selected_ids <- all_selected_ids
  } else if (reference_population == "flavoproteins_only") {
    universe_ids <- intersect(all_universe_ids, flavoprotein_ids)
    selected_ids <- intersect(all_selected_ids, flavoprotein_ids)
  } else {
    stop_with_log(
      "Unknown class reference population: ",
      reference_population
    )
  }

  background_ids <- setdiff(universe_ids, selected_ids)
  universe_class_ids <- intersect(universe_ids, class_ids)
  selected_class_ids <- intersect(selected_ids, universe_class_ids)
  background_class_ids <- intersect(background_ids, universe_class_ids)
  selected_nonclass_ids <- setdiff(selected_ids, universe_class_ids)
  background_nonclass_ids <- setdiff(background_ids, universe_class_ids)

  contingency_table <- matrix(
    c(
      length(selected_class_ids),
      length(selected_nonclass_ids),
      length(background_class_ids),
      length(background_nonclass_ids)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      selection_status = c("selected", "not_selected"),
      feature_status = c("class_member", "not_class_member")
    )
  )

  fisher_result <- run_fisher_test(contingency_table)
  universe_size <- length(universe_ids)
  selected_size <- length(selected_ids)
  universe_class_size <- length(universe_class_ids)
  selected_class_size <- length(selected_class_ids)

  expected_class_genes <- if (universe_size > 0) {
    selected_size * universe_class_size / universe_size
  } else {
    NA_real_
  }

  fold_enrichment <- if (
    !is.na(expected_class_genes) && expected_class_genes > 0
  ) {
    selected_class_size / expected_class_genes
  } else {
    NA_real_
  }

  summary <- data.frame(
    dataset = dataset_name,
    dataset_label = dataset_labels[[dataset_name]],
    direction = direction,
    reference_population = reference_population,
    flavin_class = class_name,
    test_status = fisher_result$status,
    universe_genes = universe_size,
    universe_class_genes = universe_class_size,
    selected_genes = selected_size,
    selected_class_genes = selected_class_size,
    expected_class_genes = expected_class_genes,
    fold_enrichment = fold_enrichment,
    odds_ratio = fisher_result$odds_ratio,
    odds_ratio_ci_low = fisher_result$odds_ratio_ci_low,
    odds_ratio_ci_high = fisher_result$odds_ratio_ci_high,
    p_value = fisher_result$p_value,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    contingency_table = contingency_table,
    contingency_long = make_long_contingency_table(
      contingency_table,
      dataset_name,
      direction,
      test_type = "flavin_class",
      flavin_class = class_name,
      reference_population = reference_population
    ),
    class_hits = selected_class_ids
  )
}

class_enrichment_tests <- list()

for (dataset_name in names(classified_results)) {
  for (direction in directions_to_test) {
    for (reference_population in class_reference_populations) {
      for (class_name in flavin_classes) {
        class_ids <- flavoproteins$Rv_ID[
          flavoproteins$flavin_class == class_name
        ]
        test_key <- paste(
          dataset_name,
          direction,
          reference_population,
          class_name,
          sep = "__"
        )

        class_enrichment_tests[[test_key]] <- run_class_enrichment(
          classified_results[[dataset_name]],
          dataset_name,
          direction,
          class_name,
          class_ids,
          reference_population
        )
      }
    }
  }
}

class_enrichment_summary <- do.call(
  rbind,
  lapply(class_enrichment_tests, function(test) test$summary)
)
row.names(class_enrichment_summary) <- NULL
class_enrichment_summary$class_enrichment_FDR <- NA_real_

# Each reference population is treated as its own multiple-testing family.
for (reference_population in class_reference_populations) {
  family_row <- class_enrichment_summary$reference_population ==
    reference_population
  valid_row <- family_row & !is.na(class_enrichment_summary$p_value)
  class_enrichment_summary$class_enrichment_FDR[valid_row] <- p.adjust(
    class_enrichment_summary$p_value[valid_row],
    method = "BH"
  )
}

class_enrichment_summary$class_hit_ids <- vapply(
  class_enrichment_tests,
  function(test) collapse_values(test$class_hits),
  character(1)
)
class_enrichment_summary$class_hit_genes <- vapply(
  class_enrichment_tests,
  function(test) get_hit_annotation(test$class_hits, "Gene_name"),
  character(1)
)

class_contingency_tables <- do.call(
  rbind,
  lapply(class_enrichment_tests, function(test) {
    test$contingency_long
  })
)
row.names(class_contingency_tables) <- NULL


# -----------------------------------------------------------------------------
# 9. Annotation-rich gene-level summaries
# -----------------------------------------------------------------------------

make_annotated_flavoprotein_summary <- function(
  results,
  dataset_name
) {
  result_match <- match(flavoproteins$Rv_ID, results$gene_id)
  matched_deseq2 <- results[result_match, , drop = FALSE]
  names(matched_deseq2) <- paste0("DESeq2_", names(matched_deseq2))

  annotated_summary <- data.frame(
    dataset = rep(dataset_name, nrow(flavoproteins)),
    dataset_label = rep(
      dataset_labels[[dataset_name]],
      nrow(flavoproteins)
    ),
    present_in_DESeq2 = !is.na(result_match),
    flavoproteins,
    matched_deseq2,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  row.names(annotated_summary) <- NULL
  annotated_summary
}

annotated_flavoprotein_summary <- do.call(
  rbind,
  lapply(names(classified_results), function(dataset_name) {
    make_annotated_flavoprotein_summary(
      classified_results[[dataset_name]],
      dataset_name
    )
  })
)
row.names(annotated_flavoprotein_summary) <- NULL

significant_hit_row <- !is.na(
  annotated_flavoprotein_summary$DESeq2_is_significant
) & annotated_flavoprotein_summary$DESeq2_is_significant

significant_flavoprotein_hits <- annotated_flavoprotein_summary[
  significant_hit_row,
  ,
  drop = FALSE
]

if (nrow(significant_flavoprotein_hits) > 0) {
  significant_flavoprotein_hits <- significant_flavoprotein_hits[
    order(
      significant_flavoprotein_hits$dataset,
      significant_flavoprotein_hits$DESeq2_expression_class,
      significant_flavoprotein_hits$flavin_class,
      significant_flavoprotein_hits$DESeq2_padj,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}


# -----------------------------------------------------------------------------
# 10. Write tables
# -----------------------------------------------------------------------------

write_csv <- function(data, filename, directory = table_directory) {
  output_path <- file.path(directory, filename)
  write.csv(
    data,
    output_path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
  log_message("Wrote table: ", output_path)
  invisible(output_path)
}

write_csv(
  classification_summary,
  "differential_expression_classification_summary.csv"
)
write_csv(matching_summary, "input_gene_matching_summary.csv")
write_csv(unmatched_flavoproteins, "unmatched_flavoproteins.csv")
write_csv(
  overall_enrichment_summary,
  "overall_flavoprotein_enrichment_summary.csv"
)
write_csv(
  overall_contingency_tables,
  "overall_flavoprotein_contingency_tables.csv"
)
write_csv(
  class_enrichment_summary,
  "flavin_class_enrichment_summary.csv"
)
write_csv(
  class_contingency_tables,
  "flavin_class_contingency_tables.csv"
)
write_csv(
  annotated_flavoprotein_summary,
  "all_flavoproteins_annotated_by_DESeq2.csv"
)
write_csv(
  significant_flavoprotein_hits,
  "significant_flavoprotein_hits_annotated.csv"
)

# Export every classified DESeq2 result with its flavoprotein status.
for (dataset_name in names(classified_results)) {
  result <- classified_results[[dataset_name]]
  annotation_match <- match(result$gene_id, flavoproteins$Rv_ID)
  uniprot_annotation <- flavoproteins[
    annotation_match,
    ,
    drop = FALSE
  ]
  names(uniprot_annotation) <- paste0(
    "UniProt_",
    names(uniprot_annotation)
  )

  annotated_result <- data.frame(
    result,
    uniprot_annotation,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write_csv(
    annotated_result,
    paste0(dataset_name, "_classified_annotated.csv"),
    directory = classified_directory
  )
}


# -----------------------------------------------------------------------------
# 11. Plotting helpers
# -----------------------------------------------------------------------------

save_plot <- function(plot, filename_stem, width, height) {
  width <- min(max(width, 7), 30)
  height <- min(max(height, 4.5), 30)

  png_path <- file.path(plot_directory, paste0(filename_stem, ".png"))
  ggplot2::ggsave(
    png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  log_message("Wrote plot: ", png_path)

  if (isTRUE(config$save_pdf_plots)) {
    pdf_path <- file.path(
      plot_directory,
      paste0(filename_stem, ".pdf")
    )
    ggplot2::ggsave(
      pdf_path,
      plot = plot,
      width = width,
      height = height,
      bg = "white"
    )
    log_message("Wrote plot: ", pdf_path)
  }
}

add_significance_label <- function(fdr_values) {
  ifelse(
    is.na(fdr_values),
    "Not testable",
    ifelse(
      fdr_values <= config$enrichment_fdr_cutoff,
      paste0("FDR <= ", config$enrichment_fdr_cutoff),
      paste0("FDR > ", config$enrichment_fdr_cutoff)
    )
  )
}

direction_colors <- c(
  downregulated = "#2166AC",
  upregulated = "#B2182B"
)

significance_shapes <- setNames(
  c(19, 1, 4),
  c(
    paste0("FDR <= ", config$enrichment_fdr_cutoff),
    paste0("FDR > ", config$enrichment_fdr_cutoff),
    "Not testable"
  )
)


# -----------------------------------------------------------------------------
# 12. Overall enrichment plots
# -----------------------------------------------------------------------------

overall_plot_data <- overall_enrichment_summary[
  !is.na(overall_enrichment_summary$fold_enrichment),
  ,
  drop = FALSE
]

if (nrow(overall_plot_data) > 0) {
  overall_plot_data$significance <- add_significance_label(
    overall_plot_data$enrichment_FDR
  )
  overall_plot_data$dataset_label <- factor(
    overall_plot_data$dataset_label,
    levels = rev(unique(overall_plot_data$dataset_label))
  )

  overall_enrichment_plot <- ggplot2::ggplot(
    overall_plot_data,
    ggplot2::aes(
      x = dataset_label,
      y = fold_enrichment,
      color = direction,
      shape = significance,
      size = selected_flavoproteins
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::geom_point(
      alpha = 0.9,
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = direction_colors) +
    ggplot2::scale_shape_manual(values = significance_shapes) +
    ggplot2::scale_size_continuous(range = c(2.5, 9)) +
    ggplot2::labs(
      title = "Flavoprotein enrichment",
      subtitle = paste0(
        "DE genes: adjusted p-value <= ",
        config$padj_cutoff,
        " and |log2 fold change| >= ",
        config$log2fc_cutoff
      ),
      x = NULL,
      y = "Fold enrichment",
      color = "Direction",
      shape = "Enrichment test",
      size = "Observed\nflavoproteins"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  save_plot(
    overall_enrichment_plot,
    "overall_flavoprotein_enrichment",
    width = 10,
    height = 3.5 + 0.5 * length(unique(overall_plot_data$dataset))
  )

  observed_expected_plot_data <- rbind(
    data.frame(
      dataset_label = overall_plot_data$dataset_label,
      direction = overall_plot_data$direction,
      count_type = "Observed",
      flavoprotein_count = overall_plot_data$selected_flavoproteins
    ),
    data.frame(
      dataset_label = overall_plot_data$dataset_label,
      direction = overall_plot_data$direction,
      count_type = "Expected",
      flavoprotein_count = overall_plot_data$expected_flavoproteins
    )
  )

  observed_expected_plot <- ggplot2::ggplot(
    observed_expected_plot_data,
    ggplot2::aes(
      x = dataset_label,
      y = flavoprotein_count,
      fill = count_type
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.72
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~direction,
      scales = "free_y"
    ) +
    ggplot2::scale_fill_manual(
      values = c(Observed = "#5B8DB8", Expected = "#CFCFCF")
    ) +
    ggplot2::labs(
      title = "Observed and expected significant flavoproteins",
      x = NULL,
      y = "Flavoprotein count",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )

  save_plot(
    observed_expected_plot,
    "overall_flavoproteins_observed_vs_expected",
    width = 11,
    height = 4 + 0.55 * length(unique(overall_plot_data$dataset))
  )
} else {
  log_message("No testable overall enrichment results were available to plot.")
}


# -----------------------------------------------------------------------------
# 13. Flavin-class enrichment plots
# -----------------------------------------------------------------------------

reference_plot_titles <- c(
  all_testable_genes = "Flavin-class enrichment among all testable genes",
  flavoproteins_only = "Flavin-class enrichment within flavoproteins"
)

reference_plot_stems <- c(
  all_testable_genes = "flavin_class_enrichment_all_genes",
  flavoproteins_only = "flavin_class_enrichment_within_flavoproteins"
)

for (reference_population in class_reference_populations) {
  class_plot_data <- class_enrichment_summary[
    class_enrichment_summary$reference_population ==
      reference_population &
      class_enrichment_summary$universe_class_genes >=
        config$minimum_class_size_for_plot &
      !is.na(class_enrichment_summary$expected_class_genes),
    ,
    drop = FALSE
  ]

  if (nrow(class_plot_data) == 0) {
    log_message(
      "No flavin-class results met the plotting threshold for ",
      reference_population,
      "."
    )
    next
  }

  class_plot_data$comparison <- paste(
    class_plot_data$dataset_label,
    class_plot_data$direction,
    sep = " | "
  )

  # A 0.5 continuity correction is used only for plot colors so that classes
  # with zero observed hits remain visible. Reported tests are unchanged.
  class_plot_data$plot_log2_enrichment <- log2(
    (class_plot_data$selected_class_genes + 0.5) /
      (class_plot_data$expected_class_genes + 0.5)
  )
  class_plot_data$significance <- add_significance_label(
    class_plot_data$class_enrichment_FDR
  )

  class_enrichment_plot <- ggplot2::ggplot(
    class_plot_data,
    ggplot2::aes(
      x = flavin_class,
      y = comparison,
      color = plot_log2_enrichment,
      size = selected_class_genes,
      shape = significance
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "grey92",
      high = "#B2182B",
      midpoint = 0
    ) +
    ggplot2::scale_shape_manual(values = significance_shapes) +
    ggplot2::scale_size_continuous(range = c(2.5, 10)) +
    ggplot2::labs(
      title = reference_plot_titles[[reference_population]],
      subtitle = paste0(
        "Classes with at least ",
        config$minimum_class_size_for_plot,
        " genes in the reference population"
      ),
      x = "Flavin class",
      y = NULL,
      color = "Stabilized log2\nenrichment",
      size = "Observed\nclass genes",
      shape = "Enrichment test"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  save_plot(
    class_enrichment_plot,
    reference_plot_stems[[reference_population]],
    width = 7 + 0.65 * length(unique(class_plot_data$flavin_class)),
    height = 3 + 0.5 * length(unique(class_plot_data$comparison))
  )
}


# -----------------------------------------------------------------------------
# 14. Publication plots of up- and downregulated flavoprotein genes
# -----------------------------------------------------------------------------

# A diverging count plot summarizes the number of significant up- and
# downregulated genes in every flavin class. The underlying table also retains
# observed, expected, fold-enrichment, and FDR values.
flavin_class_deg_counts <- class_enrichment_summary[
  class_enrichment_summary$reference_population == "all_testable_genes",
  c(
    "dataset",
    "dataset_label",
    "direction",
    "flavin_class",
    "universe_genes",
    "universe_class_genes",
    "selected_genes",
    "selected_class_genes",
    "expected_class_genes",
    "fold_enrichment",
    "odds_ratio",
    "p_value",
    "class_enrichment_FDR",
    "class_hit_ids",
    "class_hit_genes"
  ),
  drop = FALSE
]

names(flavin_class_deg_counts)[
  names(flavin_class_deg_counts) == "selected_class_genes"
] <- "differentially_expressed_class_genes"

flavin_class_deg_counts$signed_gene_count <- ifelse(
  flavin_class_deg_counts$direction == "downregulated",
  -flavin_class_deg_counts$differentially_expressed_class_genes,
  flavin_class_deg_counts$differentially_expressed_class_genes
)

write_csv(
  flavin_class_deg_counts,
  "flavin_class_up_down_DEG_summary.csv"
)

class_count_totals <- aggregate(
  differentially_expressed_class_genes ~ flavin_class,
  data = flavin_class_deg_counts,
  FUN = sum
)

classes_with_hits <- class_count_totals$flavin_class[
  class_count_totals$differentially_expressed_class_genes > 0
]

class_count_plot_data <- flavin_class_deg_counts[
  flavin_class_deg_counts$flavin_class %in% classes_with_hits,
  ,
  drop = FALSE
]

if (nrow(class_count_plot_data) > 0) {
  class_order <- class_count_totals$flavin_class[
    order(class_count_totals$differentially_expressed_class_genes)
  ]
  class_count_plot_data$flavin_class <- factor(
    class_count_plot_data$flavin_class,
    levels = class_order
  )

  class_direction_count_plot <- ggplot2::ggplot(
    class_count_plot_data,
    ggplot2::aes(
      x = signed_gene_count,
      y = flavin_class,
      fill = direction
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey30",
      linewidth = 0.45
    ) +
    ggplot2::geom_col(
      width = 0.72,
      color = "grey25",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(
          differentially_expressed_class_genes > 0,
          differentially_expressed_class_genes,
          ""
        ),
        hjust = ifelse(signed_gene_count >= 0, -0.2, 1.2)
      ),
      size = 3.2,
      color = "grey15"
    ) +
    ggplot2::facet_wrap(
      ~dataset_label,
      scales = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = direction_colors,
      breaks = c("downregulated", "upregulated"),
      labels = c(
        downregulated = "DOWN",
        upregulated = "UP"
      )
    ) +
    ggplot2::scale_x_continuous(
      labels = function(value) abs(value),
      expand = ggplot2::expansion(mult = c(0.12, 0.12))
    ) +
    ggplot2::labs(
      title = "Significant flavoprotein DEGs by flavin class",
      subtitle = paste0(
        "Adjusted p-value <= ",
        config$padj_cutoff,
        " and |log2 fold change| >= ",
        config$log2fc_cutoff
      ),
      x = "Number of genes (DOWN <- 0 -> UP)",
      y = "Flavin class",
      fill = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top"
    )

  save_plot(
    class_direction_count_plot,
    "flavin_class_up_down_DEG_counts",
    width = 6 + 4.5 * length(unique(class_count_plot_data$dataset)),
    height = 3.5 + 0.42 * length(classes_with_hits)
  )
} else {
  log_message(
    "No significant flavoprotein DEGs were available for the class count plot."
  )
}


# Ranked diverging log2-fold-change plots reproduce the visual style of a
# candidate-gene bar chart. One overview and one flavin-class-faceted figure
# are created for each DESeq2 comparison.
ranked_gene_tables <- lapply(
  names(classified_results),
  function(dataset_name) {
    hits <- significant_flavoprotein_hits[
      significant_flavoprotein_hits$dataset == dataset_name &
        !is.na(significant_flavoprotein_hits$DESeq2_log2FoldChange),
      ,
      drop = FALSE
    ]

    if (nrow(hits) == 0) {
      return(NULL)
    }

    hits$absolute_log2_fold_change <- abs(
      hits$DESeq2_log2FoldChange
    )
    hits <- hits[
      order(
        -hits$absolute_log2_fold_change,
        hits$DESeq2_padj,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    hits$publication_rank <- seq_len(nrow(hits))
    hits$included_in_ranked_plot <- hits$publication_rank <=
      config$maximum_genes_per_ranked_plot

    has_gene_name <- !is.na(hits$Gene_name) &
      trimws(hits$Gene_name) != ""
    hits$plot_gene_label <- ifelse(
      has_gene_name,
      paste0(hits$Gene_name, " (", hits$Rv_ID, ")"),
      hits$Rv_ID
    )

    hits
  }
)

ranked_gene_tables <- Filter(Negate(is.null), ranked_gene_tables)

if (length(ranked_gene_tables) > 0) {
  ranked_flavoprotein_degs <- do.call(rbind, ranked_gene_tables)
  row.names(ranked_flavoprotein_degs) <- NULL

  write_csv(
    ranked_flavoprotein_degs,
    "publication_ranked_flavoprotein_DEGs.csv"
  )

  sanitize_filename <- function(value) {
    value <- gsub("[^A-Za-z0-9._-]+", "_", value)
    substr(value, 1, 120)
  }

  for (dataset_name in names(classified_results)) {
    plot_data <- ranked_flavoprotein_degs[
      ranked_flavoprotein_degs$dataset == dataset_name &
        ranked_flavoprotein_degs$included_in_ranked_plot,
      ,
      drop = FALSE
    ]

    if (nrow(plot_data) == 0) {
      next
    }

    plot_data <- plot_data[
      order(plot_data$DESeq2_log2FoldChange),
      ,
      drop = FALSE
    ]
    plot_data$plot_gene_label <- factor(
      plot_data$plot_gene_label,
      levels = plot_data$plot_gene_label
    )

    genes_shown <- nrow(plot_data)
    all_significant_hits <- sum(
      ranked_flavoprotein_degs$dataset == dataset_name
    )
    subtitle_text <- paste0(
      "Significant flavoproteins ranked by |log2 fold change|; ",
      genes_shown,
      " of ",
      all_significant_hits,
      " shown"
    )

    ranked_gene_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = DESeq2_log2FoldChange,
        y = plot_gene_label,
        fill = DESeq2_expression_class
      )
    ) +
      ggplot2::geom_vline(
        xintercept = 0,
        color = "grey25",
        linewidth = 0.5
      ) +
      ggplot2::geom_col(
        width = 0.78,
        color = "grey25",
        linewidth = 0.2
      ) +
      ggplot2::scale_fill_manual(
        values = direction_colors,
        breaks = c("downregulated", "upregulated"),
        labels = c(
          downregulated = "DOWN",
          upregulated = "UP"
        ),
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.04, 0.06))
      ) +
      ggplot2::labs(
        title = paste0(
          "Significant flavoprotein genes: ",
          dataset_labels[[dataset_name]]
        ),
        subtitle = subtitle_text,
        x = "log2 fold change",
        y = NULL,
        fill = "Direction"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    filename_label <- sanitize_filename(dataset_labels[[dataset_name]])
    ranked_plot_height <- 3.5 + 0.29 * genes_shown

    save_plot(
      ranked_gene_plot,
      paste0("ranked_flavoprotein_log2fc_", filename_label),
      width = 12,
      height = ranked_plot_height
    )

    ranked_gene_class_plot <- ranked_gene_plot +
      ggplot2::facet_grid(
        flavin_class ~ .,
        scales = "free_y",
        space = "free_y",
        switch = "y"
      ) +
      ggplot2::labs(
        title = paste0(
          "Significant flavoprotein genes by flavin class: ",
          dataset_labels[[dataset_name]]
        ),
        subtitle = subtitle_text
      ) +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        strip.text.y.left = ggplot2::element_text(
          angle = 0,
          face = "bold"
        )
      )

    save_plot(
      ranked_gene_class_plot,
      paste0(
        "ranked_flavoprotein_log2fc_by_class_",
        filename_label
      ),
      width = 13,
      height = ranked_plot_height +
        0.35 * length(unique(plot_data$flavin_class))
    )
  }
} else {
  ranked_flavoprotein_degs <- significant_flavoprotein_hits
  log_message(
    "No significant flavoprotein DEGs were available for ranked gene plots."
  )
}


# -----------------------------------------------------------------------------
# 15. Reproducibility files and output guide
# -----------------------------------------------------------------------------

parameters <- data.frame(
  parameter = c(
    "input_directory",
    "deseq2_pattern",
    "flavoprotein_pattern",
    "padj_cutoff",
    "log2fc_cutoff",
    "enrichment_fdr_cutoff",
    "minimum_class_size_for_plot",
    "maximum_genes_per_ranked_plot",
    "class_test_alternative",
    "multiple_testing_method"
  ),
  value = c(
    input_directory,
    config$deseq2_pattern,
    config$flavoprotein_pattern,
    config$padj_cutoff,
    config$log2fc_cutoff,
    config$enrichment_fdr_cutoff,
    config$minimum_class_size_for_plot,
    config$maximum_genes_per_ranked_plot,
    "greater",
    "Benjamini-Hochberg"
  ),
  stringsAsFactors = FALSE
)
write_csv(parameters, "analysis_parameters.csv")

input_manifest <- data.frame(
  input_type = c(
    rep("DESeq2", length(deseq2_files)),
    "flavoprotein_annotations"
  ),
  filename = c(basename(deseq2_files), basename(flavoprotein_file)),
  absolute_path = c(
    normalizePath(deseq2_files),
    normalizePath(flavoprotein_file)
  ),
  file_size_bytes = c(
    file.info(deseq2_files)$size,
    file.info(flavoprotein_file)$size
  ),
  stringsAsFactors = FALSE
)
write_csv(input_manifest, "input_file_manifest.csv")

session_info_path <- file.path(output_directory, "sessionInfo.txt")
session_connection <- file(session_info_path, open = "wt")
writeLines(capture.output(sessionInfo()), session_connection)
close(session_connection)
log_message("Wrote reproducibility file: ", session_info_path)

results_guide <- c(
  "FLAVOPROTEIN ENRICHMENT OUTPUT GUIDE",
  "",
  "Primary result tables:",
  "- tables/overall_flavoprotein_enrichment_summary.csv",
  "  Overall enrichment of flavoproteins among up- and downregulated genes.",
  "- tables/flavin_class_enrichment_summary.csv",
  "  Class-specific tests using both all genes and flavoproteins-only references.",
  "- tables/significant_flavoprotein_hits_annotated.csv",
  "  Significant hits with every UniProt, KEGG, and GO annotation retained.",
  "- tables/all_flavoproteins_annotated_by_DESeq2.csv",
  "  All flavoproteins, including absent and nonsignificant genes, with DESeq2 data.",
  "- tables/flavin_class_up_down_DEG_summary.csv",
  "  Up/down class counts with expected counts, enrichment statistics, and hit IDs.",
  "- tables/publication_ranked_flavoprotein_DEGs.csv",
  "  Ranked gene-level source data for the diverging log2-fold-change figures.",
  "",
  "Supporting tables:",
  "- tables/input_gene_matching_summary.csv",
  "- tables/unmatched_flavoproteins.csv",
  "- tables/differential_expression_classification_summary.csv",
  "- tables/overall_flavoprotein_contingency_tables.csv",
  "- tables/flavin_class_contingency_tables.csv",
  "- classified_deseq2/*_classified_annotated.csv",
  "",
  "Publication figures:",
  "- plots/flavin_class_up_down_DEG_counts.*",
  "  Diverging up/down gene counts broken down by flavin class.",
  "- plots/ranked_flavoprotein_log2fc_*.*",
  "  Image-style ranked diverging log2-fold-change plots for significant hits.",
  "- plots/ranked_flavoprotein_log2fc_by_class_*.*",
  "  Ranked gene plots separated into flavin-class panels.",
  "- plots/overall_flavoprotein_enrichment.*",
  "- plots/overall_flavoproteins_observed_vs_expected.*",
  "- plots/flavin_class_enrichment_all_genes.*",
  "- plots/flavin_class_enrichment_within_flavoproteins.*",
  "",
  "Interpretation:",
  "- Fold enrichment > 1 indicates more hits than expected.",
  "- Odds ratio > 1 supports enrichment.",
  "- enrichment_FDR and class_enrichment_FDR are BH-adjusted p-values.",
  "- all_testable_genes asks whether a class is enriched among all DE genes.",
  "- flavoproteins_only asks whether a class is enriched within DE flavoproteins.",
  "- Flavin class labels are used exactly as provided; they are not consolidated.",
  "- Classes below the plot-size threshold remain in the CSV result tables.",
  "",
  paste0(
    "Thresholds: padj <= ",
    config$padj_cutoff,
    ", |log2FC| >= ",
    config$log2fc_cutoff,
    ", enrichment FDR <= ",
    config$enrichment_fdr_cutoff
  )
)

guide_path <- file.path(output_directory, "README_results.txt")
writeLines(results_guide, guide_path, useBytes = TRUE)
log_message("Wrote output guide: ", guide_path)

run_finished <- Sys.time()
log_message(
  "Analysis completed successfully in ",
  round(as.numeric(difftime(run_finished, run_started, units = "secs")), 2),
  " seconds."
)
log_message("Results directory: ", normalizePath(output_directory))

run_log_path <- file.path(output_directory, "run_log.txt")
writeLines(run_log, run_log_path, useBytes = TRUE)
