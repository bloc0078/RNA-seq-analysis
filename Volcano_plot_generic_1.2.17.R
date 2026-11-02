# ==============================================================================
# SCRIPT: Volcano_plot_generic_1.2.17.R
# DESCRIPTION:
#   Generic volcano plotter for standardized DESeq2 RNA-seq and compatible
#   volcano-style input tables.
#
# VERSION:
#   v1.2.17
#
# CHANGES IN v1.2.17:
#   - Bold subscript digits and explicitly use Arial for the spaced plot title.
#   - Use the macOS Quartz PDF device to support the Arial title font.
#   - Bold the y-axis mathematical label; use plain axis tick values.
#   - Hide the manual-label caption by default with a dedicated setting.
#   - Reduce gray background dots to 1.4; keep DEG dots at 1.8.
#   - Retain the manual-label significance toggle and user-defined title.
#
# PURPOSE:
#   Generate publication-ready volcano plots with:
#     - automatic top upregulated and downregulated DEG labeling
#     - optional manually selected labels
#     - grouped manual label coloring for enrichment-linked genes
#     - support for KEGG ORA, GO ORA, KEGG + GO ORA, and Gene of interest labels
#     - labeled-point highlighting
#     - exported plotting-ready, significant-hit, and labeled-feature tables
#     - run logs and session information
#
# CURRENT LABELING DESIGN:
#   Automatic labels:
#     - Top N upregulated DEGs
#     - Top N downregulated DEGs
#     - selected using the ranking method defined in Section 1.6
#
#   Manual labels:
#     - defined in Section 1.6B using manual_label_group_table
#     - match_value should preferably use gene_id from the DESeq2 input table
#     - display_label controls the text shown on the plot
#     - label_group controls the manual label color/category
#
#   Recommended manual label groups:
#     - KEGG ORA
#     - GO ORA
#     - KEGG + GO ORA
#     - Gene of interest
#
#   Notes:
#     - Up/down direction does not need to be encoded in the manual label group
#       because volcano plot position already indicates direction.
#     - Genes found in both KEGG and GO ORA should be listed once as
#       "KEGG + GO ORA".
#     - Genes that are biologically important but not KEGG/GO ORA-linked can be
#       listed as "Gene of interest".
#
# MANUAL LABEL MATCHING:
#   Manual labels are matched against available standardized columns including:
#     - gene_id
#     - gene_name
#     - feature_label
#
#   Best practice:
#     - use gene_id as match_value for reliable matching
#     - use display_label for the clean plotted label
#
# EXAMPLE:
#   match_value   = "ACP3PV_RS07545"
#   display_label = "ribA2"
#   label_group   = "Gene of interest"
#
# VISUAL DESIGN:
#   - Volcano point colors show significance category.
#   - Labeled genes are slightly enlarged/highlighted.
#   - Automatic DEG labels use standard black text.
#   - Manual grouped labels use fixed text colors by label_group.
#   - Manual label color definitions are described in the plot caption.
#   - Manual label colors are applied without adding a second ggplot color scale,
#     preserving the original volcano point coloring.
#
# IMPORTANT CHANGE CONTROL NOTES:
#   - Section 1 contains user-editable settings only.
#   - Manual label groups should be edited in Section 1.6B.
#   - Do not change output filenames, folder names, or plot theme behavior unless
#     intentionally versioning the script.
#
# ==============================================================================

# ============================================================================
# SECTION 1: USER SETTINGS
# ============================================================================
#
# Section 1 defines user-editable settings only.
# No files are detected, imported, standardized, or plotted in this section.
#
# Subsections:
#   1.1 Input file detection settings
#   1.2 Project, comparison, and output naming
#   1.3 Significance and effect-size thresholds
#   1.4 Plot appearance settings
#   1.5 Color settings
#   1.6 Automatic labeling settings
#   1.6B Manual grouped labeling settings
#   1.7 Export settings
#
# ============================================================================
cat("\n--------------------------------------------------------\n")
cat("SECTION 1: USER SETTINGS\n")
cat("--------------------------------------------------------\n")

script_version <- "v1.2.17"

# ----------------------------------------------------------------------------
# 1.1 Input file detection
# ----------------------------------------------------------------------------

# Define the input directory first because input detection and output paths use it.

input_dir <- getwd()

input_file <- NULL
# input_file <- "DESeq2_Full_Results_Treatment_vs_Control.csv"

use_auto_input_detection <- is.null(input_file)

cat("[SETTINGS] Input directory:", input_dir, "\n")

if (use_auto_input_detection) {
  cat("[SETTINGS] Input file: auto-detect in working directory\n")
} else {
  cat("[SETTINGS] Input file:", input_file, "\n")
}


# ----------------------------------------------------------------------------
# Console logging
# ----------------------------------------------------------------------------

log_full_console_output <- TRUE

if (!dir.exists(input_dir)) {
  stop(
    "input_dir does not exist:\n",
    input_dir,
    "\n\nCheck the working directory before running the script."
  )
}

# The final log path cannot be created until the run-specific output directory
# and its logs subdirectory exist. It is initialized here for the Section 1
# checkpoint and assigned its final path in Section 10.
console_log_file <- NA_character_

# ----------------------------------------------------------------------------
# 1.2 Project, comparison, and output naming
# ----------------------------------------------------------------------------

# These names are used for plot titles, output folders, and exported filenames.
#
# You may set project_name, comparison_name, and analysis_type manually.
# If comparison_name or analysis_type are left as NULL, later sections can infer
# them from the detected input file and standardized input type.

project_name <- "ribA2"

comparison_name <- NULL
# comparison_name <- "Treatment_vs_Control"

analysis_type <- "RNAseq"
# analysis_type <- "RNAseq"
# analysis_type <- "TnSeq"

# Main output directory.
# A run-specific subfolder will be created later.
output_root <- file.path(input_dir, "Volcano_Plot_Results")

# TRUE prevents accidental overwriting by creating a timestamped run folder.
use_timestamped_output <- TRUE

# Optional custom tag added to output filenames.
custom_output_tag <- NULL
# custom_output_tag <- "ribA2_IFNg"

cat("[SETTINGS] Project name:", project_name, "\n")

if (is.null(comparison_name)) {
  cat("[SETTINGS] Comparison name: infer later from input filename\n")
} else {
  cat("[SETTINGS] Comparison name:", comparison_name, "\n")
}

if (is.null(analysis_type)) {
  cat("[SETTINGS] Analysis type: infer later from detected input type\n")
} else {
  cat("[SETTINGS] Analysis type:", analysis_type, "\n")
}

cat("[SETTINGS] Output root:", output_root, "\n")
cat("[SETTINGS] Timestamped output folder:", use_timestamped_output, "\n")


# ----------------------------------------------------------------------------
# 1.3 Significance and effect-size thresholds
# ----------------------------------------------------------------------------

# Primary significance threshold.
# For DESeq2 RNA-seq results, padj is usually preferred.
# For Tn-seq results, this may correspond to q-value, FDR, or adjusted p-value
# after standardization.

padj_cutoff <- 0.05

# Raw p-value fallback.
# This is used only if an adjusted p-value column is unavailable.

pvalue_cutoff <- 0.05

# Effect-size threshold.
# For RNA-seq, this is log2 fold change.
# For Tn-seq, this may represent log2 fitness, log2 enrichment/depletion,
# or another standardized log2-scale effect size.

log2fc_cutoff <- 1

# Which significance column should be preferred?
# Options:
#   "auto"   = use adjusted p-value if available; otherwise use raw p-value
#   "padj"   = require adjusted p-value
#   "pvalue" = use raw p-value

significance_column_preference <- "auto"

# Should points with significant p-values but small effect sizes be shown as a
# separate category?

classify_significant_low_effect <- TRUE

cat("[SETTINGS] Adjusted p-value cutoff:", padj_cutoff, "\n")
cat("[SETTINGS] Raw p-value cutoff:", pvalue_cutoff, "\n")
cat("[SETTINGS] Absolute log2FC/effect cutoff:", log2fc_cutoff, "\n")
cat("[SETTINGS] Significance column preference:", significance_column_preference, "\n")
cat("[SETTINGS] Classify significant low-effect points:", classify_significant_low_effect, "\n")

# ----------------------------------------------------------------------------
# 1.4 Plot appearance settings
# ----------------------------------------------------------------------------

# Main plot labels.
# Leave plot_title or plot_subtitle as NULL to build them automatically later.

plot_title <- "WT +EtOH / WT -Atc"
plot_subtitle <- ""
clean_inferred_plot_title <- TRUE

auto_italicize_title_gene <- TRUE
title_gene_to_italic <- "ribA2"

# Optional custom display name for the comparison.
# If NULL, Section 10 will infer the comparison name from the input filename.
comparison_display_name <- "WT +EtOH / WT -Atc"
# comparison_display_name <- "ribA2 Tet-OFF + RF, day 0 vs WT"

# Caption/footnote settings.

show_input_file_in_caption <- FALSE
show_rows_plotted_in_caption <- FALSE
show_color_explanation_in_caption <- FALSE
show_manual_label_key_in_caption <- FALSE

caption_color_explanation <- 
  "Colored points pass both cutoffs; gray points do not."

x_axis_label <- expression(bold(log["2"]~"fold change"))
y_axis_label <- expression(
  bold(-log["10"]("adjusted " * bolditalic(p) * "-value"))
)

# Point appearance

point_size <- 1.8
# Gray points that do not pass both cutoffs. Labeled highlights remain larger.
background_point_size <- 1.2
point_alpha <- 0.8

background_point_alpha <- 0.25
deg_point_alpha <- 0.85

# Label appearance

label_size <- 5
label_fontface <- "italic"
label_box_padding <- 0.55
label_point_padding <- 0.55
label_max_overlaps <- Inf
label_text_color <- "black"
label_segment_color <- "gray40"
label_force <- 10
label_force_pull <- 0.05
label_min_segment_length <- 0

# Manual grouped label appearance.
# These apply only to manually selected grouped labels such as KEGG ORA,
# GO ORA, KEGG + GO ORA, and Gene of interest.

manual_label_size <- label_size + 0.3
manual_label_fontface <- "italic"

# Labeled-point highlighting
highlight_labeled_points <- TRUE
labeled_point_size <- 2.8
labeled_point_stroke <- 0.4
labeled_point_outline_color <- "black"
labeled_point_fill_mode <- "category"


# Label connector lines

label_min_segment_length <- 0

# Threshold line appearance

show_threshold_lines <- TRUE
threshold_line_type <- "dashed"
threshold_line_width <- 0.4

# Plot theme

base_font_size <- 12

# Optional x-axis limits.
# Leave as NULL to let ggplot choose automatically.
# Example:
# x_axis_limits <- c(-6, 6)

x_axis_limits <- c(-6, 6)

# Optional y-axis display cap.
# Leave y_axis_limit as NULL to let ggplot choose automatically.
#
# If y_axis_cap_display_values is TRUE and y_axis_limit is not NULL,
# points above y_axis_limit are displayed at the cap value.
# This does not alter pvalue, padj, significance_value_raw, or
# neg_log10_significance in the exported tables.
#
# Examples:
#   y_axis_limit <- NULL  # no cap
#   y_axis_limit <- 15    # useful for lower-significance plots
#   y_axis_limit <- 20
#   y_axis_limit <- 150   # useful when padj values are extremely small

y_axis_limit <- 150
y_axis_cap_display_values <- TRUE

# Plot theme

base_font_size <- 15

# Axis line and tick mark appearance

show_axis_lines <- TRUE
show_axis_ticks <- TRUE

axis_line_color <- "black"
axis_line_width <- 0.5

axis_tick_color <- "black"
axis_tick_width <- 0.8
axis_tick_length <- 7

cat("[SETTINGS] Plot title:", ifelse(is.null(plot_title), "build automatically later", plot_title), "\n")
cat("[SETTINGS] Plot subtitle:", ifelse(is.null(plot_subtitle), "build automatically later", plot_subtitle), "\n")
cat("[SETTINGS] Point size:", point_size, "\n")
cat("[SETTINGS] Point alpha:", point_alpha, "\n")
cat("[SETTINGS] Label size:", label_size, "\n")
cat("[SETTINGS] Show threshold lines:", show_threshold_lines, "\n")

# ----------------------------------------------------------------------------
# 1.5 Color settings
# ----------------------------------------------------------------------------

# Colors used for volcano plot significance categories.
#
# Default categories expected after classification:
#   Up
#   Down
#   Significant_low_effect
#   Not_significant
#
# These colors can be changed later if needed for publication or accessibility.

# ATC-treated comparison palette
color_up <- "#D33682"      # orchid / rose-magenta
color_down <- "#008CBA"    # teal-blue

# Fixed gray categories.
# Only Up and Down colors are intended to be user-edited in this section.
color_significant_low_effect <- "gray60"
color_not_significant <- "gray72"

# Untreated / no-ATC comparison palette
# color_up <- "#F28E2B"    # warm orange
# color_down <- "#6F4E7C"  # muted plum

# Store colors in a named vector for ggplot scale_color_manual() later.

volcano_colors <- c(
  "Up" = color_up,
  "Down" = color_down,
  "Significant_low_effect" = color_significant_low_effect,
  "Not_significant" = color_not_significant
)

# Should the detailed color legend be shown?
# TRUE  = show Up, Down, Significant_low_effect, Not_significant legend
# FALSE = keep the colors on the plot but hide the detailed color legend

show_detailed_category_legend <- FALSE

cat("[SETTINGS] Show detailed category legend:", show_detailed_category_legend, "\n")

# ----------------------------------------------------------------------------
# 1.6 Labeling settings
# ----------------------------------------------------------------------------

# Labeling strategy:
#   1. By default, label top significant features only.
#   2. Do not force-label any dataset-specific genes/features.
#   3. Keep labels limited enough to avoid overcrowding.

# Should the script automatically label top significant hits?

label_top_hits <- TRUE

# Number of top upregulated/enriched features to label.

top_n_up <- 5

# Number of top downregulated/depleted features to label.

top_n_down <- 5

# Ranking method for automatic labels.
# Options:
#   "padj_then_log2FC" = prioritize smallest adjusted p-value, then largest effect
#   "log2FC_then_padj" = prioritize largest directional log2FC effect, then smallest adjusted p-value
#   "combined_score"  = prioritize features with both strong significance and effect
#
# Recommended for volcano labels:
#   "log2FC_then_padj"
#
# Rationale:
#   The volcano y-axis already displays significance. For labels, the most
#   visually useful genes are usually the strongest effect-size hits among DEGs.

label_ranking_method <- "log2FC_then_padj"

# Should labels be restricted to features that pass both the significance and
# effect-size thresholds?

label_only_significant <- TRUE

# Should labels use gene symbols/names when available instead of locus IDs?

prefer_gene_name_labels <- TRUE

cat("[SETTINGS] Label top hits:", label_top_hits, "\n")
cat("[SETTINGS] Top up labels:", top_n_up, "\n")
cat("[SETTINGS] Top down labels:", top_n_down, "\n")
cat("[SETTINGS] Label ranking method:", label_ranking_method, "\n")
cat("[SETTINGS] Label only significant hits:", label_only_significant, "\n")
cat("[SETTINGS] Prefer gene-name labels:", prefer_gene_name_labels, "\n")

# ----------------------------------------------------------------------------
# 1.6B Manual grouped label settings
# ----------------------------------------------------------------------------

# Optional user-defined labels grouped by biological/analysis category.
#
# Use this section when specific genes/features should be forced onto the plot
# because they connect to enrichment results, pathway plots, candidate gene
# triage, or another analysis-specific reason.
#
# Manual labels are matched against:
#   gene_id
#   gene_name
#   feature_label
#
# match_value:
#   The gene_id, gene_name, or feature_label to search for.
#
# display_label:
#   The label text shown on the volcano plot.
#
# label_group:
#   The biological/analysis group used for label styling.
#
# Examples of useful label groups:
#   "KEGG ORA down"
#   "GO ORA down"
#   "KEGG ORA up"
#   "GO ORA up"
#   "KEGG GSEA core down"
#   "GO GSEA core up"

manual_labels_enabled <- TRUE

manual_label_group_table <- data.frame(
  match_value = c(
    # Gene of interest
    "ACP3PV_RS00890",
    "ACP3PV_RS00895",
    "ACP3PV_RS00915",
    "ACP3PV_RS01125",
    "ACP3PV_RS04015",
    "ACP3PV_RS04020",
    "ACP3PV_RS06055",
    "ACP3PV_RS06060",
    "ACP3PV_RS06065",
    "ACP3PV_RS06070",
    "ACP3PV_RS07545",
    "ACP3PV_RS08280",
    "ACP3PV_RS08535",
    "ACP3PV_RS08540",
    "ACP3PV_RS08545",
    "ACP3PV_RS08550",
    "ACP3PV_RS12270",
    "ACP3PV_RS12395",
    "ACP3PV_RS12400",
    "ACP3PV_RS15605",
    "ACP3PV_RS15815",
    "ACP3PV_RS15820",
    "ACP3PV_RS16585",
    "ACP3PV_RS20255",
    "ACP3PV_RS20260",
    "ACP3PV_RS07550",
    "ACP3PV_RS07410",
    "ACP3PV_RS03690",
    "ACP3PV_RS03685"
  ),
  display_label = c(
    # Gene of interest
    "mce1R",
    "fadD5",
    "mce1B",
    "pckA",
    "fadE9",
    "mmsA",
    "Rv1128c",
    "Rv1129c",
    "prpD",
    "prpC",
    "ribA2",
    "frdB",
    "cydC",
    "cydD",
    "cydB",
    "cydA",
    "narK1",
    "plcC",
    "plcB",
    "fadD29",
    "leuD",
    "leuC",
    "PPE51",
    "mmpL8",
    "papA1",
    "ribH",
    "PPE20",
    "pqqE",
    "rv0692"
  ),
  label_group = c(
    # Gene of interest, not necessarily KEGG/GO ORA
    "Gene of interest"  
  ),
  
  stringsAsFactors = FALSE
)

# Label colors/fills for each manual label group.
# Edit these as needed for the specific analysis being displayed.

manual_label_group_fills <- c(
  "Gene of interest" = "white",
  "Manual highlight" = "#FFF2B2"
)

manual_label_group_text_colors <- c(
  "Gene of interest" = "black",
  "Manual highlight" = "black"
)

manual_label_group_color_names <- c(
  "Gene of interest" = "black",
  "Manual highlight" = "black"
)

# Should manually selected labels be shown even if they do not pass the
# significance/effect-size thresholds?

manual_labels_even_if_not_significant <- FALSE

# Should manual labels be added in addition to automatic top-hit labels?
# TRUE  = automatic top DEG labels + manual grouped labels
# FALSE = manual grouped labels only, if manual_labels_enabled is TRUE

manual_labels_add_to_auto_labels <- TRUE

# Optional label exclusion pattern.
# This removes labels based on the final plotted label text.
#
# This is useful when locus tags are automatically selected as labels but you
# only want named genes shown on the plot.
#
# Set exclude_acp3pv_labels <- FALSE to allow ACP3PV labels if selected
# automatically.

exclude_acp3pv_labels <- TRUE

if (exclude_acp3pv_labels) {
  exclude_label_pattern <- "^ACP3PV"
} else {
  exclude_label_pattern <- NULL
}

# Compatibility objects used by downstream manual-label code.
# Do not edit these unless the downstream script logic changes.

manual_label_features <- manual_label_group_table$match_value

manual_label_display_names <- stats::setNames(
  manual_label_group_table$display_label,
  manual_label_group_table$match_value
)

# ----------------------------------------------------------------------------
# 1.7 Export settings
# ----------------------------------------------------------------------------

# Plot export settings.
# Width and height are in inches.
#
# PNG is the default output format for easy viewing, sharing, and reuse.
# TIFF is useful for journal/publication-quality raster output.
# SVG is useful for editing vector graphics in Inkscape.
# PDF is useful as a vector review/archive copy.
#
# A 3 x 3 inch PNG or TIFF at 600 dpi exports as 1800 x 1800 pixels.

export_png <- TRUE
export_tiff <- TRUE
export_svg <- TRUE
export_pdf <- TRUE

png_width <- 8
png_height <- 8
png_dpi <- 600

tiff_width <- 8
tiff_height <- 8
tiff_dpi <- 600

svg_width <- 8
svg_height <- 8

pdf_width <- 8
pdf_height <- 8

# Table export settings.
# These are useful for checking exactly what was plotted and labeled.

save_standardized_table <- TRUE
save_significant_table <- TRUE
save_labeled_table <- TRUE

# Run documentation exports.

save_run_log <- TRUE
save_session_info <- TRUE

cat("[SETTINGS] Export PNG:", export_png, "\n")
cat("[SETTINGS] Export TIFF:", export_tiff, "\n")
cat("[SETTINGS] Export SVG:", export_svg, "\n")
cat("[SETTINGS] Export PDF:", export_pdf, "\n")

cat("[SETTINGS] PNG size:", png_width, "x", png_height, "inches\n")
cat("[SETTINGS] PNG DPI:", png_dpi, "\n")

cat("[SETTINGS] TIFF size:", tiff_width, "x", tiff_height, "inches\n")
cat("[SETTINGS] TIFF DPI:", tiff_dpi, "\n")

cat("[SETTINGS] SVG size:", svg_width, "x", svg_height, "inches\n")

cat("[SETTINGS] PDF size:", pdf_width, "x", pdf_height, "inches\n")

cat("[SETTINGS] Save standardized table:", save_standardized_table, "\n")
cat("[SETTINGS] Save significant table:", save_significant_table, "\n")
cat("[SETTINGS] Save labeled table:", save_labeled_table, "\n")
cat("[SETTINGS] Save run log:", save_run_log, "\n")
cat("[SETTINGS] Save session info:", save_session_info, "\n")

# ----------------------------------------------------------------------------
# 1.8 Section 1 checkpoint
# ----------------------------------------------------------------------------

required_section1_settings <- c(
  "script_version",
  "log_full_console_output",
  "console_log_file",
  "input_dir",
  "input_file",
  "use_auto_input_detection",
  "project_name",
  "comparison_name",
  "analysis_type",
  "output_root",
  "use_timestamped_output",
  "custom_output_tag",
  "padj_cutoff",
  "pvalue_cutoff",
  "log2fc_cutoff",
  "significance_column_preference",
  "classify_significant_low_effect",
  "plot_title",
  "auto_italicize_title_gene",
  "title_gene_to_italic",
  "plot_subtitle",
  "x_axis_label",
  "y_axis_label",
  "point_size",
  "point_alpha",
  "label_size",
  "label_fontface",
  "manual_label_size",
  "manual_label_fontface",
  "label_box_padding",
  "label_point_padding",
  "label_max_overlaps",
  "show_threshold_lines",
  "threshold_line_type",
  "threshold_line_width",
  "base_font_size",
  "x_axis_limits",
  "y_axis_limit",
  "color_up",
  "color_down",
  "color_significant_low_effect",
  "color_not_significant",
  "volcano_colors",
  "label_top_hits",
  "top_n_up",
  "top_n_down",
  "label_ranking_method",
  "label_only_significant",
  "prefer_gene_name_labels",
  "export_png",
  "export_tiff",
  "export_svg",
  "export_pdf",
  "png_width",
  "png_height",
  "png_dpi",
  "tiff_width",
  "tiff_height",
  "tiff_dpi",
  "svg_width",
  "svg_height",
  "pdf_width",
  "pdf_height",
  "save_standardized_table",
  "save_significant_table",
  "save_labeled_table",
  "save_run_log",
  "save_session_info",
  "show_detailed_category_legend",
  "comparison_display_name",
  "show_input_file_in_caption",
  "show_rows_plotted_in_caption",
  "show_color_explanation_in_caption",
  "caption_color_explanation",
  "manual_labels_enabled",
  "manual_label_group_table",
  "manual_label_group_fills",
  "manual_label_group_text_colors",
  "manual_label_group_color_names",
  "exclude_acp3pv_labels",
  "exclude_label_pattern",
  "manual_label_features",
  "manual_labels_even_if_not_significant",
  "manual_labels_add_to_auto_labels",
  "manual_label_display_names",
  "label_text_color",
  "label_segment_color",
  "label_force",
  "label_force_pull",
  "highlight_labeled_points",
  "labeled_point_size",
  "labeled_point_stroke",
  "labeled_point_outline_color",
  "label_min_segment_length",
  "clean_inferred_plot_title",
  "show_axis_lines",
  "show_axis_ticks",
  "axis_line_color",
  "axis_line_width",
  "axis_tick_color",
  "axis_tick_width",
  "axis_tick_length",
  "y_axis_cap_display_values"
)

missing_section1_settings <- required_section1_settings[
  !vapply(required_section1_settings, exists, logical(1))
]

if (length(missing_section1_settings) > 0) {
  stop(
    "SECTION 1 checkpoint failed. Missing settings:\n",
    paste(missing_section1_settings, collapse = "\n")
  )
}

cat("[CHECKPOINT] Section 1 completed successfully.\n")
cat("[CHECKPOINT] All expected Section 1 settings are defined.\n")
# Validate manual grouped label settings.

if (manual_labels_enabled) {
  
  required_manual_label_columns <- c(
    "match_value",
    "display_label",
    "label_group"
  )
  
  missing_manual_label_columns <- setdiff(
    required_manual_label_columns,
    names(manual_label_group_table)
  )
  
  if (length(missing_manual_label_columns) > 0) {
    stop(
      "Manual label group table is missing required columns:\n",
      paste(missing_manual_label_columns, collapse = "\n")
    )
  }
  
  if (any(is.na(manual_label_group_table$match_value)) ||
      any(manual_label_group_table$match_value == "")) {
    stop("Manual label group table contains blank or NA match_value entries.")
  }
  
  if (any(is.na(manual_label_group_table$display_label)) ||
      any(manual_label_group_table$display_label == "")) {
    stop("Manual label group table contains blank or NA display_label entries.")
  }
  
  if (any(is.na(manual_label_group_table$label_group)) ||
      any(manual_label_group_table$label_group == "")) {
    stop("Manual label group table contains blank or NA label_group entries.")
  }
  
  manual_groups_missing_fills <- setdiff(
    unique(manual_label_group_table$label_group),
    names(manual_label_group_fills)
  )
  
  if (length(manual_groups_missing_fills) > 0) {
    stop(
      "The following manual label groups are missing from manual_label_group_fills:\n",
      paste(manual_groups_missing_fills, collapse = "\n")
    )
  }
  
  manual_groups_missing_text_colors <- setdiff(
    unique(manual_label_group_table$label_group),
    names(manual_label_group_text_colors)
  )
  
  if (length(manual_groups_missing_text_colors) > 0) {
    stop(
      "The following manual label groups are missing from manual_label_group_text_colors:\n",
      paste(manual_groups_missing_text_colors, collapse = "\n")
    )
  }
}

# ============================================================================
# SECTION 2: PACKAGE LOADING
# ============================================================================
#
# Section 2 loads required R packages.
# No input files are detected, imported, standardized, or plotted in this section.
#
# Required packages:
#   ggplot2  = plotting
#   ggrepel  = non-overlapping plot labels
#   dplyr    = data manipulation
#   readr    = robust CSV/TSV import
#   stringr  = string cleaning and pattern detection
#   janitor  = column-name cleaning
#   tibble   = tidy data frame handling
#   svglite  = SVG export support
#
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 2: PACKAGE LOADING\n")
cat("--------------------------------------------------------\n")

# ----------------------------------------------------------------------------
# 2.1 Define required packages
# ----------------------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "dplyr",
  "readr",
  "stringr",
  "janitor",
  "tibble",
  "svglite"
)

# ----------------------------------------------------------------------------
# 2.2 Load required packages
# ----------------------------------------------------------------------------

load_required_packages <- function(packages) {
  
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      "The following required packages are not installed:\n",
      paste(missing_packages, collapse = "\n"),
      "\n\nPlease install them before running this script.\n",
      "Example:\n",
      "install.packages(c(",
      paste(shQuote(missing_packages), collapse = ", "),
      "))"
    )
  }
  
  for (pkg in packages) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
    cat("[LOADED]", pkg, "\n")
  }
  
  invisible(TRUE)
}

load_required_packages(required_packages)

# ----------------------------------------------------------------------------
# 2.3 Section 2 checkpoint
# ----------------------------------------------------------------------------

loaded_status <- vapply(
  required_packages,
  function(pkg) pkg %in% loadedNamespaces(),
  logical(1)
)

if (!all(loaded_status)) {
  stop(
    "SECTION 2 checkpoint failed. Packages not loaded:\n",
    paste(required_packages[!loaded_status], collapse = "\n")
  )
}

cat("[CHECKPOINT] Section 2 completed successfully.\n")
cat("[CHECKPOINT] All required packages are available and loaded.\n")

# ============================================================================
# SECTION 3: HELPER FUNCTIONS
# ============================================================================
#
# Section 3 defines reusable helper functions used throughout the script.
# No input files are detected, imported, standardized, or plotted in this section.
#
# Helper function groups:
#   3.1 String and filename helpers
#   3.2 Timestamp helper
#   3.3 Directory helper
#   3.4 Value-checking helpers
#   3.5 Section 3 checkpoint
#
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 3: HELPER FUNCTIONS\n")
cat("--------------------------------------------------------\n")

# ----------------------------------------------------------------------------
# 3.1 String and filename helpers
# ----------------------------------------------------------------------------

# Convert a character string into a safe filename component.
# This removes or replaces characters that commonly cause problems in filenames.

safe_name <- function(x) {
  
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) {
    return("unnamed")
  }
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[/\\\\:*?\"<>|]", "_")
  x <- stringr::str_replace_all(x, "\\s+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  
  if (!nzchar(x)) {
    x <- "unnamed"
  }
  
  return(x)
}

infer_comparison_name_from_filename <- function(input_file_name) {
  
  comparison <- tools::file_path_sans_ext(basename(input_file_name))
  
  # Remove common volcano-input prefixes.
  comparison <- stringr::str_remove(
    comparison,
    stringr::regex("^DESeq2_apeglm_Shrunken_Results_", ignore_case = TRUE)
  )
  comparison <- stringr::str_remove(
    comparison,
    stringr::regex("^DESeq2_Unshrunken_Results_", ignore_case = TRUE)
  )
  comparison <- stringr::str_remove(
    comparison,
    stringr::regex("^DESeq2_Full_Results_", ignore_case = TRUE)
  )
  comparison <- stringr::str_remove(
    comparison,
    stringr::regex("^DESeq2_", ignore_case = TRUE)
  )
  comparison <- stringr::str_remove(
    comparison,
    stringr::regex("^Tnseq_output_", ignore_case = TRUE)
  )
  
  # If filename includes project text before an RNA-seq token, keep the
  # comparison-like tail so titles do not repeat project/comparison text.
  if (stringr::str_detect(comparison, "_RNAseq_|_RNA_seq_|_RNA-seq_")) {
    comparison_parts <- stringr::str_split(
      comparison,
      "_RNAseq_|_RNA_seq_|_RNA-seq_",
      simplify = TRUE
    )
    
    if (ncol(comparison_parts) > 1) {
      comparison <- comparison_parts[1, ncol(comparison_parts)]
    }
  }
  
  comparison <- stringr::str_replace_all(comparison, "\\s+", "_")
  comparison <- stringr::str_replace_all(comparison, "_+", "_")
  comparison <- stringr::str_replace_all(comparison, "^_|_$", "")
  
  if (!nzchar(comparison)) {
    comparison <- "comparison"
  }
  
  return(comparison)
}

infer_plot_title_from_filename <- function(input_file_name) {
  
  title <- infer_comparison_name_from_filename(input_file_name)
  
  # Clean separators
  title <- stringr::str_replace_all(title, "_", " ")
  title <- stringr::str_squish(title)
  
  # Cosmetic replacements
  title <- stringr::str_replace_all(title, "Tet OFF", "Tet-OFF")
  title <- stringr::str_replace_all(title, "\\bd0\\b", "d0")
  title <- stringr::str_replace_all(title, "\\bd1\\b", "d1")
  title <- stringr::str_replace_all(title, "\\bd3\\b", "d3")
  title <- stringr::str_replace_all(title, "\\bd7\\b", "d7")
  title <- stringr::str_replace_all(title, "\\batc\\b", "+Atc")
  title <- stringr::str_replace_all(title, "\\bATC\\b", "+Atc")
  title <- stringr::str_replace_all(title, "\\brf\\b", "+RF")
  title <- stringr::str_replace_all(title, "\\bRF\\b", "+RF")
  
  # Remove accidental double plus signs
  title <- stringr::str_replace_all(title, "\\+\\+", "+")
  
  return(title)
}

italicize_leading_gene_in_title <- function(title_text, gene_to_italic) {
  
  if (is.null(title_text) || !is.character(title_text)) {
    return(title_text)
  }
  
  if (is.null(gene_to_italic) || !nzchar(gene_to_italic)) {
    return(title_text)
  }
  
  title_text <- stringr::str_squish(title_text)
  
  if (identical(title_text, gene_to_italic)) {
    return(
      as.expression(
        bquote(italic(.(gene_to_italic)))
      )
    )
  }
  
  gene_prefix <- paste0(gene_to_italic, " ")
  
  if (!stringr::str_starts(title_text, gene_prefix)) {
    return(title_text)
  }
  
  title_remainder <- stringr::str_remove(
    title_text,
    paste0("^", gene_to_italic, "\\s+")
  )
  
  title_remainder <- stringr::str_squish(title_remainder)
  
  as.expression(
    bquote(italic(.(gene_to_italic)) ~ .(title_remainder))
  )
}

# Collapse multiple name components into one safe tag.
# NULL or empty components are ignored.

make_file_tag <- function(...) {
  
  parts <- list(...)
  
  parts <- unlist(parts)
  parts <- parts[!is.na(parts)]
  parts <- parts[nzchar(parts)]
  
  if (length(parts) == 0) {
    return("volcano_plot")
  }
  
  tag <- paste(parts, collapse = "_")
  tag <- safe_name(tag)
  
  return(tag)
}


# ----------------------------------------------------------------------------
# 3.2 Timestamp helper
# ----------------------------------------------------------------------------

# Create a timestamp for output folders and logs.

make_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
}


# ----------------------------------------------------------------------------
# 3.3 Directory helper
# ----------------------------------------------------------------------------

# Create a directory if it does not already exist.

create_dir_if_missing <- function(dir_path) {
  
  if (is.null(dir_path) || length(dir_path) == 0 || is.na(dir_path) || !nzchar(dir_path)) {
    stop("Directory path is missing or empty.")
  }
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("[CREATED DIR]", dir_path, "\n")
  } else {
    cat("[DIR EXISTS]", dir_path, "\n")
  }
  
  invisible(dir_path)
}


# ----------------------------------------------------------------------------
# 3.4 Value-checking helpers
# ----------------------------------------------------------------------------

# Check that a value is TRUE or FALSE.

check_logical_single <- function(x, setting_name) {
  
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(setting_name, " must be a single TRUE or FALSE value.")
  }
  
  invisible(TRUE)
}


# Check that a numeric setting is a single non-missing number.

check_numeric_single <- function(x, setting_name, allow_zero = TRUE) {
  
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(setting_name, " must be a single numeric value.")
  }
  
  if (!allow_zero && x == 0) {
    stop(setting_name, " must be non-zero.")
  }
  
  invisible(TRUE)
}


# Check that a setting is one of a defined set of allowed values.

check_allowed_value <- function(x, setting_name, allowed_values) {
  
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(setting_name, " must be a single character value.")
  }
  
  if (!x %in% allowed_values) {
    stop(
      setting_name, " must be one of: ",
      paste(allowed_values, collapse = ", "),
      "\nCurrent value: ", x
    )
  }
  
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# 3.5 Section 3 checkpoint
# ----------------------------------------------------------------------------

required_section3_functions <- c(
  "safe_name",
  "infer_comparison_name_from_filename",
  "infer_plot_title_from_filename",
  "italicize_leading_gene_in_title",
  "make_file_tag",
  "make_timestamp",
  "create_dir_if_missing",
  "check_logical_single",
  "check_numeric_single",
  "check_allowed_value"
)

missing_section3_functions <- required_section3_functions[
  !vapply(required_section3_functions, exists, logical(1))
]

if (length(missing_section3_functions) > 0) {
  stop(
    "SECTION 3 checkpoint failed. Missing helper functions:\n",
    paste(missing_section3_functions, collapse = "\n")
  )
}

# Basic function tests.

test_safe_name <- safe_name("Treatment vs Control: DESeq2/results?")
test_comparison_name <- infer_comparison_name_from_filename(
  "DESeq2_Full_Results_Project_RNAseq_Treatment_vs_Control.csv"
)
test_plot_title <- infer_plot_title_from_filename(
  "DESeq2_Full_Results_Project_RNAseq_Treatment_vs_Control.csv"
)
test_file_tag <- make_file_tag("Volcano Plot", "Treatment vs Control")
test_timestamp <- make_timestamp()

if (!identical(test_safe_name, "Treatment_vs_Control_DESeq2_results")) {
  stop("SECTION 3 checkpoint failed. safe_name() did not return the expected result.")
}

if (!identical(test_comparison_name, "Treatment_vs_Control")) {
  stop("SECTION 3 checkpoint failed. infer_comparison_name_from_filename() did not return the expected result.")
}

if (!identical(test_plot_title, "Treatment vs Control")) {
  stop("SECTION 3 checkpoint failed. infer_plot_title_from_filename() did not return the expected result.")
}

if (!nzchar(test_file_tag)) {
  stop("SECTION 3 checkpoint failed. make_file_tag() returned an empty value.")
}

if (!nzchar(test_timestamp)) {
  stop("SECTION 3 checkpoint failed. make_timestamp() returned an empty value.")
}

cat("[CHECKPOINT] Section 3 completed successfully.\n")
cat("[CHECKPOINT] Helper functions are defined and passed basic tests.\n")

# ============================================================================
# SECTION 4: INPUT FILE DETECTION AND IMPORT
# ============================================================================
#
# Section 4 detects and imports the volcano input file.
#
# Expected workflow:
#   - The working directory contains exactly one DESeq2 or Tn-seq result file.
#   - If input_file is NULL, the script auto-detects the file.
#   - If input_file is manually provided, the script uses that file instead.
#
# This section imports the raw table and cleans column names.
# Input type standardization happens later in Section 5.
#
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 4: INPUT FILE DETECTION AND IMPORT\n")
cat("--------------------------------------------------------\n")

# ----------------------------------------------------------------------------
# 4.1 Detect or validate input file
# ----------------------------------------------------------------------------

if (use_auto_input_detection) {
  
  candidate_files <- list.files(
    path = input_dir,
    pattern = "\\.(csv|tsv|txt)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  volcano_candidates <- candidate_files[
    grepl(
      pattern = "DESeq2|Tnseq|TnSeq|Tn-seq|Tn_seq",
      x = basename(candidate_files),
      ignore.case = TRUE
    )
  ]
  
  cat("[SCAN] Tabular files found:", length(candidate_files), "\n")
  cat("[SCAN] Candidate volcano input files found:", length(volcano_candidates), "\n")
  
  if (length(volcano_candidates) == 0) {
    stop(
      "No candidate volcano input file was detected in the input directory.\n\n",
      "Input directory:\n",
      input_dir,
      "\n\nExpected one file with one of these filename patterns:\n",
      "  DESeq2\n",
      "  Tnseq\n",
      "  TnSeq\n",
      "  Tn-seq\n",
      "  Tn_seq\n\n",
      "Either place one matching result file in the working directory or set input_file manually in Section 1.1."
    )
  }
  
  if (length(volcano_candidates) > 1) {
    stop(
      "Multiple candidate volcano input files were detected:\n",
      paste(basename(volcano_candidates), collapse = "\n"),
      "\n\nPlease keep only one candidate result file in the working directory, ",
      "or set input_file manually in Section 1.1."
    )
  }
  
  input_file_path <- volcano_candidates[1]
  cat("[DETECTED] Input file:", basename(input_file_path), "\n")
  
} else {
  
  input_file_path <- file.path(input_dir, input_file)
  
  if (!file.exists(input_file_path)) {
    stop(
      "Manual input_file was specified but could not be found:\n",
      input_file_path
    )
  }
  
  cat("[MANUAL] Input file:", basename(input_file_path), "\n")
}

# ----------------------------------------------------------------------------
# 4.2 Determine file extension
# ----------------------------------------------------------------------------

input_extension <- tools::file_ext(input_file_path)
input_extension <- tolower(input_extension)

if (!input_extension %in% c("csv", "tsv", "txt")) {
  stop(
    "Unsupported input file extension: ", input_extension, "\n",
    "Supported extensions are: csv, tsv, txt"
  )
}

cat("[IMPORT] Input extension:", input_extension, "\n")

# ----------------------------------------------------------------------------
# 4.3 Import input table
# ----------------------------------------------------------------------------

if (input_extension == "csv") {
  
  results_raw <- readr::read_csv(
    file = input_file_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  
} else if (input_extension %in% c("tsv", "txt")) {
  
  results_raw <- readr::read_tsv(
    file = input_file_path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

cat("[IMPORT] Rows imported:", nrow(results_raw), "\n")
cat("[IMPORT] Columns imported:", ncol(results_raw), "\n")

# ----------------------------------------------------------------------------
# 4.4 Preserve original column names and clean names
# ----------------------------------------------------------------------------

original_column_names <- names(results_raw)

results_raw <- janitor::clean_names(results_raw)

cleaned_column_names <- names(results_raw)

column_name_map <- tibble::tibble(
  original_name = original_column_names,
  cleaned_name = cleaned_column_names
)

cat("[IMPORT] Column names cleaned with janitor::clean_names().\n")

# ----------------------------------------------------------------------------
# 4.5 Basic imported-data checks
# ----------------------------------------------------------------------------

if (nrow(results_raw) == 0) {
  stop("Imported input table has zero rows.")
}

if (ncol(results_raw) == 0) {
  stop("Imported input table has zero columns.")
}

if (anyDuplicated(names(results_raw)) > 0) {
  duplicated_names <- names(results_raw)[duplicated(names(results_raw))]
  stop(
    "Duplicate column names detected after cleaning:\n",
    paste(unique(duplicated_names), collapse = "\n")
  )
}

cat("[CHECK] Input table is non-empty.\n")
cat("[CHECK] No duplicate cleaned column names detected.\n")

# ----------------------------------------------------------------------------
# 4.6 Create preliminary input metadata
# ----------------------------------------------------------------------------

input_metadata <- list(
  input_file_path = input_file_path,
  input_file_name = basename(input_file_path),
  input_extension = input_extension,
  input_dir = input_dir,
  n_rows_raw = nrow(results_raw),
  n_cols_raw = ncol(results_raw),
  original_column_names = original_column_names,
  cleaned_column_names = cleaned_column_names
)

cat("[METADATA] Input metadata object created.\n")

if (is.null(comparison_name)) {
  comparison_name_final <- infer_comparison_name_from_filename(
    input_metadata$input_file_name
  )
} else {
  comparison_name_final <- comparison_name
}

cat("[METADATA] Comparison name:", comparison_name_final, "\n")

# ----------------------------------------------------------------------------
# 4.7 Section 4 checkpoint
# ----------------------------------------------------------------------------

required_section4_objects <- c(
  "input_file_path",
  "input_extension",
  "results_raw",
  "original_column_names",
  "cleaned_column_names",
  "column_name_map",
  "input_metadata",
  "comparison_name_final"
)

missing_section4_objects <- required_section4_objects[
  !vapply(required_section4_objects, exists, logical(1))
]

if (length(missing_section4_objects) > 0) {
  stop(
    "SECTION 4 checkpoint failed. Missing objects:\n",
    paste(missing_section4_objects, collapse = "\n")
  )
}

if (!file.exists(input_file_path)) {
  stop("SECTION 4 checkpoint failed. input_file_path does not exist.")
}

if (!inherits(results_raw, c("tbl_df", "tbl", "data.frame"))) {
  stop("SECTION 4 checkpoint failed. results_raw is not a data frame or tibble.")
}

cat("[CHECKPOINT] Section 4 completed successfully.\n")
cat("[CHECKPOINT] Input file detected/imported and column names cleaned.\n")

# ============================================================================
# SECTION 5: INPUT TYPE DETECTION AND STANDARDIZATION
# ============================================================================
#
# Section 5 detects the imported result-table type and converts it into a
# standardized volcano-ready data frame.
#
# Standardized output object:
#   volcano_df
#
# Standardized columns:
#   gene_id
#   gene_name
#   feature_label
#   log2FC
#   pvalue
#   padj
#   abundance
#   comparison
#   source_type
#
# Notes:
#   - DESeq2 standardization is implemented first.
#   - Tn-seq standardization will be added later after we define expected
#     Tn-seq output columns.
#
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 5: INPUT TYPE DETECTION AND STANDARDIZATION\n")
cat("--------------------------------------------------------\n")

# ----------------------------------------------------------------------------
# 5.1 Detect input type from cleaned column names
# ----------------------------------------------------------------------------

detect_input_type <- function(df) {
  
  cols <- names(df)
  
  is_deseq2 <- all(c("log2fold_change", "pvalue", "padj") %in% cols)
  
  if (is_deseq2) {
    return("DESeq2")
  }
  
  possible_tnseq_effect_cols <- c(
    "log2fc",
    "log2_fc",
    "log2_fold_change",
    "fitness",
    "log2_fitness",
    "effect_size"
  )
  
  possible_tnseq_sig_cols <- c(
    "pvalue",
    "p_value",
    "padj",
    "qvalue",
    "q_value",
    "fdr"
  )
  
  has_tnseq_effect <- any(possible_tnseq_effect_cols %in% cols)
  has_tnseq_significance <- any(possible_tnseq_sig_cols %in% cols)
  
  if (has_tnseq_effect && has_tnseq_significance) {
    return("TnSeq_candidate")
  }
  
  return("unknown")
}

detected_input_type <- detect_input_type(results_raw)

cat("[DETECTED] Input type:", detected_input_type, "\n")

if (!is.null(analysis_type)) {
  
  cat("[SETTINGS] Manual analysis_type was provided:", analysis_type, "\n")
  
  if (analysis_type %in% c("RNAseq", "DESeq2") && detected_input_type != "DESeq2") {
    warning(
      "Manual analysis_type suggests DESeq2/RNAseq, but detected input type is: ",
      detected_input_type
    )
  }
}

# ----------------------------------------------------------------------------
# 5.2 Standardize DESeq2 input
# ----------------------------------------------------------------------------

standardize_deseq2 <- function(df, comparison_name = NULL) {
  
  required_cols <- c(
    "gene_id",
    "log2fold_change",
    "pvalue",
    "padj"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "DESeq2 standardization failed. Missing required columns:\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  if ("gene" %in% names(df)) {
    gene_name <- df$gene
  } else if ("gene_name" %in% names(df)) {
    gene_name <- df$gene_name
  } else {
    gene_name <- rep(NA_character_, nrow(df))
  }
  
  if ("base_mean" %in% names(df)) {
    abundance <- df$base_mean
  } else {
    abundance <- rep(NA_real_, nrow(df))
  }
  
  # Optional display_label column.
  # If present in the input CSV, this is used first for plot labels.
  # This is useful for curated labels such as ACP3PV IDs corrected to Rv numbers.
  
  display_label <- if ("display_label" %in% names(df)) {
    as.character(df$display_label)
  } else {
    rep(NA_character_, nrow(df))
  }
  
  # feature_label controls what is displayed on the plot.
  # Priority:
  #   1. display_label from input CSV, if provided
  #   2. gene_name, if prefer_gene_name_labels is TRUE
  #   3. gene_id fallback
  
  if (prefer_gene_name_labels) {
    
    feature_label <- dplyr::case_when(
      !is.na(display_label) & display_label != "" ~ display_label,
      !is.na(gene_name) & gene_name != "" ~ gene_name,
      !is.na(df$gene_id) & df$gene_id != "" ~ df$gene_id,
      TRUE ~ NA_character_
    )
    
  } else {
    
    feature_label <- dplyr::case_when(
      !is.na(display_label) & display_label != "" ~ display_label,
      !is.na(df$gene_id) & df$gene_id != "" ~ df$gene_id,
      !is.na(gene_name) & gene_name != "" ~ gene_name,
      TRUE ~ NA_character_
    )
  }
  
  standardized_df <- tibble::tibble(
    gene_id = as.character(df$gene_id),
    gene_name = as.character(gene_name),
    feature_label = as.character(feature_label),
    log2FC = as.numeric(df$log2fold_change),
    pvalue = as.numeric(df$pvalue),
    padj = as.numeric(df$padj),
    abundance = as.numeric(abundance),
    comparison = ifelse(is.null(comparison_name), NA_character_, comparison_name),
    source_type = "DESeq2"
  )
  
  return(standardized_df)
}

# ----------------------------------------------------------------------------
# 5.3 Placeholder for Tn-seq standardization
# ----------------------------------------------------------------------------

standardize_tnseq <- function(df, comparison_name = NULL) {
  
  stop(
    "Tn-seq standardization has not been fully defined yet.\n",
    "Before enabling this branch, define the expected Tn-seq columns for:\n",
    "  - feature/gene identifier\n",
    "  - feature/gene label\n",
    "  - log2-scale effect size or fitness score\n",
    "  - p-value and/or adjusted p-value/q-value/FDR\n"
  )
}

# ----------------------------------------------------------------------------
# 5.4 Apply standardization
# ----------------------------------------------------------------------------

if (detected_input_type == "DESeq2") {
  
  volcano_df <- standardize_deseq2(
    df = results_raw,
    comparison_name = comparison_name_final
  )
  
} else if (detected_input_type == "TnSeq_candidate") {
  
  volcano_df <- standardize_tnseq(
    df = results_raw,
    comparison_name = comparison_name_final
  )
  
} else {
  
  stop(
    "Input type could not be detected or is not yet supported.\n\n",
    "Detected input type: ", detected_input_type, "\n\n",
    "Cleaned input columns:\n",
    paste(names(results_raw), collapse = "\n")
  )
}

cat("[STANDARDIZED] Rows:", nrow(volcano_df), "\n")
cat("[STANDARDIZED] Columns:", ncol(volcano_df), "\n")

# ----------------------------------------------------------------------------
# 5.5 Standardized-data checks
# ----------------------------------------------------------------------------

required_volcano_cols <- c(
  "gene_id",
  "gene_name",
  "feature_label",
  "log2FC",
  "pvalue",
  "padj",
  "abundance",
  "comparison",
  "source_type"
)

missing_volcano_cols <- setdiff(required_volcano_cols, names(volcano_df))

if (length(missing_volcano_cols) > 0) {
  stop(
    "Standardized volcano table is missing required columns:\n",
    paste(missing_volcano_cols, collapse = "\n")
  )
}

if (nrow(volcano_df) == 0) {
  stop("Standardized volcano table has zero rows.")
}

if (all(is.na(volcano_df$log2FC))) {
  stop("All standardized log2FC values are NA.")
}

if (all(is.na(volcano_df$pvalue)) && all(is.na(volcano_df$padj))) {
  stop("Both standardized pvalue and padj columns are entirely NA.")
}

cat("[CHECK] Required standardized columns are present.\n")
cat("[CHECK] Standardized volcano table is non-empty.\n")
cat("[CHECK] log2FC contains usable values.\n")
cat("[CHECK] pvalue/padj contains usable values.\n")

# ----------------------------------------------------------------------------
# 5.6 Report standardization summary
# ----------------------------------------------------------------------------

cat("[REPORT] Source type:", unique(volcano_df$source_type), "\n")
cat("[REPORT] Rows with gene_id:", sum(!is.na(volcano_df$gene_id) & volcano_df$gene_id != ""), "\n")
cat("[REPORT] Rows with gene_name:", sum(!is.na(volcano_df$gene_name) & volcano_df$gene_name != ""), "\n")
cat("[REPORT] Rows with feature_label:", sum(!is.na(volcano_df$feature_label) & volcano_df$feature_label != ""), "\n")
cat("[REPORT] Rows with NA log2FC:", sum(is.na(volcano_df$log2FC)), "\n")
cat("[REPORT] Rows with NA pvalue:", sum(is.na(volcano_df$pvalue)), "\n")
cat("[REPORT] Rows with NA padj:", sum(is.na(volcano_df$padj)), "\n")

# ----------------------------------------------------------------------------
# 5.7 Section 5 checkpoint
# ----------------------------------------------------------------------------

required_section5_objects <- c(
  "detected_input_type",
  "comparison_name_final",
  "volcano_df",
  "required_volcano_cols"
)

missing_section5_objects <- required_section5_objects[
  !vapply(required_section5_objects, exists, logical(1))
]

if (length(missing_section5_objects) > 0) {
  stop(
    "SECTION 5 checkpoint failed. Missing objects:\n",
    paste(missing_section5_objects, collapse = "\n")
  )
}

cat("[CHECKPOINT] Section 5 completed successfully.\n")
cat("[CHECKPOINT] Input table standardized into volcano_df.\n")

# ============================================================================
# SECTION 6: PLOTTING VALUE PREPARATION AND QUALITY CHECKS
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 6: PLOTTING VALUE PREPARATION AND QUALITY CHECKS\n")
cat("--------------------------------------------------------\n")

available_significance_cols <- c()

if ("padj" %in% names(volcano_df) && any(!is.na(volcano_df$padj))) {
  available_significance_cols <- c(available_significance_cols, "padj")
}

if ("pvalue" %in% names(volcano_df) && any(!is.na(volcano_df$pvalue))) {
  available_significance_cols <- c(available_significance_cols, "pvalue")
}

if (length(available_significance_cols) == 0) {
  stop("No usable significance column found. Both padj and pvalue are missing or entirely NA.")
}

if (significance_column_preference == "auto") {
  
  if ("padj" %in% available_significance_cols) {
    significance_col <- "padj"
    significance_cutoff <- padj_cutoff
    y_axis_label_final <- expression(-log[10]("adjusted p-value"))
  } else {
    significance_col <- "pvalue"
    significance_cutoff <- pvalue_cutoff
    y_axis_label_final <- expression(-log[10]("p-value"))
  }
  
} else if (significance_column_preference == "padj") {
  
  if (!"padj" %in% available_significance_cols) {
    stop("significance_column_preference is 'padj', but no usable padj column is available.")
  }
  
  significance_col <- "padj"
  significance_cutoff <- padj_cutoff
  y_axis_label_final <- expression(-log[10]("adjusted p-value"))
  
} else if (significance_column_preference == "pvalue") {
  
  if (!"pvalue" %in% available_significance_cols) {
    stop("significance_column_preference is 'pvalue', but no usable pvalue column is available.")
  }
  
  significance_col <- "pvalue"
  significance_cutoff <- pvalue_cutoff
  y_axis_label_final <- expression(-log[10]("p-value"))
  
} else {
  
  stop(
    "Invalid significance_column_preference: ", significance_column_preference,
    "\nAllowed values are: auto, padj, pvalue"
  )
}

cat("[PLOT PREP] Significance column selected:", significance_col, "\n")
cat("[PLOT PREP] Significance cutoff selected:", significance_cutoff, "\n")

volcano_plot_df <- volcano_df %>%
  dplyr::mutate(
    significance_value_raw = .data[[significance_col]]
  )

cat("[PLOT PREP] Plotting data frame created.\n")

n_na_significance <- sum(is.na(volcano_plot_df$significance_value_raw))
n_total_before_filter <- nrow(volcano_plot_df)

cat("[PLOT PREP] Rows before significance filtering:", n_total_before_filter, "\n")
cat("[PLOT PREP] Rows with NA selected significance value:", n_na_significance, "\n")

volcano_plot_df <- volcano_plot_df %>%
  dplyr::filter(!is.na(significance_value_raw))

n_total_after_filter <- nrow(volcano_plot_df)

cat("[PLOT PREP] Rows retained for plotting:", n_total_after_filter, "\n")
cat("[PLOT PREP] Rows removed from plotting due to NA significance:", 
    n_total_before_filter - n_total_after_filter, "\n")

n_nonpositive_significance <- sum(
  volcano_plot_df$significance_value_raw <= 0,
  na.rm = TRUE
)

if (n_nonpositive_significance > 0) {
  
  smallest_positive_significance <- min(
    volcano_plot_df$significance_value_raw[
      volcano_plot_df$significance_value_raw > 0
    ],
    na.rm = TRUE
  )
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::mutate(
      significance_value_for_plot = dplyr::if_else(
        significance_value_raw <= 0,
        smallest_positive_significance,
        significance_value_raw
      )
    )
  
  cat("[PLOT PREP] Non-positive significance values detected:", 
      n_nonpositive_significance, "\n")
  cat("[PLOT PREP] Replaced non-positive values with smallest positive value:", 
      smallest_positive_significance, "\n")
  
} else {
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::mutate(
      significance_value_for_plot = significance_value_raw
    )
  
  cat("[PLOT PREP] No non-positive significance values detected.\n")
}

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    neg_log10_significance = -log10(significance_value_for_plot)
  )

if (any(!is.finite(volcano_plot_df$neg_log10_significance))) {
  stop("Non-finite -log10 significance values detected after plot preparation.")
}

cat("[PLOT PREP] -log10 significance values calculated.\n")

n_missing_log2fc <- sum(is.na(volcano_plot_df$log2FC))

if (n_missing_log2fc > 0) {
  
  cat("[PLOT PREP] Rows with NA log2FC removed from plotting:", 
      n_missing_log2fc, "\n")
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::filter(!is.na(log2FC))
  
} else {
  
  cat("[PLOT PREP] No missing log2FC values detected.\n")
}

cat("[REPORT] log2FC range:", 
    round(min(volcano_plot_df$log2FC, na.rm = TRUE), 3), "to",
    round(max(volcano_plot_df$log2FC, na.rm = TRUE), 3), "\n")

cat("[REPORT] -log10 significance range:", 
    round(min(volcano_plot_df$neg_log10_significance, na.rm = TRUE), 3), "to",
    round(max(volcano_plot_df$neg_log10_significance, na.rm = TRUE), 3), "\n")

# ----------------------------------------------------------------------------
# 6.1 Optional display-only y-axis cap
# ----------------------------------------------------------------------------

if (!is.null(y_axis_limit) && y_axis_cap_display_values) {
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::mutate(
      neg_log10_significance_for_plot = pmin(
        neg_log10_significance,
        y_axis_limit
      ),
      y_value_capped_for_plot = neg_log10_significance > y_axis_limit
    )
  
  cat("[PLOT PREP] Display-only y-axis cap applied at:", y_axis_limit, "\n")
  cat("[PLOT PREP] Points displayed at y-axis cap:", 
      sum(volcano_plot_df$y_value_capped_for_plot, na.rm = TRUE), "\n")
  
} else {
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::mutate(
      neg_log10_significance_for_plot = neg_log10_significance,
      y_value_capped_for_plot = FALSE
    )
  
  cat("[PLOT PREP] No display-only y-axis cap applied.\n")
}

required_section6_objects <- c(
  "significance_col",
  "significance_cutoff",
  "y_axis_label_final",
  "volcano_plot_df"
)

missing_section6_objects <- required_section6_objects[
  !vapply(required_section6_objects, exists, logical(1))
]

if (length(missing_section6_objects) > 0) {
  stop(
    "SECTION 6 checkpoint failed. Missing objects:\n",
    paste(missing_section6_objects, collapse = "\n")
  )
}

required_plot_cols <- c(
  "significance_value_raw",
  "significance_value_for_plot",
  "neg_log10_significance",
  "neg_log10_significance_for_plot",
  "y_value_capped_for_plot"
)

missing_plot_cols <- setdiff(required_plot_cols, names(volcano_plot_df))

if (length(missing_plot_cols) > 0) {
  stop(
    "SECTION 6 checkpoint failed. Missing plot-preparation columns:\n",
    paste(missing_plot_cols, collapse = "\n")
  )
}

if (nrow(volcano_plot_df) == 0) {
  stop("SECTION 6 checkpoint failed. volcano_plot_df has zero rows.")
}



cat("[CHECKPOINT] Section 6 completed successfully.\n")
cat("[CHECKPOINT] Plotting values prepared in volcano_plot_df.\n")

# ============================================================================
# SECTION 7: SIGNIFICANCE CLASSIFICATION
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 7: SIGNIFICANCE CLASSIFICATION\n")
cat("--------------------------------------------------------\n")

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    significance_category = dplyr::case_when(
      significance_value_raw <= significance_cutoff &
        log2FC >= log2fc_cutoff ~ "Up",
      
      significance_value_raw <= significance_cutoff &
        log2FC <= -log2fc_cutoff ~ "Down",
      
      classify_significant_low_effect &
        significance_value_raw <= significance_cutoff &
        abs(log2FC) < log2fc_cutoff ~ "Significant_low_effect",
      
      TRUE ~ "Not_significant"
    )
  )

cat("[CLASSIFY] Significance categories assigned.\n")

significance_category_levels <- c(
  "Up",
  "Down",
  "Significant_low_effect",
  "Not_significant"
)

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    significance_category = factor(
      significance_category,
      levels = significance_category_levels
    )
  )

cat("[CLASSIFY] Significance category order set.\n")

significance_category_counts <- volcano_plot_df %>%
  dplyr::count(significance_category, name = "n") %>%
  dplyr::arrange(factor(significance_category, levels = significance_category_levels))

cat("[REPORT] Significance category counts:\n")
print(significance_category_counts)

significant_hits_df <- volcano_plot_df %>%
  dplyr::filter(
    significance_category %in% c("Up", "Down", "Significant_low_effect")
  ) %>%
  dplyr::arrange(
    significance_value_raw,
    dplyr::desc(abs(log2FC))
  )

cat("[REPORT] Significant hits:", nrow(significant_hits_df), "\n")
cat("[REPORT] Up hits:", sum(volcano_plot_df$significance_category == "Up", na.rm = TRUE), "\n")
cat("[REPORT] Down hits:", sum(volcano_plot_df$significance_category == "Down", na.rm = TRUE), "\n")
cat("[REPORT] Significant low-effect hits:", 
    sum(volcano_plot_df$significance_category == "Significant_low_effect", na.rm = TRUE), "\n")

observed_categories <- as.character(unique(volcano_plot_df$significance_category))
observed_categories <- observed_categories[!is.na(observed_categories)]

missing_color_categories <- setdiff(observed_categories, names(volcano_colors))

if (length(missing_color_categories) > 0) {
  stop(
    "Some significance categories do not have colors defined in volcano_colors:\n",
    paste(missing_color_categories, collapse = "\n")
  )
}

cat("[CHECK] All observed significance categories have assigned colors.\n")

required_section7_objects <- c(
  "significance_category_levels",
  "significance_category_counts",
  "significant_hits_df"
)

missing_section7_objects <- required_section7_objects[
  !vapply(required_section7_objects, exists, logical(1))
]

if (length(missing_section7_objects) > 0) {
  stop(
    "SECTION 7 checkpoint failed. Missing objects:\n",
    paste(missing_section7_objects, collapse = "\n")
  )
}

if (!"significance_category" %in% names(volcano_plot_df)) {
  stop("SECTION 7 checkpoint failed. significance_category column was not created.")
}

if (any(is.na(volcano_plot_df$significance_category))) {
  stop("SECTION 7 checkpoint failed. Some rows have NA significance_category.")
}

cat("[CHECKPOINT] Section 7 completed successfully.\n")
cat("[CHECKPOINT] Volcano plot significance categories assigned.\n")

# ============================================================================
# SECTION 8: AUTOMATIC AND MANUAL GROUPED LABEL SELECTION
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 8: AUTOMATIC AND MANUAL GROUPED LABEL SELECTION\n")
cat("--------------------------------------------------------\n")

# Initialize label-tracking columns.

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    plot_label = "",
    selected_auto_label = FALSE,
    selected_manual_label = FALSE,
    manual_label_group = NA_character_,
    manual_label_display = NA_character_,
    label_category = NA_character_
  )

cat("[LABELS] Initialized label columns.\n")

# ----------------------------------------------------------------------------
# 8.1 Label ranking helper
# ----------------------------------------------------------------------------

rank_label_candidates <- function(df, ranking_method, direction) {
  
  if (!direction %in% c("Up", "Down")) {
    stop("direction must be either 'Up' or 'Down'.")
  }
  
  if (ranking_method == "padj_then_log2FC") {
    
    ranked_df <- df %>%
      dplyr::arrange(
        significance_value_raw,
        dplyr::desc(abs(log2FC))
      )
    
  } else if (ranking_method == "log2FC_then_padj") {
    
    if (direction == "Up") {
      
      ranked_df <- df %>%
        dplyr::arrange(
          dplyr::desc(log2FC),
          significance_value_raw
        )
      
    } else if (direction == "Down") {
      
      ranked_df <- df %>%
        dplyr::arrange(
          log2FC,
          significance_value_raw
        )
    }
    
  } else if (ranking_method == "combined_score") {
    
    ranked_df <- df %>%
      dplyr::mutate(
        label_score = abs(log2FC) * neg_log10_significance
      ) %>%
      dplyr::arrange(
        dplyr::desc(label_score)
      )
    
  } else {
    
    stop(
      "Invalid label_ranking_method: ", ranking_method,
      "\nAllowed values are: padj_then_log2FC, log2FC_then_padj, combined_score"
    )
  }
  
  return(ranked_df)
}

cat("[LABELS] Label ranking helper defined.\n")

# ----------------------------------------------------------------------------
# 8.2 Select automatic top Up labels
# ----------------------------------------------------------------------------

if (label_top_hits && top_n_up > 0) {
  
  up_label_candidates <- volcano_plot_df %>%
    dplyr::filter(significance_category == "Up")
  
  if (label_only_significant) {
    up_label_candidates <- up_label_candidates %>%
      dplyr::filter(
        significance_value_raw <= significance_cutoff,
        abs(log2FC) >= log2fc_cutoff
      )
  }
  
  up_labels_df <- up_label_candidates %>%
    rank_label_candidates(
      ranking_method = label_ranking_method,
      direction = "Up"
    ) %>%
    dplyr::slice_head(n = top_n_up)
  
} else {
  
  up_labels_df <- volcano_plot_df[0, ]
}

cat("[LABELS] Automatic Up labels selected:", nrow(up_labels_df), "\n")

# ----------------------------------------------------------------------------
# 8.3 Select automatic top Down labels
# ----------------------------------------------------------------------------

if (label_top_hits && top_n_down > 0) {
  
  down_label_candidates <- volcano_plot_df %>%
    dplyr::filter(significance_category == "Down")
  
  if (label_only_significant) {
    down_label_candidates <- down_label_candidates %>%
      dplyr::filter(
        significance_value_raw <= significance_cutoff,
        abs(log2FC) >= log2fc_cutoff
      )
  }
  
  down_labels_df <- down_label_candidates %>%
    rank_label_candidates(
      ranking_method = label_ranking_method,
      direction = "Down"
    ) %>%
    dplyr::slice_head(n = top_n_down)
  
} else {
  
  down_labels_df <- volcano_plot_df[0, ]
}

cat("[LABELS] Automatic Down labels selected:", nrow(down_labels_df), "\n")

# ----------------------------------------------------------------------------
# 8.4 Apply automatic top DEG labels
# ----------------------------------------------------------------------------

auto_labels_df <- dplyr::bind_rows(
  up_labels_df,
  down_labels_df
) %>%
  dplyr::distinct(gene_id, feature_label, .keep_all = TRUE)

cat("[LABELS] Total automatic labels selected:", nrow(auto_labels_df), "\n")

if (nrow(auto_labels_df) > 0) {
  
  auto_label_lookup <- auto_labels_df %>%
    dplyr::select(gene_id, feature_label) %>%
    dplyr::mutate(selected_auto_label_lookup = TRUE)
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::left_join(
      auto_label_lookup,
      by = c("gene_id", "feature_label")
    ) %>%
    dplyr::mutate(
      selected_auto_label = !is.na(selected_auto_label_lookup),
      plot_label = dplyr::if_else(
        selected_auto_label,
        feature_label,
        plot_label
      ),
      label_category = dplyr::if_else(
        selected_auto_label,
        "Top DEG",
        label_category
      )
    ) %>%
    dplyr::select(-selected_auto_label_lookup)
}

cat("[LABELS] Automatic labels applied to volcano_plot_df.\n")

# ----------------------------------------------------------------------------
# 8.5 Select and apply manual grouped labels
# ----------------------------------------------------------------------------

manual_labels_df <- volcano_plot_df[0, ]

if (manual_labels_enabled && nrow(manual_label_group_table) > 0) {
  
  manual_match_columns <- intersect(
    c("gene_id", "gene", "gene_name", "gene_label", "feature_label"),
    names(volcano_plot_df)
  )
  
  if (length(manual_match_columns) == 0) {
    stop(
      "Manual labels cannot be matched because none of these columns exist:\n",
      "gene_id, gene, gene_name, gene_label, feature_label"
    )
  }
  
  cat(
    "[LABELS] Manual labels will be matched against:",
    paste(manual_match_columns, collapse = ", "),
    "\n"
  )
  
  for (i in seq_len(nrow(manual_label_group_table))) {
    
    current_match_value <- manual_label_group_table$match_value[i]
    current_display_label <- manual_label_group_table$display_label[i]
    current_label_group <- manual_label_group_table$label_group[i]
    
    current_matches <- rep(FALSE, nrow(volcano_plot_df))
    
    for (current_column in manual_match_columns) {
      current_matches <- current_matches |
        (
          !is.na(volcano_plot_df[[current_column]]) &
            volcano_plot_df[[current_column]] == current_match_value
        )
    }
    
    if (!any(current_matches)) {
      warning(
        "Manual label match was not found: ",
        current_match_value
      )
    }
    
    already_assigned <- current_matches & volcano_plot_df$selected_manual_label
    
    if (any(already_assigned)) {
      warning(
        "Manual label match was assigned to more than one manual group. ",
        "Keeping the first assignment for: ",
        current_match_value
      )
    }
    
    rows_to_assign <- current_matches & !volcano_plot_df$selected_manual_label
    
    volcano_plot_df$selected_manual_label[rows_to_assign] <- TRUE
    volcano_plot_df$manual_label_group[rows_to_assign] <- current_label_group
    volcano_plot_df$manual_label_display[rows_to_assign] <- current_display_label
  }
  
  if (!manual_labels_even_if_not_significant) {
    
    volcano_plot_df <- volcano_plot_df %>%
      dplyr::mutate(
        selected_manual_label = selected_manual_label &
          significance_value_raw <= significance_cutoff &
          abs(log2FC) >= log2fc_cutoff
      )
  }
  
  manual_labels_df <- volcano_plot_df %>%
    dplyr::filter(selected_manual_label)
  
} else {
  
  cat("[LABELS] Manual grouped labels are disabled or manual label table is empty.\n")
}

cat("[LABELS] Manual grouped labels selected:", nrow(manual_labels_df), "\n")

# ----------------------------------------------------------------------------
# 8.6 Combine automatic and manual label decisions
# ----------------------------------------------------------------------------

# Manual labels override plotted label text for manually selected genes.
# Automatic top DEG labels remain unchanged for automatic-only genes.

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    plot_label = dplyr::case_when(
      selected_manual_label ~ manual_label_display,
      selected_auto_label ~ plot_label,
      TRUE ~ ""
    )
  )

# If manual labels should replace automatic labels, remove automatic-only labels.

if (manual_labels_enabled && !manual_labels_add_to_auto_labels) {
  
  volcano_plot_df <- volcano_plot_df %>%
    dplyr::mutate(
      plot_label = dplyr::if_else(
        selected_auto_label & !selected_manual_label,
        "",
        plot_label
      ),
      selected_auto_label = dplyr::if_else(
        selected_auto_label & !selected_manual_label,
        FALSE,
        selected_auto_label
      )
    )
}

# Assign final label category.
#
# Top DEG only:
#   "Top DEG"
#
# Manual only:
#   manual label group, such as "KEGG ORA down"
#
# Both automatic top DEG and manual enrichment-connected gene:
#   "Top DEG + KEGG ORA down"

volcano_plot_df <- volcano_plot_df %>%
  dplyr::mutate(
    label_category = dplyr::case_when(
      selected_auto_label & selected_manual_label ~ paste0(
        "Top DEG + ",
        manual_label_group
      ),
      selected_manual_label ~ manual_label_group,
      selected_auto_label ~ "Top DEG",
      TRUE ~ NA_character_
    )
  )

# Apply optional label exclusion pattern to final plotted labels.
# This happens after manual display labels are assigned, so a locus tag manually
# displayed as a gene name will not be excluded unless the displayed label itself
# matches the exclusion pattern.

if (!is.null(exclude_label_pattern)) {
  
  excluded_label_rows <- stringr::str_detect(
    volcano_plot_df$plot_label,
    stringr::regex(exclude_label_pattern)
  )
  
  excluded_label_rows[is.na(excluded_label_rows)] <- FALSE
  
  volcano_plot_df$plot_label[excluded_label_rows] <- ""
  volcano_plot_df$label_category[excluded_label_rows] <- NA_character_
  volcano_plot_df$selected_auto_label[excluded_label_rows] <- FALSE
  volcano_plot_df$selected_manual_label[excluded_label_rows] <- FALSE
  
  cat(
    "[LABELS] Labels excluded by pattern:",
    sum(excluded_label_rows),
    "\n"
  )
}

# Final labeled data frame used for reporting, export, and plotting.

labeled_features_df <- volcano_plot_df %>%
  dplyr::filter(
    !is.na(plot_label),
    plot_label != ""
  ) %>%
  dplyr::distinct(gene_id, feature_label, .keep_all = TRUE)

labeled_data <- labeled_features_df

cat("[LABELS] Total unique labels selected:", nrow(labeled_features_df), "\n")

# ----------------------------------------------------------------------------
# 8.7 Prepare manual label color groups
# ----------------------------------------------------------------------------

# These extended color vectors include both plain manual groups and combined
# groups such as "Top DEG + KEGG ORA down".

manual_label_group_fills_extended <- manual_label_group_fills
manual_label_group_text_colors_extended <- manual_label_group_text_colors

for (current_group in names(manual_label_group_fills)) {
  
  combined_group <- paste0("Top DEG + ", current_group)
  
  manual_label_group_fills_extended[combined_group] <-
    manual_label_group_fills[current_group]
  
  manual_label_group_text_colors_extended[combined_group] <-
    manual_label_group_text_colors[current_group]
}

manual_label_group_fills_extended["Top DEG"] <- "white"
manual_label_group_text_colors_extended["Top DEG"] <- label_text_color

cat("[LABELS] Manual label color groups prepared.\n")

# ----------------------------------------------------------------------------
# 8.8 Label report
# ----------------------------------------------------------------------------

cat("[REPORT] Selected labels:\n")

selected_label_report <- volcano_plot_df %>%
  dplyr::filter(plot_label != "") %>%
  dplyr::select(
    gene_id,
    gene_name,
    feature_label,
    plot_label,
    selected_auto_label,
    selected_manual_label,
    manual_label_group,
    manual_label_display,
    label_category,
    log2FC,
    significance_value_raw,
    neg_log10_significance,
    significance_category
  ) %>%
  dplyr::mutate(
    label_sort_value = dplyr::case_when(
      significance_category == "Up" ~ -log2FC,
      significance_category == "Down" ~ log2FC,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::arrange(
    label_category,
    significance_category,
    label_sort_value,
    significance_value_raw
  ) %>%
  dplyr::select(-label_sort_value)

print(selected_label_report)

cat("[REPORT] Label categories:\n")
print(table(selected_label_report$label_category, useNA = "ifany"))

# ----------------------------------------------------------------------------
# 8.9 Section 8 checkpoint
# ----------------------------------------------------------------------------

required_section8_objects <- c(
  "up_labels_df",
  "down_labels_df",
  "auto_labels_df",
  "manual_labels_df",
  "labeled_features_df",
  "labeled_data",
  "selected_label_report",
  "manual_label_group_fills_extended",
  "manual_label_group_text_colors_extended"
)

missing_section8_objects <- required_section8_objects[
  !vapply(required_section8_objects, exists, logical(1))
]

if (length(missing_section8_objects) > 0) {
  stop(
    "SECTION 8 checkpoint failed. Missing objects:\n",
    paste(missing_section8_objects, collapse = "\n")
  )
}

required_section8_columns <- c(
  "plot_label",
  "selected_auto_label",
  "selected_manual_label",
  "manual_label_group",
  "manual_label_display",
  "label_category"
)

missing_section8_columns <- setdiff(
  required_section8_columns,
  names(volcano_plot_df)
)

if (length(missing_section8_columns) > 0) {
  stop(
    "SECTION 8 checkpoint failed. Missing volcano_plot_df columns:\n",
    paste(missing_section8_columns, collapse = "\n")
  )
}

if (any(is.na(volcano_plot_df$plot_label))) {
  stop("SECTION 8 checkpoint failed. Some plot_label values are NA.")
}

actual_label_count <- sum(volcano_plot_df$plot_label != "")

cat("[LABELS] Total labels assigned:", actual_label_count, "\n")

cat("[CHECKPOINT] Section 8 completed successfully.\n")
cat("[CHECKPOINT] Volcano labels selected and applied.\n")

# ============================================================================
# SECTION 9: VOLCANO PLOT GENERATION
# ============================================================================
#
# Section 9 creates the volcano plot as a ggplot object.
#
# Input object:
#   volcano_plot_df
#
# Output object:
#   volcano_plot
#
# Notes:
#   - Plot export happens later in Section 10.
#   - This section prints the plot for visual inspection.
#
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 9: VOLCANO PLOT GENERATION\n")
cat("--------------------------------------------------------\n")

# ----------------------------------------------------------------------------
# 9.1 Build plot title and subtitle
# ----------------------------------------------------------------------------

if (is.null(plot_title)) {
  
  if (clean_inferred_plot_title) {
    plot_title_final <- infer_plot_title_from_filename(input_metadata$input_file_name)
  } else {
    plot_title_final <- tools::file_path_sans_ext(input_metadata$input_file_name)
  }
  
} else {
  
  plot_title_final <- plot_title
}
if (
  is.character(plot_title_final) &&
  auto_italicize_title_gene
) {
  
  plot_title_final <- italicize_leading_gene_in_title(
    title_text = plot_title_final,
    gene_to_italic = title_gene_to_italic
  )
}

if (is.null(plot_subtitle)) {
  
  plot_subtitle_final <- paste0(
    "|log2FC| >= ", log2fc_cutoff,
    ", ",
    significance_col,
    " <= ", significance_cutoff
  )
  
} else {
  
  plot_subtitle_final <- plot_subtitle
}

y_axis_label_final <- y_axis_label

cat(
  "[PLOT] Plot title:",
  paste(deparse(plot_title_final), collapse = ""),
  "\n"
)
cat("[PLOT] Plot subtitle:", plot_subtitle_final, "\n")

# ----------------------------------------------------------------------------
# 9.2 Build final plot caption
# ----------------------------------------------------------------------------

caption_parts <- c()

if (show_input_file_in_caption) {
  caption_parts <- c(
    caption_parts,
    paste0("Input: ", input_metadata$input_file_name)
  )
}

if (show_rows_plotted_in_caption) {
  caption_parts <- c(
    caption_parts,
    paste0("Rows plotted: ", nrow(volcano_plot_df))
  )
}

if (show_color_explanation_in_caption) {
  caption_parts <- c(
    caption_parts,
    caption_color_explanation
  )
}

# Add manual label group key to the plot caption.
# This keeps the manual-label explanation outside the plot axes.

if (show_manual_label_key_in_caption && exists("selected_label_report")) {
  
  manual_label_groups_for_caption <- selected_label_report %>%
    dplyr::filter(
      selected_manual_label,
      !is.na(manual_label_group)
    ) %>%
    dplyr::pull(manual_label_group) %>%
    unique()
  
  if (length(manual_label_groups_for_caption) > 0) {
    
    manual_label_group_key <- paste0(
      manual_label_groups_for_caption,
      " = ",
      manual_label_group_color_names[manual_label_groups_for_caption]
    )
    
    caption_parts <- c(
      caption_parts,
      paste0(
        "Manual colored labels: ",
        paste(manual_label_group_key, collapse = "; ")
      )
    )
  }
}

plot_caption_final <- paste(caption_parts, collapse = "\n")

cat("[PLOT] Plot caption:", plot_caption_final, "\n")

# ----------------------------------------------------------------------------
# 9.3 Define volcano plotting function
# ----------------------------------------------------------------------------

make_volcano_plot <- function(df) {
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log2FC,
      y = neg_log10_significance_for_plot,
      color = significance_category
    )
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = ifelse(df$significance_category %in% c("Up", "Down"),
                    point_size, background_point_size)
    ) +
    ggplot2::scale_color_manual(
      values = volcano_colors,
      drop = FALSE,
      name = "Category",
      guide = if (show_detailed_category_legend) "legend" else "none"
    ) +
    ggplot2::scale_color_manual(
      values = volcano_colors,
      drop = FALSE,
      name = "Category",
      guide = if (show_detailed_category_legend) "legend" else "none"
    ) +
    ggplot2::labs(
      title = plot_title_final,
      subtitle = plot_subtitle_final,
      x = x_axis_label,
      y = y_axis_label_final,
      caption = plot_caption_final
    ) +
    ggplot2::theme_minimal(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(family = "Arial", face = "bold", size = 24),
      plot.subtitle = ggplot2::element_text(size = base_font_size),
      plot.caption = ggplot2::element_text(
        size = base_font_size * 0.75,
        hjust = 0
      ),
      legend.position = "right",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = 14,
        face = "plain",
        color = "black"
      ),
      axis.line = if (show_axis_lines) {
        ggplot2::element_line(
          color = axis_line_color,
          linewidth = axis_line_width
        )
      } else {
        ggplot2::element_blank()
      },
      axis.ticks = if (show_axis_ticks) {
        ggplot2::element_line(
          color = axis_tick_color,
          linewidth = axis_tick_width
        )
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.length = grid::unit(axis_tick_length, "pt")
    )
  
  return(p)
}

cat("[PLOT] Volcano plotting function defined.\n")

# ----------------------------------------------------------------------------
# 9.4 Add threshold lines
# ----------------------------------------------------------------------------

add_threshold_lines <- function(p) {
  
  if (show_threshold_lines) {
    
    p <- p +
      ggplot2::geom_vline(
        xintercept = c(-log2fc_cutoff, log2fc_cutoff),
        linetype = threshold_line_type,
        linewidth = threshold_line_width
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(significance_cutoff),
        linetype = threshold_line_type,
        linewidth = threshold_line_width
      )
  }
  
  return(p)
}

cat("[PLOT] Threshold-line helper defined.\n")

# ----------------------------------------------------------------------------
# 9.5 Add plot labels
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 9.5A Highlight labeled points
# ----------------------------------------------------------------------------

add_labeled_point_highlights <- function(p, df) {
  
  if (highlight_labeled_points && any(df$plot_label != "")) {
    
    p <- p +
      ggplot2::geom_point(
        data = df %>% dplyr::filter(plot_label != ""),
        ggplot2::aes(
          x = log2FC,
          y = neg_log10_significance_for_plot,
          fill = significance_category
        ),
        inherit.aes = FALSE,
        shape = 21,
        size = labeled_point_size,
        stroke = labeled_point_stroke,
        color = labeled_point_outline_color,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = volcano_colors,
        guide = "none"
      )
  }
  
  return(p)
}

cat("[PLOT] Labeled-point highlight helper defined.\n")

# ----------------------------------------------------------------------------
# 9.5B Add plot labels
# ----------------------------------------------------------------------------
add_volcano_labels <- function(p, df) {
  
  if (any(df$plot_label != "")) {
    
    labeled_data <- df[
      !is.na(df$plot_label) &
        df$plot_label != "",
    ]
    
    # Assign fixed text colors without adding another color scale.
    # This keeps the volcano point colors intact.
    
    labeled_data$label_text_color_for_plot <- label_text_color
    
    manual_label_rows <- !is.na(labeled_data$label_category) &
      labeled_data$label_category != "Top DEG"
    
    if (any(manual_label_rows)) {
      
      manual_label_colors <- manual_label_group_text_colors_extended[
        labeled_data$label_category[manual_label_rows]
      ]
      
      manual_label_colors[is.na(manual_label_colors)] <- label_text_color
      
      labeled_data$label_text_color_for_plot[manual_label_rows] <-
        manual_label_colors
    }
    
    # Place automatic and manual labels together so they repel each other.
    p <- p +
      ggrepel::geom_text_repel(
        data = labeled_data,
        ggplot2::aes(label = plot_label),
        size = ifelse(manual_label_rows, manual_label_size, label_size),
        fontface = ifelse(manual_label_rows, manual_label_fontface, label_fontface),
        color = labeled_data$label_text_color_for_plot,
        segment.color = ifelse(manual_label_rows,
                               labeled_data$label_text_color_for_plot,
                               label_segment_color),
        box.padding = label_box_padding,
        point.padding = label_point_padding,
        max.overlaps = label_max_overlaps,
        force = label_force,
        force_pull = label_force_pull,
        seed = 42,
        max.time = 3,
        max.iter = 50000,
        min.segment.length = label_min_segment_length
      )
  }

  return(p)
}

cat("[PLOT] Label helper defined.\n")
# ----------------------------------------------------------------------------
# 9.6 Apply axis limits
# ----------------------------------------------------------------------------

apply_axis_limits <- function(p) {
  
  if (!is.null(x_axis_limits) && !is.null(y_axis_limit)) {
    
    p <- p + ggplot2::coord_cartesian(
      xlim = x_axis_limits,
      ylim = c(0, y_axis_limit)
    )
    
  } else if (!is.null(x_axis_limits)) {
    
    p <- p + ggplot2::coord_cartesian(
      xlim = x_axis_limits
    )
    
  } else if (!is.null(y_axis_limit)) {
    
    p <- p + ggplot2::coord_cartesian(
      ylim = c(0, y_axis_limit)
    )
  }
  
  return(p)
}

cat("[PLOT] Axis-limit helper defined.\n")

# ----------------------------------------------------------------------------
# 9.7 Generate volcano plot
# ----------------------------------------------------------------------------

volcano_plot <- make_volcano_plot(volcano_plot_df)
volcano_plot <- add_threshold_lines(volcano_plot)
volcano_plot <- add_labeled_point_highlights(volcano_plot, volcano_plot_df)
volcano_plot <- add_volcano_labels(volcano_plot, volcano_plot_df)
volcano_plot <- apply_axis_limits(volcano_plot)

cat("[PLOT] Volcano plot object created.\n")

# ----------------------------------------------------------------------------
# 9.8 Print plot for inspection
# ----------------------------------------------------------------------------

print(volcano_plot)

cat("[PLOT] Volcano plot printed for visual inspection.\n")

# ----------------------------------------------------------------------------
# 9.9 Section 9 checkpoint
# ----------------------------------------------------------------------------

required_section9_objects <- c(
  "plot_title_final",
  "plot_subtitle_final",
  "plot_caption_final",
  "make_volcano_plot",
  "add_threshold_lines",
  "add_labeled_point_highlights",
  "add_volcano_labels",
  "apply_axis_limits",
  "volcano_plot"
)

missing_section9_objects <- required_section9_objects[
  !vapply(required_section9_objects, exists, logical(1))
]

if (length(missing_section9_objects) > 0) {
  stop(
    "SECTION 9 checkpoint failed. Missing objects:\n",
    paste(missing_section9_objects, collapse = "\n")
  )
}

if (!inherits(volcano_plot, "ggplot")) {
  stop("SECTION 9 checkpoint failed. volcano_plot is not a ggplot object.")
}

cat("[CHECKPOINT] Section 9 completed successfully.\n")
cat("[CHECKPOINT] Volcano plot object created and printed.\n")

# ============================================================================
# SECTION 10: OUTPUT DIRECTORY CREATION AND EXPORTS
# ============================================================================

cat("\n--------------------------------------------------------\n")
cat("SECTION 10: OUTPUT DIRECTORY CREATION AND EXPORTS\n")
cat("--------------------------------------------------------\n")

run_timestamp <- make_timestamp()

if (is.null(analysis_type)) {
  analysis_type_final <- unique(volcano_plot_df$source_type)
} else {
  analysis_type_final <- analysis_type
}

file_tag <- make_file_tag(
  project_name,
  analysis_type_final,
  comparison_name_final,
  custom_output_tag
)

if (use_timestamped_output) {
  run_folder_name <- make_file_tag(file_tag, run_timestamp)
} else {
  run_folder_name <- file_tag
}

output_dir <- file.path(output_root, run_folder_name)

cat("[OUTPUT] Comparison name:", comparison_name_final, "\n")
cat("[OUTPUT] Analysis type:", analysis_type_final, "\n")
cat("[OUTPUT] File tag:", file_tag, "\n")
cat("[OUTPUT] Output directory:", output_dir, "\n")

create_dir_if_missing(output_root)
create_dir_if_missing(output_dir)

plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")

create_dir_if_missing(plots_dir)
create_dir_if_missing(tables_dir)
create_dir_if_missing(logs_dir)

cat("[OUTPUT] Output subdirectories created.\n")

# Start console logging only after the run-specific logs directory exists, so
# the log is created directly at its final destination rather than first being
# created in the working directory and copied later.
if (log_full_console_output) {
  
  console_log_file <- file.path(
    logs_dir,
    paste0(
      "GenericVolcanoPlotter_console_log_",
      format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
      ".txt"
    )
  )
  
  sink(console_log_file, split = TRUE)
  
  on.exit({
    sink()
  }, add = TRUE)
  
  cat("[LOG] Console output from this point forward is being saved to:\n")
  cat(console_log_file, "\n")
}

if (export_png) {
  
  png_file <- file.path(
    plots_dir,
    paste0(file_tag, "_volcano_plot.png")
  )
  
  ggplot2::ggsave(
    filename = png_file,
    plot = volcano_plot,
    width = png_width,
    height = png_height,
    dpi = png_dpi,
    units = "in"
  )
  
  cat("[EXPORT] PNG volcano plot saved:", png_file, "\n")
  
} else {
  
  png_file <- NA_character_
  cat("[EXPORT] PNG export skipped.\n")
}

if (export_tiff) {
  
  tiff_file <- file.path(
    plots_dir,
    paste0(file_tag, "_volcano_plot.tiff")
  )
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = volcano_plot,
    width = tiff_width,
    height = tiff_height,
    dpi = tiff_dpi,
    units = "in",
    compression = "lzw"
  )
  
  cat("[EXPORT] TIFF volcano plot saved:", tiff_file, "\n")
  
} else {
  
  tiff_file <- NA_character_
  cat("[EXPORT] TIFF export skipped.\n")
}

if (export_svg) {
  
  svg_file <- file.path(
    plots_dir,
    paste0(file_tag, "_volcano_plot.svg")
  )
  
  ggplot2::ggsave(
    filename = svg_file,
    plot = volcano_plot,
    width = svg_width,
    height = svg_height,
    units = "in"
  )
  
  cat("[EXPORT] SVG volcano plot saved:", svg_file, "\n")
  
} else {
  
  svg_file <- NA_character_
  cat("[EXPORT] SVG export skipped.\n")
}

if (export_pdf) {
  
  pdf_file <- file.path(
    plots_dir,
    paste0(file_tag, "_volcano_plot.pdf")
  )
  
  ggplot2::ggsave(
    filename = pdf_file,
    device = function(filename, width, height, ...) {
      grDevices::quartz(type = "pdf", file = filename,
                       width = width, height = height, ...)
    },
    plot = volcano_plot,
    width = pdf_width,
    height = pdf_height,
    units = "in"
  )
  
  cat("[EXPORT] PDF volcano plot saved:", pdf_file, "\n")
  
} else {
  
  pdf_file <- NA_character_
  cat("[EXPORT] PDF export skipped.\n")
}

if (save_standardized_table) {
  
  standardized_table_file <- file.path(
    tables_dir,
    paste0(file_tag, "_standardized_volcano_table.csv")
  )
  
  readr::write_csv(volcano_df, standardized_table_file)
  
  cat("[EXPORT] Standardized volcano table saved:", standardized_table_file, "\n")
  
} else {
  
  standardized_table_file <- NA_character_
  cat("[EXPORT] Standardized volcano table export skipped.\n")
}

plotting_table_file <- file.path(
  tables_dir,
  paste0(file_tag, "_plotting_ready_volcano_table.csv")
)

readr::write_csv(volcano_plot_df, plotting_table_file)

cat("[EXPORT] Plotting-ready volcano table saved:", plotting_table_file, "\n")

if (save_significant_table) {
  
  significant_table_file <- file.path(
    tables_dir,
    paste0(file_tag, "_significant_hits.csv")
  )
  
  readr::write_csv(significant_hits_df, significant_table_file)
  
  cat("[EXPORT] Significant hits table saved:", significant_table_file, "\n")
  
} else {
  
  significant_table_file <- NA_character_
  cat("[EXPORT] Significant hits table export skipped.\n")
}

if (save_labeled_table) {
  
  labeled_table_file <- file.path(
    tables_dir,
    paste0(file_tag, "_labeled_features.csv")
  )
  
  readr::write_csv(selected_label_report, labeled_table_file)
  
  cat("[EXPORT] Labeled features table saved:", labeled_table_file, "\n")
  
} else {
  
  labeled_table_file <- NA_character_
  cat("[EXPORT] Labeled features table export skipped.\n")
}

column_name_map_file <- file.path(
  tables_dir,
  paste0(file_tag, "_column_name_map.csv")
)

readr::write_csv(column_name_map, column_name_map_file)

cat("[EXPORT] Column name map saved:", column_name_map_file, "\n")

if (log_full_console_output && file.exists(console_log_file)) {
  
  # Retain this alias for run-summary compatibility. No copy is needed because
  # console_log_file already points to the final file inside logs_dir.
  console_log_file_final <- console_log_file
  
  cat("[EXPORT] Console log saved directly to:", console_log_file_final, "\n")
  
} else {
  
  console_log_file_final <- NA_character_
}

run_summary <- tibble::tibble(
  setting = c(
    "script_name",
    "script_version",
    "input_file",
    "input_type_detected",
    "analysis_type_final",
    "comparison_name_final",
    "significance_column",
    "significance_cutoff",
    "log2fc_cutoff",
    "rows_raw",
    "rows_standardized",
    "rows_plotted",
    "up_hits",
    "down_hits",
    "significant_low_effect_hits",
    "not_significant_hits",
    "labels_selected",
    "label_ranking_method",
    "y_axis_limit",
    "run_timestamp",
    "output_dir",
    "manual_labels_enabled",
    "manual_label_features",
    "manual_labels_add_to_auto_labels",
    "console_log_file",
    "console_log_file_final"
  ),
  value = as.character(c(
    "GenericVolcanoPlotter",
    script_version,
    input_metadata$input_file_name,
    detected_input_type,
    analysis_type_final,
    comparison_name_final,
    significance_col,
    significance_cutoff,
    log2fc_cutoff,
    input_metadata$n_rows_raw,
    nrow(volcano_df),
    nrow(volcano_plot_df),
    sum(volcano_plot_df$significance_category == "Up", na.rm = TRUE),
    sum(volcano_plot_df$significance_category == "Down", na.rm = TRUE),
    sum(volcano_plot_df$significance_category == "Significant_low_effect", na.rm = TRUE),
    sum(volcano_plot_df$significance_category == "Not_significant", na.rm = TRUE),
    sum(volcano_plot_df$plot_label != ""),
    label_ranking_method,
    ifelse(is.null(y_axis_limit), "none", y_axis_limit),
    run_timestamp,
    output_dir,
    manual_labels_enabled,
    paste(manual_label_features, collapse = ";"),
    manual_labels_add_to_auto_labels,
    console_log_file,
    console_log_file_final
  ))
)

run_summary_file <- file.path(
  logs_dir,
  paste0(file_tag, "_run_summary.csv")
)

readr::write_csv(run_summary, run_summary_file)

cat("[EXPORT] Run summary saved:", run_summary_file, "\n")

if (save_session_info) {
  
  session_info_file <- file.path(
    logs_dir,
    paste0(file_tag, "_session_info.txt")
  )
  
  writeLines(
    capture.output(sessionInfo()),
    con = session_info_file
  )
  
  cat("[EXPORT] Session info saved:", session_info_file, "\n")
  
} else {
  
  session_info_file <- NA_character_
  cat("[EXPORT] Session info export skipped.\n")
}

required_section10_objects <- c(
  "run_timestamp",
  "comparison_name_final",
  "analysis_type_final",
  "file_tag",
  "output_dir",
  "plots_dir",
  "tables_dir",
  "logs_dir",
  "run_summary",
  "run_summary_file"
)

missing_section10_objects <- required_section10_objects[
  !vapply(required_section10_objects, exists, logical(1))
]

if (length(missing_section10_objects) > 0) {
  stop(
    "SECTION 10 checkpoint failed. Missing objects:\n",
    paste(missing_section10_objects, collapse = "\n")
  )
}

if (!dir.exists(output_dir)) {
  stop("SECTION 10 checkpoint failed. output_dir was not created.")
}

if (export_png && !file.exists(png_file)) {
  stop("SECTION 10 checkpoint failed. PNG file was not created.")
}

if (export_tiff && !file.exists(tiff_file)) {
  stop("SECTION 10 checkpoint failed. TIFF file was not created.")
}

if (export_svg && !file.exists(svg_file)) {
  stop("SECTION 10 checkpoint failed. SVG file was not created.")
}

if (export_pdf && !file.exists(pdf_file)) {
  stop("SECTION 10 checkpoint failed. PDF file was not created.")
}

if (save_standardized_table && !file.exists(standardized_table_file)) {
  stop("SECTION 10 checkpoint failed. Standardized table was not created.")
}

if (!file.exists(plotting_table_file)) {
  stop("SECTION 10 checkpoint failed. Plotting-ready table was not created.")
}

if (save_significant_table && !file.exists(significant_table_file)) {
  stop("SECTION 10 checkpoint failed. Significant hits table was not created.")
}

if (save_labeled_table && !file.exists(labeled_table_file)) {
  stop("SECTION 10 checkpoint failed. Labeled features table was not created.")
}

if (!file.exists(column_name_map_file)) {
  stop("SECTION 10 checkpoint failed. Column name map was not created.")
}

if (!file.exists(run_summary_file)) {
  stop("SECTION 10 checkpoint failed. Run summary was not created.")
}

if (save_session_info && !file.exists(session_info_file)) {
  stop("SECTION 10 checkpoint failed. Session info file was not created.")
}

cat("[CHECKPOINT] Section 10 completed successfully.\n")
cat("[CHECKPOINT] Output files exported successfully.\n")
cat("[OUTPUT] Final output directory:", output_dir, "\n")

# ============================================================================
# SECTION 11: FINAL RUN SUMMARY
# ============================================================================

cat("\n========================================================\n")
cat("GENERIC VOLCANO PLOTTER RUN COMPLETE\n")
cat("========================================================\n")

cat("[SUMMARY] Input file:", input_metadata$input_file_name, "\n")
cat("[SUMMARY] Detected input type:", detected_input_type, "\n")
cat("[SUMMARY] Comparison:", comparison_name_final, "\n")
cat("[SUMMARY] Significance column:", significance_col, "\n")
cat("[SUMMARY] Significance cutoff:", significance_cutoff, "\n")
cat("[SUMMARY] log2FC cutoff:", log2fc_cutoff, "\n")
cat("[SUMMARY] Rows imported:", input_metadata$n_rows_raw, "\n")
cat("[SUMMARY] Rows plotted:", nrow(volcano_plot_df), "\n")
cat("[SUMMARY] Up hits:", sum(volcano_plot_df$significance_category == "Up", na.rm = TRUE), "\n")
cat("[SUMMARY] Down hits:", sum(volcano_plot_df$significance_category == "Down", na.rm = TRUE), "\n")
cat("[SUMMARY] Significant low-effect hits:", 
    sum(volcano_plot_df$significance_category == "Significant_low_effect", na.rm = TRUE), "\n")
cat("[SUMMARY] Not significant:", 
    sum(volcano_plot_df$significance_category == "Not_significant", na.rm = TRUE), "\n")
cat("[SUMMARY] Labels selected:", sum(volcano_plot_df$plot_label != ""), "\n")
cat("[SUMMARY] Output directory:", output_dir, "\n")

cat("========================================================\n")
cat("Run completed successfully.\n")
cat("========================================================\n")
