# scReportDE: Main API — build_screport_de() -----------------------------------
#
# Single entry-point for the DE report generator.
# Accepts either a Seurat object (runs FindMarkers/FindAllMarkers internally)
# or a pre-computed DE data.frame, and produces an interactive HTML report
# with volcano plot, DE table, dot plot, and violin plot.


#' Generate a Differential Expression HTML Report
#'
#' Takes a Seurat object (to compute DE) or a pre-computed DE data.frame,
#' and produces a self-contained interactive HTML report with:
#' \itemize{
#'   \item Overview — summary cards (n DEGs, up/down/sig counts, comparison)
#'   \item Volcano Plot — interactive plotly with logFC vs -log10(p-value)
#'   \item DE Table — searchable, sortable, paginated DT table
#'   \item Dot Plot — top N genes expression dot plot (needs Seurat object)
#'   \item Violin Plot — top gene expression distribution (needs Seurat object)
#'   \item Method Info — all input parameters and generation metadata
#' }
#'
#' @param seurat_obj A Seurat object. Required when \code{de_df} is NULL.
#' @param de_df Optional pre-computed DE data.frame. Takes priority over
#'   \code{seurat_obj} when both are provided. Must contain at minimum a
#'   gene identifier, logFC, and p-value columns. Seurat FindMarkers /
#'   FindAllMarkers output is detected automatically.
#' @param group_col Optional metadata column name to set as identity
#'   before running DE (e.g. \code{"condition"}).
#' @param ident_1 Identity class 1 for pairwise mode (e.g. \code{"Treatment"}).
#' @param ident_2 Identity class 2 for pairwise mode (e.g. \code{"Control"}).
#' @param mode \code{"pairwise"} (default) runs \code{Seurat::FindMarkers};
#'   \code{"all_markers"} runs \code{Seurat::FindAllMarkers}. Ignored when
#'   \code{de_df} is provided.
#' @param assay Assay to pull expression data from. Passed to Seurat DE
#'   functions and used for dot/violin plots. Default: \code{NULL} (Seurat default).
#' @param slot Expression data slot. Default: \code{"data"}.
#' @param test_use Statistical test for \code{FindMarkers}/\code{FindAllMarkers}.
#'   Default: \code{"wilcox"}.
#' @param logfc_threshold LogFC threshold for DE filtering. Default: 0.25.
#' @param min_pct Minimum fraction of cells expressing a gene. Default: 0.1.
#' @param only_pos Return only positive markers. Default: FALSE.
#' @param top_n Number of top genes labelled in volcano / used in dot plot.
#'   Default: 20.
#' @param volcano_label_top_n Max gene labels on Volcano Plot. When NULL
#'   (default), auto-detects: 5 for marker-only data, 8 for bidirectional.
#'   Set explicitly to override. Full gene info available in hover.
#' @param dotplot_identity_layers Character vector of meta.data columns to use
#'   as x-axis grouping variables for the interactive DotPlot. When NULL, defaults
#'   are auto-detected (see \code{.collect_dotplot_identity_layers}).
#' @param dotplot_marker_pool_top_n Integer. Top N marker genes per identity value
#'   to include in the precomputed gene pool. Default: 50.
#' @param dotplot_pool_max_genes Integer. Maximum number of genes in the
#'   precomputed pool. Default: 500.
#' @param dotplot_top_n Integer. Default number of top markers shown per group.
#'   Default: 10.
#' @param dotplot_max_display_genes Integer. Maximum genes shown on the plot.
#'   Default: 80.
#' @param dotplot_direction Character. Default marker direction filter:
#'   \code{"up"}, \code{"down"}, or \code{"both"}. Default: \code{"up"}.
#' @param dotplot_extra_genes Character vector of extra genes to force into the
#'   precomputed gene pool. Default: \code{NULL}.
#' @param dotplot_size_min Numeric. Minimum bubble diameter (Plotly units).
#'   Default: 3.
#' @param dotplot_size_max Numeric. Maximum bubble diameter. Default: 14.
#' @param output_file Path to the output HTML file.
#'   Default: \code{"scReport_DE.html"}.
#' @param title Report title shown in the header.
#'   Default: \code{"Differential Expression Report"}.
#' @param self_contained If TRUE, attempt to embed all dependencies.
#'   v0.1.0: reserved — not yet implemented. When FALSE, a sharing notice
#'   is included in the report.
#'
#' @return Invisibly, a list with elements:
#'   \item{de_df}{The normalised DE data.frame used in the report}
#'   \item{output_file}{Normalised path to the generated HTML file}
#'   \item{metadata}{List of all parameters and computed info}
#'   \item{warnings}{Character vector of non-fatal warnings encountered}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(scReportDE)
#' library(Seurat)
#'
#' # Pairwise DE from a Seurat object
#' obj <- readRDS("my_seurat.rds")
#' build_screport_de(
#'   seurat_obj = obj,
#'   group_col  = "condition",
#'   ident_1    = "Treatment",
#'   ident_2    = "Control",
#'   mode       = "pairwise",
#'   output_file = "DE_Treatment_vs_Control.html"
#' )
#'
#' # Using a pre-computed DE data.frame
#' de_results <- FindMarkers(obj, ident.1 = "A", ident.2 = "B")
#' build_screport_de(de_df = de_results)
#' }
build_screport_de <- function(
    seurat_obj      = NULL,
    de_df           = NULL,
    group_col       = NULL,
    ident_1         = NULL,
    ident_2         = NULL,
    mode            = c("pairwise", "all_markers"),
    assay           = NULL,
    slot            = "data",
    test_use        = "wilcox",
    logfc_threshold = 0.25,
    min_pct         = 0.1,
    only_pos        = FALSE,
    top_n           = 20,
    volcano_label_top_n         = NULL,
    dotplot_identity_layers     = NULL,
    dotplot_marker_pool_top_n   = 50,
    dotplot_pool_max_genes      = 500,
    dotplot_top_n               = 10,
    dotplot_max_display_genes   = 80,
    dotplot_direction           = "up",
    dotplot_extra_genes         = NULL,
    dotplot_size_min            = 3,
    dotplot_size_max            = 14,
    output_file     = "scReport_DE.html",
    title           = "Differential Expression Report",
    self_contained  = FALSE) {

  mode <- match.arg(mode)
  generation_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  all_warnings <- character(0)

  # ---- 1. Validate inputs ----
  if (is.null(de_df) && is.null(seurat_obj)) {
    stop("Either seurat_obj or de_df must be provided")
  }

  # ---- 2. Compute or ingest DE results ----
  message("Step 1/4: Computing / ingesting DE results...")
  de_result <- compute_de(
    seurat_obj      = seurat_obj,
    de_df           = de_df,
    group_col       = group_col,
    ident_1         = ident_1,
    ident_2         = ident_2,
    mode            = mode,
    assay           = assay,
    slot            = slot,
    test_use        = test_use,
    logfc_threshold = logfc_threshold,
    min_pct         = min_pct,
    only_pos        = only_pos
  )

  de_df_norm <- de_result$de_df
  all_warnings <- c(all_warnings, de_result$warnings)

  # ---- 3. Build comparison label ----
  if (!is.null(de_result$compute_info$source) &&
      grepl("Provided", de_result$compute_info$source)) {
    comparison_label <- "Provided de_df"
  } else if (mode == "pairwise") {
    comparison_label <- paste(ident_1 %||% "?", "vs", ident_2 %||% "?")
  } else {
    comparison_label <- "All clusters (one-vs-rest)"
  }

  # ---- 4. Build method metadata ----
  method_meta <- build_method_meta(de_result$compute_info$source, list(
    mode            = mode,
    assay           = assay %||% "default",
    slot            = slot,
    test_use        = test_use,
    logfc_threshold = logfc_threshold,
    min_pct         = min_pct,
    only_pos        = only_pos,
    group_col       = group_col,
    ident_1         = ident_1,
    ident_2         = ident_2
  ))

  # ---- 5. Generate plots (each wrapped in tryCatch) ----
  message("Step 2/4: Generating plots...")

  # Volcano plot
  message("  - Volcano plot...")
  volcano_widget <- tryCatch({
    plot_volcano(de_df_norm, top_n = top_n, alpha = 0.05,
                 comparison_label = comparison_label,
                 label_top_n = volcano_label_top_n)
  }, error = function(e) {
    all_warnings <<- c(all_warnings, paste("Volcano plot:", e$message))
    NULL
  })

  # Dot plot — interactive panel
  message("  - Dot plot (interactive)...")
  dotplot_widget <- tryCatch({
    .build_interactive_marker_dotplot_panel(
      seurat_obj       = seurat_obj,
      marker_df        = de_df_norm,
      identity_layers  = dotplot_identity_layers,
      group_col        = group_col,
      assay            = assay,
      slot             = slot,
      marker_pool_top_n = dotplot_marker_pool_top_n,
      pool_max_genes   = dotplot_pool_max_genes,
      top_n            = dotplot_top_n,
      max_display_genes = dotplot_max_display_genes,
      direction        = dotplot_direction,
      extra_genes      = dotplot_extra_genes,
      size_min         = dotplot_size_min,
      size_max         = dotplot_size_max
    )
  }, error = function(e) {
    all_warnings <<- c(all_warnings, paste("Dot plot:", e$message))
    NULL
  })

  # Violin plot
  message("  - Violin plot...")
  violin_widget <- tryCatch({
    plot_violin(seurat_obj, de_df_norm, group_col = group_col,
                assay = assay, slot = slot)
  }, error = function(e) {
    all_warnings <<- c(all_warnings, paste("Violin plot:", e$message))
    NULL
  })

  # ---- 6. Build overview cards ----
  overview_cards <- build_overview_cards(de_df_norm, comparison_label)

  # ---- 7. Build DE table ----
  message("Step 3/4: Building DE table...")
  de_table_widget <- build_de_table(de_df_norm)

  # ---- 8. Assemble HTML ----
  message("Step 4/4: Assembling HTML report...")
  build_html(
    de_df             = de_df_norm,
    overview_cards    = overview_cards,
    volcano_widget    = volcano_widget,
    de_table_widget   = de_table_widget,
    dotplot_widget    = dotplot_widget,
    violin_widget     = violin_widget,
    method_meta       = method_meta,
    generation_time   = generation_time,
    output            = output_file,
    title             = title,
    self_contained    = self_contained,
    warnings          = all_warnings
  )

  # ---- 9. Console sharing notice ----
  if (!self_contained) {
    message(
      "\n[!] Sharing notice:\n",
      "    This report is NOT self-contained. To share it:\n",
      "    (1) Send the HTML file together with the generated 'lib' folder, or\n",
      "    (2) Re-generate with self_contained = TRUE (not yet available in v0.1.0).\n"
    )
  }

  # ---- 10. Return ----
  invisible(list(
    de_df       = de_df_norm,
    output_file = normalizePath(output_file, mustWork = FALSE),
    metadata    = method_meta,
    warnings    = all_warnings
  ))
}
