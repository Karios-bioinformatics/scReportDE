# scReportDE: DE Computation ---------------------------------------------------
#
# Orchestrates DE computation: decides whether to use a user-supplied de_df,
# run Seurat::FindMarkers, or run Seurat::FindAllMarkers.


#' Compute or ingest differential expression results
#'
#' If \code{de_df} is provided it is validated and normalised directly.
#' Otherwise, Seurat DE is called according to \code{mode}.
#'
#' @param seurat_obj A Seurat object (required when de_df is NULL)
#' @param de_df Optional pre-computed DE data.frame
#' @param group_col Optional identity column in seurat_obj meta.data
#' @param ident_1 Identity 1 (for pairwise mode)
#' @param ident_2 Identity 2 (for pairwise mode)
#' @param mode \code{"pairwise"} or \code{"all_markers"}
#' @param assay Assay name
#' @param slot Expression slot
#' @param test_use Statistical test for FindMarkers/FindAllMarkers
#' @param logfc_threshold LogFC threshold
#' @param min_pct Minimum detection fraction
#' @param only_pos Return only positive markers
#' @return A list with \code{de_df}, \code{warnings}, and \code{compute_info}
#' @keywords internal
compute_de <- function(seurat_obj, de_df, group_col,
                       ident_1, ident_2, mode,
                       assay, slot, test_use,
                       logfc_threshold, min_pct, only_pos) {
  warns <- character(0)
  source <- NULL

  # ---- Path 1: User-provided de_df ----
  if (!is.null(de_df)) {
    message("Using provided de_df (skipping DE computation)")
    norm <- normalize_de_df(de_df)
    warns <- c(warns, norm$warnings)
    source <- "Provided de_df"
    return(list(
      de_df        = norm$df,
      warnings     = warns,
      compute_info = list(source = source)
    ))
  }

  # ---- Path 2: Compute from Seurat ----
  if (is.null(seurat_obj)) {
    stop("Either seurat_obj or de_df must be provided")
  }

  # Set identity if group_col is given
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(seurat_obj[[]])) {
      stop("group_col '", group_col, "' not found in seurat_obj meta.data")
    }
    Seurat::Idents(seurat_obj) <- group_col
  }

  message("Running DE: mode=", mode, " ...")
  res <- safe_compute_de(
    seurat_obj      = seurat_obj,
    mode            = mode,
    ident_1         = ident_1,
    ident_2         = ident_2,
    assay           = assay,
    slot            = slot,
    test_use        = test_use,
    logfc_threshold = logfc_threshold,
    min_pct         = min_pct,
    only_pos        = only_pos
  )

  warns <- c(warns, res$warnings)
  source <- if (mode == "pairwise") {
    paste("Seurat::FindMarkers:", ident_1, "vs", ident_2)
  } else {
    "Seurat::FindAllMarkers"
  }

  list(
    de_df        = res$de_df,
    warnings     = warns,
    compute_info = list(source = source)
  )
}


#' Build a metadata list capturing all parameters for the Method Info section
#'
#' @param source Source description (from compute_de)
#' @param params Named list of user-facing parameters
#' @return A list of display-ready key-value pairs
#' @keywords internal
build_method_meta <- function(source, params) {
  list(
    "Input source"      = source %||% "Unknown",
    "Mode"              = params$mode %||% "N/A",
    "Assay"             = params$assay %||% "default",
    "Slot"              = params$slot %||% "data",
    "Test"              = params$test_use %||% "wilcox",
    "LogFC threshold"   = as.character(params$logfc_threshold %||% 0.25),
    "Min pct"           = as.character(params$min_pct %||% 0.1),
    "Only positive"     = as.character(params$only_pos %||% FALSE),
    "Group column"      = if (is.null(params$group_col)) "(none)" else params$group_col,
    "ident.1"           = params$ident_1 %||% "N/A",
    "ident.2"           = params$ident_2 %||% "N/A",
    "Package version"   = "scReportDE v0.1.0"
  )
}
