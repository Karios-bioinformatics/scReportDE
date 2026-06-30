# scReportDE: Utility Functions ------------------------------------------------


#' Null coalescing operator
#'
#' Returns \code{x} unless it is \code{NULL}, then returns \code{y}.
#' @param x A value (possibly NULL)
#' @param y Fallback value
#' @return \code{x} if not NULL, otherwise \code{y}
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Format a numeric value for display
#'
#' @param x Numeric value
#' @param digits Number of decimal places (default 0)
#' @return Character string
#' @keywords internal
fmt_num <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}


#' Format a p-value for display
#'
#' Scientific notation for values < 1e-4, otherwise 4 decimal places.
#' @param x Numeric p-value
#' @return Character string
#' @keywords internal
fmt_pval <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-4) return(format(x, digits = 2, scientific = TRUE))
  sprintf("%.4f", x)
}


#' Normalise a DE result data.frame to a standard schema
#'
#' Detects column names from Seurat::FindMarkers / FindAllMarkers output
#' (or user-supplied equivalents) and renames them to a canonical set so
#' downstream code sees predictable column names.
#'
#' Canonical output columns:
#' \itemize{
#'   \item \code{gene} — gene identifier (from \code{gene} column or rownames)
#'   \item \code{cluster} — cluster label (FindAllMarkers only; \code{NA} if absent)
#'   \item \code{avg_log2FC} — log2 fold change
#'   \item \code{p_val} — raw p-value
#'   \item \code{p_val_adj} — adjusted p-value
#'   \item \code{pct.1} — fraction of cells expressing in group 1
#'   \item \code{pct.2} — fraction of cells expressing in group 2
#' }
#'
#' @param df A data.frame of DE results
#' @return A list with elements \code{df} (normalised data.frame) and
#'   \code{warnings} (character vector of issues found)
#' @keywords internal
normalize_de_df <- function(df) {
  warns <- character(0)

  if (!inherits(df, "data.frame")) {
    stop("de_df must be a data.frame, got ", class(df)[1])
  }

  # ---- gene column ----
  if (!"gene" %in% colnames(df)) {
    rn <- rownames(df)
    if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) {
      df$gene <- rn
      rownames(df) <- NULL
    } else {
      stop("Cannot identify gene column: no 'gene' column and rownames look like row numbers")
    }
  }

  # ---- Rename logFC ----
  if ("avg_log2FC" %in% colnames(df)) {
    # already canonical
  } else if ("avg_logFC" %in% colnames(df)) {
    colnames(df)[colnames(df) == "avg_logFC"] <- "avg_log2FC"
    warns <- c(warns, "Renamed 'avg_logFC' to 'avg_log2FC'")
  } else if ("log2FoldChange" %in% colnames(df)) {
    colnames(df)[colnames(df) == "log2FoldChange"] <- "avg_log2FC"
    warns <- c(warns, "Renamed 'log2FoldChange' to 'avg_log2FC'")
  } else {
    warns <- c(warns, "No logFC column found (avg_log2FC / avg_logFC / log2FoldChange)")
  }

  # ---- Check required columns ----
  required <- c("p_val", "p_val_adj")
  for (col in required) {
    if (!col %in% colnames(df)) {
      warns <- c(warns, paste0("Missing column: '", col, "'"))
    }
  }

  # ---- Fill optional columns with NA ----
  optional <- c("cluster", "pct.1", "pct.2")
  for (col in optional) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA_real_
      warns <- c(warns, paste0("Column '", col, "' not found — filled with NA"))
    }
  }

  # ---- Ensure gene is character ----
  df$gene <- as.character(df$gene)

  list(df = df, warnings = warns)
}


#' Safe wrapper around Seurat differential expression
#'
#' Wraps \code{Seurat::FindMarkers} and \code{Seurat::FindAllMarkers} in
#' \code{tryCatch} so that upstream errors (e.g. missing ident, empty group)
#' become informative warnings instead of hard stops.
#'
#' @param seurat_obj A Seurat object
#' @param mode \code{"pairwise"} or \code{"all_markers"}
#' @param ident_1 Identity class 1 (for pairwise)
#' @param ident_2 Identity class 2 (for pairwise)
#' @param assay Assay to use
#' @param slot Slot to use
#' @param test_use Statistical test
#' @param logfc_threshold LogFC threshold
#' @param min_pct Minimum detection fraction
#' @param only_pos Return only positive markers
#' @return A list with elements \code{de_df} (data.frame or NULL) and
#'   \code{warnings} (character vector)
#' @keywords internal
safe_compute_de <- function(seurat_obj, mode, ident_1, ident_2,
                            assay, slot, test_use,
                            logfc_threshold, min_pct, only_pos) {
  warns <- character(0)
  de_df <- NULL

  result <- tryCatch({
    if (mode == "pairwise") {
      if (is.null(ident_1) || is.null(ident_2)) {
        stop("mode='pairwise' requires ident_1 and ident_2")
      }
      Seurat::FindMarkers(
        object          = seurat_obj,
        ident.1         = ident_1,
        ident.2         = ident_2,
        assay           = assay,
        slot            = slot,
        test.use        = test_use,
        logfc.threshold = logfc_threshold,
        min.pct         = min_pct,
        only.pos        = only_pos
      )
    } else if (mode == "all_markers") {
      Seurat::FindAllMarkers(
        object          = seurat_obj,
        assay           = assay,
        slot            = slot,
        test.use        = test_use,
        logfc.threshold = logfc_threshold,
        min.pct         = min_pct,
        only.pos        = only_pos
      )
    } else {
      stop("Unknown mode: ", mode)
    }
  }, error = function(e) {
    warns <<- c(warns, paste("DE computation failed:", e$message))
    NULL
  })

  if (is.null(result) || nrow(result) == 0) {
    warns <- c(warns, "DE result is empty — no differentially expressed genes found")
    return(list(de_df = NULL, warnings = warns))
  }

  # Normalise column names
  norm <- normalize_de_df(result)
  list(de_df = norm$df, warnings = c(warns, norm$warnings))
}
