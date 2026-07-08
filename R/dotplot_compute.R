# scReportDE: Marker DotPlot — Compute Functions -------------------------------
#
# Identity layer collection, marker normalisation, gene pool construction,
# dotplot expression data precomputation, bubble size scaling.
#
# Design principle (Karios, 2026-07-04):
#   - First version is marker-driven, not full-transcriptome.
#   - Gene pool is built from top K markers per identity value + extra_genes.
#   - Only this controlled gene pool is precomputed and embedded as JSON.
#   - Report-side interactivity is limited to the precomputed gene pool.


# ---- Identity Layer Collection -----------------------------------------------

#' Collect candidate DotPlot identity layers from a Seurat object
#'
#' Returns a character vector of meta.data column names that could serve as
#' x-axis grouping variables for DotPlot.
#'
#' Priority:
#'   1. \code{user_layers} (explicit user input) always takes precedence.
#'   2. If \code{user_layers} is NULL or empty, auto-detect from meta.data.
#'
#' Auto-detection rules:
#'   - Current \code{Idents(seurat_obj)}
#'   - \code{seurat_clusters} if present
#'   - Columns matching cluster / annotation / celltype / cell_type / sample / condition
#'   - Columns with the prefix \code{subcluster_}
#'
#' Subcluster fields are NOT merged with the main cluster; each subcluster
#' field becomes its own independent identity layer.
#'
#' @param seurat_obj A Seurat object.
#' @param user_layers Optional character vector of meta.data column names.
#'   When provided, these are used verbatim (validation is applied).
#' @return Character vector of usable identity layer column names.
#' @keywords internal
.collect_dotplot_identity_layers <- function(seurat_obj, user_layers = NULL) {
  if (is.null(seurat_obj)) return(character(0))

  # Explicit user input overrides everything
  if (!is.null(user_layers) && length(user_layers) > 0) {
    valid <- intersect(user_layers, colnames(seurat_obj[[]]))
    if (length(valid) == 0) {
      warning("dotplot_identity_layers: none of the specified columns found in meta.data")
      return(character(0))
    }
    missing <- setdiff(user_layers, valid)
    if (length(missing) > 0) {
      warning("dotplot_identity_layers: columns not found in meta.data: ",
              paste(missing, collapse = ", "))
    }
    return(valid)
  }

  # Auto-detect
  meta_cols <- colnames(seurat_obj[[]])
  candidates <- character(0)

  # Current Idents
  idents_name <- "seurat_clusters"  # fallback name
  tryCatch({
    idents_name <- as.character(SeuratObject::Idents(seurat_obj))
    # Actually we need the *column name*, not the values.
    # Seurat stores active.ident name differently. Just include seurat_clusters if present.
  }, error = function(e) NULL)

  # seurat_clusters
  if ("seurat_clusters" %in% meta_cols) {
    candidates <- c(candidates, "seurat_clusters")
  }

  # Keyword-based detection
  keyword_patterns <- c("cluster", "annotation", "celltype", "cell_type",
                         "sample", "condition")
  for (pat in keyword_patterns) {
    hits <- grep(pat, meta_cols, ignore.case = TRUE, value = TRUE)
    for (h in hits) {
      if (!h %in% candidates) {
        candidates <- c(candidates, h)
      }
    }
  }

  # subcluster_ prefix
  sub_hits <- grep("^subcluster_", meta_cols, value = TRUE)
  for (h in sub_hits) {
    if (!h %in% candidates) {
      candidates <- c(candidates, h)
    }
  }

  # Always include active identity if we can determine it
  active_ident <- NULL
  tryCatch({
    active_ident <- SeuratObject::DefaultDimRed(seurat_obj)
  }, error = function(e) NULL)
  # Fallback: use active.assay or just return what we found

  unique(candidates)
}


# ---- Marker DF Normalisation ------------------------------------------------

#' Normalise a marker/DE data.frame to the canonical format
#'
#' Ensures the marker data.frame contains at minimum:
#'   gene, identity_layer, identity_value, avg_log2FC, p_val_adj
#'
#' Compatibility rules (from spec):
#'   - If rownames are gene names but no 'gene' column → add gene column
#'   - If no identity_layer column → set to "seurat_clusters"
#'   - If no identity_value but 'cluster' exists → use cluster
#'   - If no identity_value but 'group' exists → use group
#'   - If no p_val_adj but p_val exists → use p_val as p_val_adj
#'   - If no avg_log2FC but avg_logFC exists → rename
#'
#' @param marker_df A data.frame of marker/DE results.
#' @param default_identity Character. Identity layer name to use when
#'   the data.frame lacks an identity_layer column.
#' @return A normalised data.frame, or NULL if the input is NULL/empty.
#' @keywords internal
.normalize_marker_df <- function(marker_df, default_identity = "seurat_clusters") {
  if (is.null(marker_df) || nrow(marker_df) == 0) return(NULL)

  df <- marker_df

  # rownames → gene
  if (!"gene" %in% colnames(df)) {
    rn <- rownames(df)
    if (!is.null(rn) && length(rn) == nrow(df) && !all(grepl("^[0-9]+$", rn))) {
      df$gene <- rn
    }
  }

  if (!"gene" %in% colnames(df)) {
    warning("marker_df has no 'gene' column and rownames cannot be used as gene names")
    return(NULL)
  }

  # identity_layer
  if (!"identity_layer" %in% colnames(df)) {
    df$identity_layer <- default_identity
  }

  # identity_value
  if (!"identity_value" %in% colnames(df)) {
    if ("cluster" %in% colnames(df)) {
      df$identity_value <- as.character(df$cluster)
    } else if ("group" %in% colnames(df)) {
      df$identity_value <- as.character(df$group)
    } else {
      warning("marker_df has no identity_value / cluster / group column")
      return(NULL)
    }
  }

  # avg_log2FC
  if (!"avg_log2FC" %in% colnames(df) && "avg_logFC" %in% colnames(df)) {
    df$avg_log2FC <- df$avg_logFC
  }

  # p_val_adj
  if (!"p_val_adj" %in% colnames(df) && "p_val" %in% colnames(df)) {
    df$p_val_adj <- df$p_val
  }

  # Ensure identity_value is character
  df$identity_value <- as.character(df$identity_value)

  df
}


# ---- Gene Pool Construction --------------------------------------------------

#' Build the controlled gene pool for DotPlot precomputation
#'
#' For each identity_layer × identity_value, extracts the top K marker genes.
#' Merges and deduplicates across all identity values. Optionally adds
#' user-specified extra genes. Caps at \code{pool_max_genes}.
#'
#' @param marker_df Normalised marker data.frame (from \code{.normalize_marker_df}).
#' @param identity_layers Character vector of identity layer names.
#' @param extra_genes Optional character vector of extra genes to force into pool.
#' @param pool_top_n Integer. Top N markers to take per identity value. Default 50.
#' @param pool_max_genes Integer. Absolute cap on gene pool size. Default 500.
#' @return Character vector of gene names (the gene pool).
#' @keywords internal
.build_gene_pool <- function(marker_df, identity_layers, extra_genes = NULL,
                              pool_top_n = 50, pool_max_genes = 500) {
  gene_set <- character(0)

  if (!is.null(marker_df) && nrow(marker_df) > 0) {
    # For each identity_layer, for each identity_value, take top N by rank or p_val_adj
    for (lyr in identity_layers) {
      sub <- marker_df[marker_df$identity_layer == lyr, , drop = FALSE]
      if (nrow(sub) == 0) next

      # If marker_rank exists, use it; otherwise rank by p_val_adj
      if ("marker_rank" %in% colnames(sub)) {
        sub <- sub[order(sub$marker_rank), , drop = FALSE]
      } else if ("p_val_adj" %in% colnames(sub)) {
        sub <- sub[order(sub$p_val_adj), , drop = FALSE]
      } else if ("avg_log2FC" %in% colnames(sub)) {
        sub <- sub[order(-abs(sub$avg_log2FC)), , drop = FALSE]
      }

      # Top N per identity_value
      vals <- unique(sub$identity_value)
      for (v in vals) {
        sub_val <- sub[sub$identity_value == v, , drop = FALSE]
        top_genes <- head(unique(sub_val$gene), pool_top_n)
        gene_set <- union(gene_set, top_genes)
      }
    }
  }

  # Add extra genes
  if (!is.null(extra_genes) && length(extra_genes) > 0) {
    gene_set <- union(gene_set, extra_genes)
  }

  # Cap
  if (length(gene_set) > pool_max_genes) {
    # Keep markers first, then extras; cut from the tail
    gene_set <- head(gene_set, pool_max_genes)
  }

  gene_set
}


# ---- DotPlot Expression Computation ------------------------------------------

#' Precompute DotPlot expression data for all identity layers
#'
#' For each identity layer and each gene in the gene pool, computes:
#'   - avg_expr: mean expression of the gene in cells of that identity value
#'   - pct_expr: fraction of cells with expression > 0
#'   - avg_expr_scaled: min-max scaled avg_expr (per identity_layer, for colour)
#'
#' Also merges marker statistics (avg_log2FC, p_val_adj, marker_rank) when
#' available from the marker data.frame.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_pool Character vector of genes to compute.
#' @param identity_layers Character vector of identity layer names.
#' @param assay Assay name. Default: \code{NULL} (Seurat default assay).
#' @param slot Slot name for expression data. Default: \code{"data"}.
#' @param marker_df Optional normalised marker data.frame for merging stats.
#' @return A data.frame in the canonical dotplot_data format, or NULL on failure.
#' @keywords internal
.compute_dotplot_data <- function(seurat_obj, gene_pool, identity_layers,
                                   assay = NULL, slot = "data",
                                   marker_df = NULL) {
  if (is.null(seurat_obj) || length(gene_pool) == 0 || length(identity_layers) == 0) {
    return(NULL)
  }

  # Get expression matrix
  expr_mat <- tryCatch({
    if (!is.null(assay)) {
      SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = slot)
    } else {
      SeuratObject::GetAssayData(seurat_obj, layer = slot)
    }
  }, error = function(e) {
    warning("Cannot extract expression data: ", e$message)
    NULL
  })

  if (is.null(expr_mat)) return(NULL)

  # Filter to genes present in the matrix
  found_genes <- intersect(gene_pool, rownames(expr_mat))
  skipped <- setdiff(gene_pool, found_genes)

  if (length(found_genes) == 0) {
    warning("None of the gene pool genes are present in the expression matrix")
    return(NULL)
  }

  expr_mat <- expr_mat[found_genes, , drop = FALSE]

  # Build rows
  rows <- list()

  for (lyr in identity_layers) {
    if (!lyr %in% colnames(seurat_obj[[]])) next

    # Get identity values for this layer
    group_vec <- as.character(seurat_obj[[lyr]][, 1])

    # Identify parent info for subcluster fields
    parent_lyr <- NA_character_
    parent_val <- NA_character_
    if (grepl("^subcluster_", lyr)) {
      parent_lyr <- "seurat_clusters"
    }

    for (gene in found_genes) {
      gene_expr <- as.numeric(expr_mat[gene, , drop = TRUE])

      for (grp in sort(unique(group_vec))) {
        cells_in_group <- which(group_vec == grp)
        if (length(cells_in_group) == 0) next

        vals <- gene_expr[cells_in_group]
        avg_expr <- mean(vals)
        pct_expr <- mean(vals > 0)

        # Parent identity value for subclusters
        pv <- if (is.na(parent_lyr)) NA_character_ else grp

        rows[[length(rows) + 1]] <- data.frame(
          identity_layer  = lyr,
          identity_value  = as.character(grp),
          gene            = gene,
          avg_expr        = avg_expr,
          pct_expr        = pct_expr,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) return(NULL)

  dot_data <- do.call(rbind, rows)

  # Scale avg_expr per identity_layer (min-max to [0,1])
  dot_data$avg_expr_scaled <- NA_real_
  for (lyr in unique(dot_data$identity_layer)) {
    idx <- which(dot_data$identity_layer == lyr)
    vals <- dot_data$avg_expr[idx]
    if (length(vals) == 0) next
    min_v <- min(vals, na.rm = TRUE)
    max_v <- max(vals, na.rm = TRUE)
    if (max_v > min_v) {
      dot_data$avg_expr_scaled[idx] <- (vals - min_v) / (max_v - min_v)
    } else {
      dot_data$avg_expr_scaled[idx] <- 0.5
    }
  }

  # Merge marker stats
  # Do NOT pre-create marker_avg_log2FC / marker_p_val_adj / marker_rank on dot_data
  # before merge — that would cause .x / .y suffixes when marker_sub also has them.
  # Instead let merge add them naturally from the renamed marker_sub.
  dot_data$parent_identity_layer <- NA_character_
  dot_data$parent_identity_value <- NA_character_

  if (!is.null(marker_df) && nrow(marker_df) > 0) {
    mark_cols <- intersect(c("gene", "identity_layer", "identity_value",
                              "avg_log2FC", "p_val_adj", "marker_rank"),
                           colnames(marker_df))
    if (all(c("gene", "identity_layer", "identity_value") %in% mark_cols)) {
      marker_sub <- marker_df[, mark_cols, drop = FALSE]

      # Rename marker stats to dot_data contract names before merge
      # so JS sees marker_avg_log2FC / marker_p_val_adj / marker_rank
      if ("avg_log2FC" %in% colnames(marker_sub)) {
        colnames(marker_sub)[colnames(marker_sub) == "avg_log2FC"] <- "marker_avg_log2FC"
      }
      if ("p_val_adj" %in% colnames(marker_sub)) {
        colnames(marker_sub)[colnames(marker_sub) == "p_val_adj"] <- "marker_p_val_adj"
      }

      dot_data <- merge(
        dot_data, marker_sub,
        by = c("gene", "identity_layer", "identity_value"),
        all.x = TRUE, sort = FALSE
      )
    }
  }

  # Ensure canonical marker columns exist after merge (fill NA if missing)
  for (col in c("marker_avg_log2FC", "marker_p_val_adj", "marker_rank")) {
    if (!col %in% colnames(dot_data)) {
      dot_data[[col]] <- if (col == "marker_rank") NA_integer_ else NA_real_
    }
  }

  # Set parent info for subcluster fields
  for (lyr in identity_layers) {
    if (grepl("^subcluster_", lyr)) {
      idx <- which(dot_data$identity_layer == lyr)
      dot_data$parent_identity_layer[idx] <- "seurat_clusters"
    }
  }

  dot_data
}


# ---- Bubble Size Scaling -----------------------------------------------------

#' Scale DotPlot bubble sizes with sqrt compression and adaptive capping
#'
#' Applies sqrt(pct_expr) transformation and maps to a bounded [size_min, size_max]
#' pixel range. If the number of x-axis groups or y-axis genes triggers
#' crowding, the effective max size is reduced.
#'
#' @param pct_expr Numeric vector of pct expressed values (0-1).
#' @param n_x Integer. Number of identity values on x-axis.
#' @param n_y Integer. Number of genes on y-axis.
#' @param size_min Minimum marker size (pixels/Plotly units). Default 3.
#' @param size_max Maximum marker size. Default 14.
#' @param cap_quantile Cap pct_expr at this quantile to avoid outlier bloat. Default 0.95.
#' @return Numeric vector of bounded marker sizes, same length as pct_expr.
#' @keywords internal
.scale_dotplot_bubble_size <- function(pct_expr, n_x, n_y,
                                        size_min = 3, size_max = 14,
                                        cap_quantile = 0.95) {
  if (length(pct_expr) == 0) return(numeric(0))

  # Cap extreme values
  cap_val <- quantile(pct_expr, probs = cap_quantile, na.rm = TRUE)
  pct_capped <- pmin(pct_expr, cap_val)

  # Adaptive max size
  adaptive_max <- size_max
  if (n_x > 35) adaptive_max <- min(adaptive_max, 8)
  else if (n_x > 20) adaptive_max <- min(adaptive_max, 10)
  if (n_y > 60) adaptive_max <- min(adaptive_max, 9)

  # sqrt compression + map to range
  size_min + sqrt(pct_capped) * (adaptive_max - size_min)
}


# ---- Natural Sorting ---------------------------------------------------------

#' Natural sort for identity values
#'
#' Sorts character vectors containing numbers in natural order:
#' "0, 1, 2, 10" not "0, 1, 10, 2".
#'
#' @param x Character vector of identity values.
#' @return Sorted character vector.
#' @keywords internal
.natural_sort_identity_values <- function(x) {
  if (length(x) == 0) return(x)

  # Try numeric sort
  nums <- suppressWarnings(as.numeric(x))
  if (!anyNA(nums)) {
    return(x[order(nums)])
  }

  # Mixed: natural sort via string padding
  # Extract numeric parts, pad, then sort
  x[order(nchar(x), x)]
}


# ---- Gene Panel from Marker DF -----------------------------------------------

#' Build gene panels for each identity layer
#'
#' For a given identity layer, returns:
#'   - auto_panel: top N markers per identity value
#'   - identity_gene_map: which genes belong to which identity value
#'
#' @param marker_df Normalised marker data.frame.
#' @param identity_layer Character. The identity layer to build panels for.
#' @param top_n Integer. Genes per identity value. Default 10.
#' @param direction Character. "up", "down", or "both". Default "up".
#' @param max_genes Integer. Cap total gene count. Default 80.
#' @return A named list with \code{genes} (ordered character vector),
#'   \code{identity_gene_map} (list of identity_value → genes).
#' @keywords internal
.prepare_marker_gene_panels <- function(marker_df, identity_layer,
                                         top_n = 10, direction = "up",
                                         max_genes = 80) {
  if (is.null(marker_df) || nrow(marker_df) == 0) {
    return(list(genes = character(0), identity_gene_map = list()))
  }

  sub <- marker_df[marker_df$identity_layer == identity_layer, , drop = FALSE]
  if (nrow(sub) == 0) {
    return(list(genes = character(0), identity_gene_map = list()))
  }

  # Filter by direction
  if (direction == "up") {
    sub <- sub[sub$avg_log2FC > 0 | is.na(sub$avg_log2FC), , drop = FALSE]
  } else if (direction == "down") {
    sub <- sub[sub$avg_log2FC < 0 | is.na(sub$avg_log2FC), , drop = FALSE]
  }

  # Rank within each identity_value
  if ("marker_rank" %in% colnames(sub)) {
    sub <- sub[order(sub$identity_value, sub$marker_rank), , drop = FALSE]
  } else if ("p_val_adj" %in% colnames(sub)) {
    sub <- sub[order(sub$identity_value, sub$p_val_adj), , drop = FALSE]
  } else if ("avg_log2FC" %in% colnames(sub)) {
    sub <- sub[order(sub$identity_value, -abs(sub$avg_log2FC)), , drop = FALSE]
  }

  # Take top_n per identity_value
  id_map <- list()
  seen_genes <- character(0)

  for (val in unique(sub$identity_value)) {
    val_sub <- sub[sub$identity_value == val, , drop = FALSE]
    val_genes <- head(unique(val_sub$gene), top_n)
    id_map[[as.character(val)]] <- val_genes
  }

  # Build ordered gene list: by identity value, then rank
  all_genes <- character(0)
  for (val in names(id_map)) {
    for (g in id_map[[val]]) {
      if (!g %in% all_genes) {
        all_genes <- c(all_genes, g)
      }
    }
  }

  # Cap
  if (length(all_genes) > max_genes) {
    all_genes <- head(all_genes, max_genes)
  }

  list(genes = all_genes, identity_gene_map = id_map)
}
