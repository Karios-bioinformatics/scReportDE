# scReportDE: DotPlot HTML Builder ---------------------------------------------
#
# Assembles the interactive DotPlot panel: reads JS/CSS assets, embeds
# precomputed JSON data, and returns an htmltools tag for the report.


#' Build the interactive marker DotPlot panel
#'
#' Orchestrates the full pipeline: collect identity layers, normalise marker
#' data, build gene pool, precompute expression data, embed as interactive
#' Plotly panel.
#'
#' @param seurat_obj A Seurat object.
#' @param marker_df A DE/marker data.frame (raw, will be normalised).
#' @param identity_layers Character vector of identity layer column names.
#'   When NULL, auto-detected.
#' @param marker_pool_top_n Top N markers per identity value for gene pool.
#' @param pool_max_genes Maximum total genes in the pool.
#' @param top_n Default display top N per group.
#' @param max_display_genes Maximum genes to display on the plot.
#' @param direction Default marker direction: \"up\", \"down\", or \"both\".
#' @param extra_genes Extra genes to force into the gene pool.
#' @param size_min Minimum bubble diameter.
#' @param size_max Maximum bubble diameter.
#' @return An htmltools tag (the DotPlot section), or a no-data placeholder.
#' @keywords internal
.build_interactive_marker_dotplot_panel <- function(
    seurat_obj,
    marker_df           = NULL,
    identity_layers     = NULL,
    group_col           = NULL,
    assay               = NULL,
    slot                = "data",
    marker_pool_top_n   = 50,
    pool_max_genes      = 500,
    top_n               = 10,
    max_display_genes   = 80,
    direction           = "up",
    extra_genes         = NULL,
    size_min            = 3,
    size_max            = 14) {

  # ── 1. Validate input ──
  if (is.null(seurat_obj)) {
    return(de_section("dotplot", "Dot Plot",
      no_data_block("DotPlot requires a Seurat object for expression data.")
    ))
  }

  # ── 2. Collect identity layers ──
  id_layers <- .collect_dotplot_identity_layers(seurat_obj, identity_layers)

  # Inject group_col into identity layers ONLY when it exists in Seurat metadata.
  # A missing meta.data column cannot drive x-axis grouping — but it can still
  # serve as default_identity for marker_df normalisation (see step 3).
  if (!is.null(group_col) && nzchar(group_col)) {
    if (group_col %in% colnames(seurat_obj[[]]) && !group_col %in% id_layers) {
      id_layers <- unique(c(group_col, id_layers))
    }
  }

  if (length(id_layers) == 0) {
    return(de_section("dotplot", "Dot Plot",
      no_data_block("No usable identity layers found for DotPlot. Provide dotplot_identity_layers or ensure meta.data has cluster/annotation fields.")
    ))
  }

  # ── 3. Normalise marker data ──
  # Use group_col as default identity layer only when it exists in Seurat metadata.
  # Otherwise fall back to first auto-detected layer or "seurat_clusters".
  default_id <- if (!is.null(group_col) && nzchar(group_col) && group_col %in% colnames(seurat_obj[[]])) {
    group_col
  } else {
    id_layers[1] %||% "seurat_clusters"
  }
  norm_marker <- .normalize_marker_df(marker_df, default_identity = default_id)

  # ── 4. Build gene pool ──
  gene_pool <- .build_gene_pool(
    marker_df      = norm_marker,
    identity_layers = id_layers,
    extra_genes    = extra_genes,
    pool_top_n     = marker_pool_top_n,
    pool_max_genes = pool_max_genes
  )

  if (length(gene_pool) == 0) {
    # If no marker_df, try extra_genes only
    if (!is.null(extra_genes) && length(extra_genes) > 0) {
      gene_pool <- extra_genes
    } else {
      return(de_section("dotplot", "Dot Plot",
        no_data_block("No marker genes available for DotPlot. Provide a marker data.frame or set dotplot_extra_genes.")
      ))
    }
  }

  # ── 5. Compute dotplot expression data ──
  dot_data <- .compute_dotplot_data(
    seurat_obj       = seurat_obj,
    gene_pool        = gene_pool,
    identity_layers  = id_layers,
    assay            = assay,
    slot             = slot,
    marker_df        = norm_marker
  )

  if (is.null(dot_data) || nrow(dot_data) == 0) {
    return(de_section("dotplot", "Dot Plot",
      no_data_block("Could not compute DotPlot expression data. Check that the gene pool genes exist in the expression matrix.")
    ))
  }

  # ── 6. Build marker metadata for JS ──
  marker_meta <- list()
  for (lyr in id_layers) {
    panel <- .prepare_marker_gene_panels(
      marker_df       = norm_marker,
      identity_layer   = lyr,
      top_n            = max(marker_pool_top_n, top_n * 5),  # generous pool for JS filtering
      direction        = "both",  # store all, let JS filter
      max_genes        = pool_max_genes
    )

    # Also get identity value order
    lyr_data <- dot_data[dot_data$identity_layer == lyr, , drop = FALSE]
    group_order <- .natural_sort_identity_values(unique(lyr_data$identity_value))

    # Build per-group gene lists
    genes_by_group <- list()
    for (g in group_order) {
      if (!is.null(norm_marker) && nrow(norm_marker) > 0) {
        sub <- norm_marker[
          norm_marker$identity_layer == lyr &
          norm_marker$identity_value == g, , drop = FALSE
        ]
        if (nrow(sub) > 0) {
          genes_by_group[[g]] <- head(unique(sub$gene), marker_pool_top_n)
        }
      }
    }

    marker_meta[[lyr]] <- list(
      group_order    = group_order,
      genes_by_group = genes_by_group,
      ranked_genes   = panel$genes,
      up_genes       = if (!is.null(norm_marker) && nrow(norm_marker) > 0) {
        sub_up <- norm_marker[
          norm_marker$identity_layer == lyr &
          (norm_marker$avg_log2FC > 0 | is.na(norm_marker$avg_log2FC)), ,
          drop = FALSE
        ]
        sort(unique(sub_up$gene))
      } else character(0),
      down_genes     = if (!is.null(norm_marker) && nrow(norm_marker) > 0) {
        sub_dn <- norm_marker[
          norm_marker$identity_layer == lyr &
          norm_marker$avg_log2FC < 0, ,
          drop = FALSE
        ]
        sort(unique(sub_dn$gene))
      } else character(0)
    )
  }

  # ── 7. Prepare JSON ──
  json_data <- list(
    identity_layers = id_layers,
    default_layer   = id_layers[1],
    default_top_n   = top_n,
    default_direction = direction,
    default_max_genes = max_display_genes,
    size_min        = size_min,
    size_max        = size_max,
    marker_meta     = marker_meta,
    data            = dot_data
  )

  json_str <- jsonlite::toJSON(json_data, auto_unbox = TRUE, na = "null",
                                dataframe = "rows", pretty = FALSE)

  # ── 8. Read JS / CSS assets ──
  resolve_asset <- function(filename) {
    pkg_path <- system.file(file.path("assets", filename), package = "scReportDE")
    if (nzchar(pkg_path) && file.exists(pkg_path)) return(pkg_path)

    # Dev-mode fallback: look for inst/assets/ relative to package root or cwd
    src_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    bases <- c(getwd(), if (!is.null(src_file)) dirname(dirname(src_file)))
    for (base in bases) {
      dev_path <- file.path(base, "inst", "assets", filename)
      if (file.exists(dev_path)) return(dev_path)
    }
    return("")
  }

  js_path  <- resolve_asset("dotplot_panel.js")
  css_path <- resolve_asset("dotplot_panel.css")

  dotplot_js <- if (nzchar(js_path)) {
    paste(readLines(js_path, warn = FALSE), collapse = "\n")
  } else {
    warning("DotPlot JS asset not found — interactive features will be unavailable")
    ""
  }
  dotplot_css <- if (nzchar(css_path)) {
    paste(readLines(css_path, warn = FALSE), collapse = "\n")
  } else {
    warning("DotPlot CSS asset not found — styling may be incomplete")
    ""
  }

  # ── 9. Build HTML section ──
  section <- htmltools::tags$div(
    class = "de-section",
    id    = "section-dotplot",

    # Inline CSS
    htmltools::tags$style(htmltools::HTML(dotplot_css)),

    htmltools::tags$div(class = "de-section-title", "Marker DotPlot"),

    # Controls
    htmltools::tags$div(id = "dotplot-controls"),

    # Plot area
    htmltools::tags$div(id = "dotplot-plot",
      style = "min-height:350px; width:100%;"),

    # Embedded JSON data
    htmltools::tags$script(
      id   = "dotplot-json-data",
      type = "application/json",
      htmltools::HTML(json_str)
    ),

    # Panel JS
    htmltools::tags$script(htmltools::HTML(dotplot_js))
  )

  section
}
