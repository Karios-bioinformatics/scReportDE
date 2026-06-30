# scReportDE: Plot Functions ---------------------------------------------------
#
# Each plot function returns a plotly widget (or NULL on failure).
# All are wrapped in tryCatch so a single plot failure never kills the report.


#' Define a shared plotly theme for DE reports
#'
#' Returns a list of layout options applied to every plot.
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @return Named list for plotly::layout()
#' @keywords internal
de_plot_theme <- function(title = "", xlab = "", ylab = "") {
  list(
    title    = list(text = title, font = list(size = 13)),
    xaxis    = list(title = xlab, gridcolor = "#ecf0f1", zerolinecolor = "#dfe6e9"),
    yaxis    = list(title = ylab, gridcolor = "#ecf0f1", zerolinecolor = "#dfe6e9"),
    plot_bgcolor  = "#fafbfc",
    paper_bgcolor = "#ffffff",
    margin   = list(t = 40, r = 20, b = 50, l = 60),
    font     = list(family = "-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif",
                    color = "#2d3436"),
    hovermode = "closest"
  )
}


#' Volcano plot from DE results
#'
#' Generates an interactive plotly volcano plot with log2FC on the x-axis
#' and -log10(p-value) on the y-axis. Significant genes (p_val_adj < alpha)
#' are highlighted. Top N up- and down-regulated genes by absolute logFC
#' are labelled.
#'
#' @param de_df Normalised DE data.frame (must contain gene, avg_log2FC, p_val_adj)
#' @param top_n Number of top genes to label (default 20)
#' @param alpha Significance threshold (default 0.05)
#' @param comparison_label Display label (e.g. "A vs B")
#' @return A plotly widget, or NULL on failure
#' @keywords internal
plot_volcano <- function(de_df, top_n = 20, alpha = 0.05,
                         comparison_label = NULL) {
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)

  out <- tryCatch({
    df <- de_df

    # Ensure required columns exist
    if (!"avg_log2FC" %in% colnames(df)) return(NULL)
    y_col <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else "p_val"
    if (!y_col %in% colnames(df)) return(NULL)

    # Cap extreme values for cleaner display
    logfc <- pmax(pmin(df$avg_log2FC, 5), -5)

    # Compute -log10; clip at 20 for visual sanity
    pval <- df[[y_col]]
    pval[pval <= 0] <- .Machine$double.xmin
    neglog10 <- pmin(-log10(pval), 20)

    # Significance classification
    is_sig  <- if ("p_val_adj" %in% colnames(df)) df$p_val_adj < alpha else df$p_val < alpha
    sig_cat <- ifelse(is_sig,
      ifelse(df$avg_log2FC > 0, "Up (sig)", "Down (sig)"),
      "Not sig"
    )
    colors <- c("Up (sig)" = "#E6194B", "Down (sig)" = "#3CB44B",
                "Not sig" = "#b2bec3")

    # Hover text
    hover <- paste0(
      "<b>", df$gene, "</b><br>",
      "log2FC: ", sprintf("%.4f", df$avg_log2FC), "<br>",
      "p_val: ", fmt_pval(df$p_val), "<br>"
    )
    if ("p_val_adj" %in% colnames(df)) {
      hover <- paste0(hover, "p_adj: ", fmt_pval(df$p_val_adj), "<br>")
    }
    if ("pct.1" %in% colnames(df) && !all(is.na(df$pct.1))) {
      hover <- paste0(hover,
        "pct.1: ", sprintf("%.3f", df$pct.1),
        " / pct.2: ", sprintf("%.3f", df$pct.2))
    }

    # Top N labels (by abs logFC among significant genes)
    label_df <- df[is_sig, , drop = FALSE]
    if (nrow(label_df) > 0) {
      label_df <- label_df[order(-abs(label_df$avg_log2FC)), , drop = FALSE]
      label_df <- head(label_df, top_n)
    }

    title_text <- "Volcano Plot"
    if (!is.null(comparison_label)) {
      title_text <- paste0(title_text, ": ", comparison_label)
    }

    p <- plotly::plot_ly(
      data      = df,
      x         = ~logfc,
      y         = ~neglog10,
      type      = "scatter",
      mode      = "markers",
      color     = ~sig_cat,
      colors    = colors,
      text      = ~hover,
      hoverinfo = "text",
      marker    = list(
        size       = 4,
        opacity    = 0.7,
        line       = list(width = 0.3, color = "#ffffff")
      )
    )

    p <- do.call(plotly::layout, c(list(p), de_plot_theme(
      title = title_text,
      xlab  = "log2 Fold Change",
      ylab  = paste0("-log10(", y_col, ")")
    )))

    # Add significance threshold line
    if (y_col == "p_val_adj" || y_col == "p_val") {
      thresh <- -log10(alpha)
      p <- plotly::add_segments(p, x = -10, xend = 10, y = thresh, yend = thresh,
                                line = list(dash = "dash", color = "#636e72",
                                            width = 0.6),
                                inherit = FALSE, showlegend = FALSE)
    }

    # Add gene labels for top N
    if (nrow(label_df) > 0) {
      p <- plotly::add_annotations(
        p,
        x      = pmax(pmin(label_df$avg_log2FC, 5), -5),
        y      = pmin(-log10(pmax(label_df[[y_col]], .Machine$double.xmin)), 20),
        text   = label_df$gene,
        showarrow = TRUE,
        arrowhead = 1,
        arrowsize = 0.4,
        arrowwidth = 0.5,
        font   = list(size = 7, color = "#2d3436"),
        xanchor = "center",
        yanchor = "bottom"
      )
    }

    plotly::config(p, displayModeBar = TRUE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d"))

  }, error = function(e) {
    warning("Volcano plot generation failed: ", e$message)
    NULL
  })

  out
}


#' Dot plot of top DE genes
#'
#' Produces a plotly dot plot where rows are the top \code{n_genes}
#' significant genes and columns are identity groups. Dot size reflects
#' the fraction of cells expressing the gene (pct.exp) and dot colour
#' reflects the average expression.
#'
#' If \code{seurat_obj} is NULL, returns NULL (no data placeholder handled
#' by caller).
#'
#' @param seurat_obj A Seurat object (optional)
#' @param de_df Normalised DE data.frame
#' @param top_n Number of top genes (default 20)
#' @param n_genes Display at most this many genes (default 15)
#' @param group_col Optional group column for x-axis groups
#' @return A plotly widget, or NULL
#' @keywords internal
plot_dotplot <- function(seurat_obj, de_df, top_n = 20, n_genes = 15,
                         group_col = NULL) {
  if (is.null(seurat_obj)) return(NULL)
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)

  out <- tryCatch({
    # Select top genes by significance + effect size
    df <- de_df
    if (!"p_val_adj" %in% colnames(df) || !"avg_log2FC" %in% colnames(df)) {
      return(NULL)
    }

    df <- df[order(df$p_val_adj, -abs(df$avg_log2FC)), , drop = FALSE]
    top_genes <- head(df$gene, min(n_genes, nrow(df)))

    # Determine groups
    if (!is.null(group_col) && group_col %in% colnames(seurat_obj[[]])) {
      groups <- sort(unique(seurat_obj[[group_col]][, 1]))
    } else {
      groups <- sort(unique(Seurat::Idents(seurat_obj)))
    }
    groups <- as.character(groups)

    # Compute expression matrix for top genes
    expr_mat <- tryCatch({
      Seurat::GetAssayData(seurat_obj, slot = "data")[top_genes, , drop = FALSE]
    }, error = function(e) {
      message("Cannot extract expression; trying counts slot...")
      tryCatch({
        Seurat::GetAssayData(seurat_obj, slot = "counts")[top_genes, , drop = FALSE]
      }, error = function(e2) NULL)
    })

    if (is.null(expr_mat)) return(NULL)

    # Ensure all top_genes are in the matrix
    found_genes <- intersect(top_genes, rownames(expr_mat))
    if (length(found_genes) == 0) return(NULL)
    expr_mat <- expr_mat[found_genes, , drop = FALSE]

    # Compute pct.exp and avg.exp per group
    build_dot_data <- function() {
      meta_groups <- if (!is.null(group_col) && group_col %in% colnames(seurat_obj[[]])) {
        as.character(seurat_obj[[group_col]][, 1])
      } else {
        as.character(Seurat::Idents(seurat_obj))
      }

      rows <- list()
      for (g in found_genes) {
        for (grp in groups) {
          cells_in_group <- which(meta_groups == grp)
          if (length(cells_in_group) == 0) next
          expr_vals <- as.numeric(expr_mat[g, cells_in_group, drop = TRUE])
          pct_exp <- mean(expr_vals > 0)
          avg_exp <- mean(expr_vals)
          rows[[length(rows) + 1]] <- data.frame(
            gene    = g,
            group   = grp,
            pct_exp = pct_exp,
            avg_exp = avg_exp,
            stringsAsFactors = FALSE
          )
        }
      }
      do.call(rbind, rows)
    }

    dot_data <- build_dot_data()
    if (is.null(dot_data) || nrow(dot_data) == 0) return(NULL)

    # Scale avg_exp for colour mapping
    max_exp <- max(dot_data$avg_exp)
    min_exp <- min(dot_data$avg_exp)
    if (max_exp == min_exp) max_exp <- min_exp + 1

    p <- plotly::plot_ly(
      data        = dot_data,
      x           = ~group,
      y           = ~gene,
      type        = "scatter",
      mode        = "markers",
      marker      = list(
        size        = ~pct_exp * 20,
        sizemode    = "diameter",
        sizeref     = 0.2,
        sizemin     = 2,
        color       = ~avg_exp,
        colorscale  = list(c(0, "#dcdde1"), c(1, "#E6194B")),
        cmin        = min_exp,
        cmax        = max_exp,
        showscale   = TRUE,
        colorbar    = list(title = "Avg expr", len = 0.5),
        line        = list(width = 0.3, color = "#ffffff")
      ),
      text        = ~paste0(
        "<b>", gene, "</b> | ", group, "<br>",
        "Avg expr: ", sprintf("%.3f", avg_exp), "<br>",
        "Pct exp: ", sprintf("%.1f%%", pct_exp * 100)
      ),
      hoverinfo   = "text"
    )

    p <- do.call(plotly::layout, c(list(p), de_plot_theme(
      title = "Top DE Genes Dot Plot",
      xlab  = if (!is.null(group_col)) group_col else "Identity",
      ylab  = "Gene"
    )))

    p <- plotly::layout(p,
      xaxis = list(tickangle = -45),
      margin = list(l = 100, r = 20, t = 40, b = 80)
    )

    plotly::config(p, displayModeBar = TRUE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d"))

  }, error = function(e) {
    warning("Dot plot generation failed: ", e$message)
    NULL
  })

  out
}


#' Violin plot of a single top DE gene
#'
#' Plots expression distribution of the top DE gene (by significance)
#' across identity groups. If \code{seurat_obj} is NULL, returns NULL.
#'
#' @param seurat_obj A Seurat object (optional)
#' @param de_df Normalised DE data.frame
#' @param gene Optional specific gene to plot; defaults to top gene by p_val_adj
#' @param group_col Optional group column for x-axis
#' @return A plotly widget, or NULL
#' @keywords internal
plot_violin <- function(seurat_obj, de_df, gene = NULL, group_col = NULL) {
  if (is.null(seurat_obj)) return(NULL)
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)

  out <- tryCatch({
    # Pick top gene
    if (is.null(gene)) {
      if ("p_val_adj" %in% colnames(de_df)) {
        gene <- de_df$gene[which.min(de_df$p_val_adj)]
      } else {
        gene <- de_df$gene[1]
      }
    }

    if (is.null(gene) || is.na(gene) || !gene %in% rownames(seurat_obj)) {
      return(NULL)
    }

    # Extract expression
    expr_vec <- tryCatch({
      as.numeric(Seurat::GetAssayData(seurat_obj, slot = "data")[gene, ])
    }, error = function(e) {
      tryCatch({
        as.numeric(Seurat::GetAssayData(seurat_obj, slot = "counts")[gene, ])
      }, error = function(e2) NULL)
    })

    if (is.null(expr_vec)) return(NULL)

    # Determine groups
    if (!is.null(group_col) && group_col %in% colnames(seurat_obj[[]])) {
      groups <- as.character(seurat_obj[[group_col]][, 1])
    } else {
      groups <- as.character(Seurat::Idents(seurat_obj))
    }

    uniq_groups <- sort(unique(groups))
    n_groups    <- length(uniq_groups)
    group_colors <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
      "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
    )
    group_colors <- setNames(rep(group_colors, length.out = n_groups), uniq_groups)

    # Build violin traces per group
    p <- plotly::plot_ly()
    for (i in seq_along(uniq_groups)) {
      grp  <- uniq_groups[i]
      vals <- expr_vec[groups == grp]
      if (length(vals) == 0) next

      p <- plotly::add_trace(
        p,
        y       = ~vals,
        type    = "violin",
        name    = grp,
        side    = "positive",
        box     = list(visible = TRUE, width = 0.2),
        meanline = list(visible = TRUE, width = 1, color = "#2d3436"),
        points  = "all",
        jitter  = 0.3,
        pointpos = 0,
        marker  = list(size = 2, opacity = 0.4),
        fillcolor = group_colors[grp],
        line    = list(color = "#636e72", width = 0.6),
        hoverinfo = "y+name"
      )
    }

    p <- do.call(plotly::layout, c(list(p), de_plot_theme(
      title = paste0("Expression of ", gene),
      xlab  = if (!is.null(group_col)) group_col else "Identity",
      ylab  = "Expression"
    )))

    p <- plotly::layout(p, showlegend = FALSE)

    plotly::config(p, displayModeBar = TRUE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d"))

  }, error = function(e) {
    warning("Violin plot generation failed: ", e$message)
    NULL
  })

  out
}
