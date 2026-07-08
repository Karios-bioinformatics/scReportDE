# scReportDE: Interactive Violin Panel -----------------------------------------
#
# Precomputes top N genes' expression data, embeds as JSON, and provides
# gene/group selectors for on-the-fly violin plot switching without recompilation.


#' Build the interactive Violin expression panel
#'
#' Precomputes expression values for top N marker genes across available
#' identity layers, embeds the data as JSON, and renders selectors so the
#' user can switch genes and grouping variables without regenerating the
#' full report.
#'
#' @param seurat_obj A Seurat object.
#' @param de_df Normalised DE data.frame (for gene ranking).
#' @param group_col Optional primary group column (passed from build_screport_de).
#' @param assay Assay name for expression extraction.
#' @param slot Expression data slot.
#' @param n_genes Number of top genes to precompute (default 10).
#' @return An htmltools tag (the Violin section), or a no-data placeholder.
#' @keywords internal
.build_interactive_violin_panel <- function(
    seurat_obj,
    de_df,
    group_col = NULL,
    assay     = NULL,
    slot      = "data",
    n_genes   = 10) {

  if (is.null(seurat_obj)) {
    return(de_section("violin", "Violin Plot",
      no_data_block("Violin plot requires a Seurat object for expression data.")
    ))
  }

  if (is.null(de_df) || nrow(de_df) == 0) {
    return(de_section("violin", "Violin Plot",
      no_data_block("No DE/marker data available for gene selection.")
    ))
  }

  # ── 1. Pick top N genes by significance ──
  df <- de_df
  if ("p_val_adj" %in% colnames(df)) {
    df <- df[order(df$p_val_adj), , drop = FALSE]
  } else if ("avg_log2FC" %in% colnames(df)) {
    df <- df[order(-abs(df$avg_log2FC)), , drop = FALSE]
  }
  top_genes <- head(unique(df$gene), n_genes)
  if (length(top_genes) == 0) {
    return(de_section("violin", "Violin Plot",
      no_data_block("No valid genes found in DE results.")
    ))
  }

  default_gene <- top_genes[1]

  # ── 2. Collect identity layers for Group by selector ──
  id_layers <- .collect_dotplot_identity_layers(seurat_obj, NULL)
  if (!is.null(group_col) && nzchar(group_col) &&
      group_col %in% colnames(seurat_obj[[]]) &&
      !group_col %in% id_layers) {
    id_layers <- unique(c(group_col, id_layers))
  }
  if (length(id_layers) == 0) {
    id_layers <- "seurat_clusters"
  }
  default_layer <- if (!is.null(group_col) && nzchar(group_col) &&
                       group_col %in% colnames(seurat_obj[[]])) {
    group_col
  } else {
    id_layers[1]
  }

  # ── 3. Precompute expression by gene × layer × group ──
  violin_data <- list()
  for (gene in top_genes) {
    if (!gene %in% rownames(seurat_obj)) next

    gene_expr <- tryCatch({
      if (!is.null(assay)) {
        as.numeric(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = slot)[gene, ])
      } else {
        as.numeric(get_expr_data(seurat_obj, slot)[gene, ])
      }
    }, error = function(e) {
      # fallback to counts
      tryCatch({
        if (!is.null(assay)) {
          as.numeric(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = "counts")[gene, ])
        } else {
          as.numeric(get_expr_data(seurat_obj, "counts")[gene, ])
        }
      }, error = function(e2) NULL)
    })

    if (is.null(gene_expr)) next

    gene_entry <- list()
    for (lyr in id_layers) {
      if (!lyr %in% colnames(seurat_obj[[]])) next
      groups <- as.character(seurat_obj[[lyr]][, 1])
      lyr_entry <- list()
      for (grp in sort(unique(groups))) {
        vals <- gene_expr[groups == grp]
        lyr_entry[[grp]] <- vals
      }
      gene_entry[[lyr]] <- lyr_entry
    }
    violin_data[[gene]] <- gene_entry
  }

  if (length(violin_data) == 0) {
    return(de_section("violin", "Violin Plot",
      no_data_block("Could not extract expression data for any top gene.")
    ))
  }

  # Update top_genes to only those we successfully extracted
  top_genes <- names(violin_data)
  if (!default_gene %in% top_genes) {
    default_gene <- top_genes[1]
  }

  # Get marker stats for subtitle
  gene_info <- list()
  for (g in top_genes) {
    g_row <- de_df[de_df$gene == g, , drop = FALSE]
    if (nrow(g_row) > 0) {
      info <- list()
      if ("avg_log2FC" %in% colnames(g_row))
        info$log2FC <- g_row$avg_log2FC[1]
      if ("p_val_adj" %in% colnames(g_row))
        info$p_adj  <- g_row$p_val_adj[1]
      gene_info[[g]] <- info
    }
  }

  # ── 4. Build JSON ──
  json_data <- list(
    genes          = top_genes,
    default_gene   = default_gene,
    group_layers   = id_layers,
    default_layer  = default_layer,
    gene_info      = gene_info,
    data           = violin_data
  )

  json_str <- jsonlite::toJSON(json_data, auto_unbox = TRUE, na = "null",
                                pretty = FALSE)

  # ── 5. Build embedded JS ──
  violin_js <- paste0('
(function(){
"use strict";
var vData = JSON.parse(document.getElementById("violin-json-data").textContent);
var colors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
              "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78"];

function renderViolin() {
  var gene  = document.getElementById("violin-gene").value;
  var layer = document.getElementById("violin-layer").value;
  var gData = vData.data[gene];
  if (!gData || !gData[layer]) return;

  var groups = Object.keys(gData[layer]).sort(function(a,b) {
    var na=parseFloat(a), nb=parseFloat(b);
    return (!isNaN(na)&&!isNaN(nb)) ? na-nb : String(a).localeCompare(String(b));
  });

  var traces = [];
  groups.forEach(function(grp, i) {
    traces.push({
      y: gData[layer][grp],
      type: "violin",
      name: grp,
      side: "positive",
      box: {visible:true, width:0.2},
      meanline: {visible:true, width:1, color:"#2d3436"},
      points: "all",
      jitter: 0.3,
      pointpos: 0,
      marker: {size:2, opacity:0.4},
      fillcolor: colors[i % colors.length],
      line: {color:"#636e72", width:0.6},
      hoverinfo: "y+name"
    });
  });

  // Update subtitle
  var info = vData.gene_info[gene];
  var subtitle = "Gene: " + gene;
  if (info) {
    var parts = [];
    if (info.log2FC != null) parts.push("log2FC=" + parseFloat(info.log2FC).toFixed(3));
    if (info.p_adj  != null) {
      var p = parseFloat(info.p_adj);
      parts.push("p_adj=" + (p < 1e-4 ? p.toExponential(2) : p.toFixed(4)));
    }
    if (parts.length > 0) subtitle += " (" + parts.join(", ") + ")";
  }
  document.getElementById("violin-subtitle").textContent = subtitle;

  var layout = {
    title: {text: "Expression by Group", font:{size:13}},
    xaxis: {title: layer, gridcolor:"#ecf0f1", zerolinecolor:"#dfe6e9"},
    yaxis: {title: "Expression", gridcolor:"#ecf0f1", zerolinecolor:"#dfe6e9"},
    plot_bgcolor: "#fafbfc",
    paper_bgcolor: "#ffffff",
    margin: {t:40,r:20,b:50,l:60},
    font: {family:"-apple-system,BlinkMacSystemFont,sans-serif", color:"#2d3436"},
    showlegend: false
  };

  Plotly.react("violin-plot", traces, layout, {
    displayModeBar: true,
    modeBarButtonsToRemove: ["lasso2d","select2d"],
    responsive: true
  });
}

document.getElementById("violin-gene").addEventListener("change", renderViolin);
document.getElementById("violin-layer").addEventListener("change", renderViolin);

// Initial render
renderViolin();

})();
')

  # ── 6. Build HTML section ──
  gene_opts <- paste0(
    sprintf('<option value="%s"%s>%s</option>',
            top_genes,
            ifelse(top_genes == default_gene, ' selected', ''),
            top_genes),
    collapse = "\n"
  )

  layer_opts <- paste0(
    sprintf('<option value="%s"%s>%s</option>',
            id_layers,
            ifelse(id_layers == default_layer, ' selected', ''),
            id_layers),
    collapse = "\n"
  )

  section <- htmltools::tags$div(
    class = "de-section",
    id    = "section-violin",

    htmltools::tags$div(class = "de-section-title", "Expression by Group"),

    # Controls
    htmltools::tags$div(
      style = "display:flex;flex-wrap:wrap;gap:14px;align-items:end;margin-bottom:12px;padding:10px 14px;background:#fafbfc;border:1px solid #dfe6e9;border-radius:6px;",

      htmltools::tags$div(
        style = "display:flex;flex-direction:column;gap:3px;",
        htmltools::tags$label("Gene", style = "font-size:0.68em;font-weight:700;color:#636e72;text-transform:uppercase;"),
        htmltools::tags$select(
          id = "violin-gene",
          style = "padding:4px 8px;border:1px solid #dfe6e9;border-radius:4px;font-size:0.82em;",
          htmltools::HTML(gene_opts)
        )
      ),

      htmltools::tags$div(
        style = "display:flex;flex-direction:column;gap:3px;",
        htmltools::tags$label("Group by", style = "font-size:0.68em;font-weight:700;color:#636e72;text-transform:uppercase;"),
        htmltools::tags$select(
          id = "violin-layer",
          style = "padding:4px 8px;border:1px solid #dfe6e9;border-radius:4px;font-size:0.82em;",
          htmltools::HTML(layer_opts)
        )
      )
    ),

    # Subtitle line
    htmltools::tags$div(
      id = "violin-subtitle",
      style = "font-size:0.72em;color:#636e72;margin-bottom:10px;padding:0 2px;"
    ),

    # Plot area
    htmltools::tags$div(
      id = "violin-plot",
      style = "min-height:350px;width:100%;"
    ),

    # Embedded JSON
    htmltools::tags$script(
      id   = "violin-json-data",
      type = "application/json",
      htmltools::HTML(json_str)
    ),

    # Panel JS
    htmltools::tags$script(htmltools::HTML(violin_js))
  )

  section
}
