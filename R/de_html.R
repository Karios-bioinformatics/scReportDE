# scReportDE: HTML Builder -----------------------------------------------------
#
# CSS, JavaScript, summary cards, section builders, and full-page assembly.
# Visual design mirrors scReportLite / scReportComposition conventions.


# ---- CSS ----------------------------------------------------------------------

report_css <- function() {
'/* === scReportDE v0.1.0 === */

:root {
  --sr-accent: #00b894;
  --sr-accent-dark: #00997a;
  --sr-accent-soft: rgba(0, 184, 148, 0.12);
  --sr-accent-border: rgba(0, 184, 148, 0.35);
  --sr-border: #dfe6e9;
  --sr-text: #2d3436;
  --sr-muted: #95a5a6;
  --sr-radius-sm: 6px;
  --sr-radius-md: 8px;
}

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
               "Helvetica Neue", Arial, sans-serif;
  background: #f5f6fa; color: #2d3436; line-height: 1.5;
  height: 100vh; overflow: hidden;
}

.container {
  height: 100vh; display: flex; flex-direction: column;
}

/* ---- Header ---- */
.report-header {
  display: flex; align-items: center; justify-content: space-between;
  background: #fff; border-bottom: 1px solid var(--sr-border);
  padding: 14px 24px; flex-shrink: 0;
}
.report-title {
  font-size: 1.25em; font-weight: 600; color: var(--sr-text);
}
.report-meta {
  font-size: 0.8em; color: var(--sr-muted);
}

/* ---- Body layout ---- */
.report-body {
  flex: 1; min-height: 0;
  display: block;
  overflow: hidden;
}

/* ---- Top Navigation Bar ---- */
.de-nav {
  flex-shrink: 0;
  background: #fff;
  border-bottom: 1px solid var(--sr-border);
  padding: 0 24px;
  display: flex;
  flex-direction: row;
  align-items: center;
  gap: 6px;
  overflow-x: auto;
  overflow-y: hidden;
  min-height: 44px;
}
.de-nav-label {
  display: none;
}
.de-nav-item {
  flex: 0 0 auto;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-height: 34px;
  padding: 7px 14px;
  cursor: pointer;
  border-left: none;
  border-bottom: 3px solid transparent;
  border-radius: var(--sr-radius-sm) var(--sr-radius-sm) 0 0;
  font-size: 0.82em;
  font-weight: 500;
  color: #636e72;
  white-space: nowrap;
  transition: background 0.15s, color 0.15s, border-color 0.15s;
  user-select: none;
  gap: 6px;
}
.de-nav-item:hover { background: var(--sr-accent-soft); }
.de-nav-item.active {
  background: var(--sr-accent-soft);
  border-bottom-color: var(--sr-accent);
  color: var(--sr-text);
  font-weight: 700;
}
.de-nav-dot {
  display: none;
}

/* ---- Content area ---- */
.de-content {
  min-width: 0; min-height: 0;
  overflow-y: auto; padding: 20px 24px 40px;
}

/* ---- Section ---- */
.de-section {
  display: none;
}
.de-section.de-visible { display: block; }

.de-section-title {
  font-size: 1.05em; font-weight: 600; color: var(--sr-text);
  margin-bottom: 16px; padding-bottom: 8px;
  border-bottom: 1px solid var(--sr-border);
}

/* ---- Summary Cards ---- */
.summary-cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 14px; margin-bottom: 20px;
}
.summary-card {
  background: #fff; border-radius: var(--sr-radius-md);
  padding: 16px 18px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  border-top: 3px solid var(--sr-accent);
}
.summary-card.warn-card {
  border-top-color: #f39c12;
}
.summary-card-value {
  font-size: 1.7em; font-weight: 700; color: var(--sr-text); line-height: 1.1;
}
.summary-card-label {
  font-size: 0.72em; font-weight: 600; color: var(--sr-muted);
  text-transform: uppercase; letter-spacing: 0.5px; margin-top: 4px;
}
.summary-card-detail {
  font-size: 0.78em; color: #636e72; margin-top: 6px; line-height: 1.5;
  word-break: break-all;
}

/* ---- Plots ---- */
.de-plot-block {
  background: #fff; border-radius: var(--sr-radius-md);
  padding: 16px; margin-bottom: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}
.de-plot-block .de-plot-title {
  font-size: 0.78em; font-weight: 600; color: #636e72;
  text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;
}
.de-plot-body { min-height: 350px; }
.de-plot-body > *,
.de-plot-body .html-widget,
.de-plot-body .plotly,
.de-plot-body .js-plotly-plot {
  width: 100% !important;
}

/* ---- DE Table ---- */
.de-table-wrapper {
  background: #fff; border-radius: var(--sr-radius-md);
  padding: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  overflow-x: auto;
}

/* ---- Method Info Table ---- */
.method-table {
  width: 100%; max-width: 600px; border-collapse: collapse;
  font-size: 0.84em; background: #fff;
  border-radius: var(--sr-radius-md);
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  overflow: hidden;
}
.method-table td {
  padding: 8px 14px; border-bottom: 1px solid #f0f1f5;
}
.method-table td:first-child {
  font-weight: 600; color: #636e72; width: 160px;
  background: #f8f9fc;
}
.method-table td:last-child { color: #2d3436; }

/* ---- Warning Banner ---- */
.warning-banner {
  background: #fff3cd; border: 1px solid #ffc107;
  border-radius: var(--sr-radius-sm); padding: 12px 16px;
  margin-bottom: 16px; font-size: 0.82em; color: #856404;
}
.warning-banner .warn-icon {
  font-weight: 700; margin-right: 6px;
}

/* ---- Violin Plot ---- */
#violin-plot {
  min-height: 350px;
  width: 100%;
}
#violin-plot .js-plotly-plot,
#violin-plot .plotly,
#violin-plot .plot-container {
  width: 100% !important;
}

/* ---- No-data placeholder ---- */
.no-data {
  color: var(--sr-muted); font-style: italic; padding: 20px 0;
  text-align: center; background: #fff; border-radius: var(--sr-radius-md);
}

/* ---- Footer ---- */
.report-footer {
  flex-shrink: 0;
  text-align: center; padding: 10px 0;
  font-size: 0.72em; color: var(--sr-muted);
  border-top: 1px solid var(--sr-border);
}

/* ---- Responsive ---- */
@media (max-width: 800px) {
  .report-header {
    flex-direction: column;
    align-items: flex-start;
    gap: 4px;
  }
}
'
}


# ---- JavaScript ---------------------------------------------------------------

report_js <- function() {
'
// === scReportDE v0.1.0 ===

// ── Helpers ────────────────────────────────────────────────

function renderHtmlWidgets() {
  if (window.HTMLWidgets && typeof window.HTMLWidgets.staticRender === "function") {
    try { window.HTMLWidgets.staticRender(); } catch(e) { console.warn(e); }
  }
}

function resizeVisibleWidgets(section) {
  if (!section) return;

  // Plotly resize
  if (window.Plotly) {
    var plots = section.querySelectorAll(".js-plotly-plot");
    plots.forEach(function(p) {
      try { Plotly.Plots.resize(p); } catch(e) {}
    });
  }

  // DT table adjust — hidden-tab init requires explicit redraw
  if (window.jQuery && window.$.fn && window.$.fn.dataTable) {
    var tables = section.querySelectorAll(
      "table.dataTable, .dataTables_wrapper table, table.display"
    );
    tables.forEach(function(tbl) {
      try {
        if ($.fn.dataTable.isDataTable(tbl)) {
          $(tbl).DataTable().columns.adjust().draw(false);
        }
      } catch(e) { console.warn(e); }
    });
  }
}

// ── Section Switching ──────────────────────────────────────

function switchSection(name) {
  var items = document.querySelectorAll(".de-nav-item");
  items.forEach(function(el) { el.classList.remove("active"); });

  var target = document.getElementById("nav-" + name);
  if (target) target.classList.add("active");

  var sections = document.querySelectorAll(".de-section");
  sections.forEach(function(s) { s.classList.remove("de-visible"); });

  var targetSection = document.getElementById("section-" + name);
  if (targetSection) {
    targetSection.classList.add("de-visible");
  }

  // Hidden-tab widgets need a beat to measure, then another to finalise
  setTimeout(function()  { renderHtmlWidgets(); resizeVisibleWidgets(targetSection); }, 100);
  setTimeout(function()  { resizeVisibleWidgets(targetSection); }, 350);
}

// ── Initial Render ─────────────────────────────────────────

window.addEventListener("DOMContentLoaded", function() {
  renderHtmlWidgets();
});

// ── Window Resize ──────────────────────────────────────────

window.addEventListener("resize", function() {
  var visible = document.querySelector(".de-section.de-visible");
  resizeVisibleWidgets(visible);
});
'
}


# ---- Navigation Helpers -------------------------------------------------------

nav_item <- function(id, label) {
  htmltools::tags$div(
    class = "de-nav-item",
    id    = paste0("nav-", id),
    onclick = paste0("switchSection(", shQuote(id), ")"),
    htmltools::tags$span(class = "de-nav-dot"),
    label
  )
}

de_section <- function(id, title, ..., visible = FALSE) {
  htmltools::tags$div(
    class = paste("de-section", if (visible) "de-visible" else ""),
    id = paste0("section-", id),
    htmltools::tags$div(class = "de-section-title", title),
    ...
  )
}

plot_block <- function(title, widget) {
  htmltools::tags$div(
    class = "de-plot-block",
    htmltools::tags$div(class = "de-plot-title", title),
    htmltools::tags$div(class = "de-plot-body",
                        htmltools::as.tags(widget))
  )
}

no_data_block <- function(message = "No data available for this section.") {
  htmltools::tags$div(class = "no-data", message)
}

warning_banner <- function(warnings) {
  if (length(warnings) == 0 || all(nchar(warnings) == 0)) return(NULL)
  items <- lapply(warnings, function(w) {
    htmltools::tags$li(paste0("\u26a0\ufe0f ", w))
  })
  htmltools::tags$div(
    class = "warning-banner",
    htmltools::tags$strong("Warnings:"),
    htmltools::tags$ul(style = "margin:4px 0 0 18px; padding:0;", items)
  )
}


# ---- Summary Card Helpers -----------------------------------------------------

card <- function(label, value, detail = NULL, warn = FALSE) {
  children <- list(
    htmltools::tags$div(class = "summary-card-value", value),
    htmltools::tags$div(class = "summary-card-label", label)
  )
  if (!is.null(detail) && nzchar(detail)) {
    children <- c(children, list(
      htmltools::tags$div(class = "summary-card-detail", detail)
    ))
  }
  card_class <- if (warn) "summary-card warn-card" else "summary-card"
  htmltools::tags$div(class = card_class, children)
}

build_overview_cards <- function(de_df, comparison_label) {
  if (is.null(de_df) || nrow(de_df) == 0) {
    return(list(
      card("DE Genes", "0", "No differentially expressed genes found", warn = TRUE),
      card("Comparison", comparison_label)
    ))
  }

  n_total <- nrow(de_df)

  # Up/down
  if ("avg_log2FC" %in% colnames(de_df)) {
    n_up   <- sum(de_df$avg_log2FC > 0, na.rm = TRUE)
    n_down <- sum(de_df$avg_log2FC < 0, na.rm = TRUE)
  } else {
    n_up   <- NA_integer_
    n_down <- NA_integer_
  }

  # Significant
  if ("p_val_adj" %in% colnames(de_df)) {
    n_sig <- sum(de_df$p_val_adj < 0.05, na.rm = TRUE)
  } else if ("p_val" %in% colnames(de_df)) {
    n_sig <- sum(de_df$p_val < 0.05, na.rm = TRUE)
  } else {
    n_sig <- NA_integer_
  }

  list(
    card("DE Genes", fmt_num(n_total),
         paste(fmt_num(n_up), "up  | ", fmt_num(n_down), "down")),
    card("Up-regulated", if (is.na(n_up)) "N/A" else fmt_num(n_up)),
    card("Down-regulated", if (is.na(n_down)) "N/A" else fmt_num(n_down)),
    card("Significant (p<0.05)", if (is.na(n_sig)) "N/A" else fmt_num(n_sig)),
    card("Comparison", comparison_label)
  )
}


# ---- DE Table Builder ---------------------------------------------------------

#' Build an interactive DT table of DE results
#'
#' @param de_df Normalised DE data.frame
#' @return An htmltools tag (DT widget), or a no-data placeholder
#' @keywords internal
build_de_table <- function(de_df) {
  if (is.null(de_df) || nrow(de_df) == 0) {
    return(no_data_block("No DE results to display."))
  }

  out <- tryCatch({
    # Select and order display columns
    display_cols <- intersect(
      c("gene", "cluster", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2"),
      colnames(de_df)
    )
    if (length(display_cols) == 0) {
      return(no_data_block("No displayable columns in DE results."))
    }

    tbl <- de_df[, display_cols, drop = FALSE]

    # Rename for nicer headers
    nice_names <- c(
      gene        = "Gene",
      cluster     = "Cluster",
      avg_log2FC  = "log2FC",
      p_val       = "p-value",
      p_val_adj   = "adj. p-value",
      pct.1       = "pct.1",
      pct.2       = "pct.2"
    )
    colnames(tbl) <- nice_names[colnames(tbl)]

    # Numeric rounding helper columns (keep original for display formatting)
    cap <- paste0("DE Results (", nrow(tbl), " genes)")

    dt <- DT::datatable(
      tbl,
      caption     = cap,
      rownames    = FALSE,
      filter      = "top",
      extensions  = "Buttons",
      options     = list(
        pageLength   = 25,
        lengthMenu   = c(10, 25, 50, 100),
        scrollX      = TRUE,
        dom          = "Btlipr",
        buttons      = c("copy", "csv"),
        columnDefs   = list(
          list(targets = "_all", className = "dt-right")
        )
      ),
      class       = "display compact"
    )

    # Round numeric columns
    if ("log2FC" %in% colnames(tbl)) {
      dt <- DT::formatRound(dt, "log2FC", digits = 4)
    }
    for (col in c("p-value", "adj. p-value")) {
      if (col %in% colnames(tbl)) {
        # Use formatSignif so tiny p-values show as e.g. 1.23e-120 instead of 0.0000
        dt <- DT::formatSignif(dt, col, digits = 3)
      }
    }
    for (col in c("pct.1", "pct.2")) {
      if (col %in% colnames(tbl)) {
        dt <- DT::formatRound(dt, col, digits = 3)
      }
    }

    # Style p-value columns
    if ("adj. p-value" %in% colnames(tbl)) {
      dt <- DT::formatStyle(
        dt, "adj. p-value",
        color = DT::styleInterval(c(0.01, 0.05), c("#E6194B", "#F58231", "#636e72")),
        fontWeight = DT::styleInterval(c(0.05), c("bold", "normal"))
      )
    }

    htmltools::tags$div(class = "de-table-wrapper", htmltools::as.tags(dt))

  }, error = function(e) {
    warning("DE table generation failed: ", e$message)
    no_data_block(paste("Table error:", e$message))
  })

  out
}


# ---- Method Info Table --------------------------------------------------------

build_method_info <- function(method_meta, generation_time) {
  rows <- lapply(names(method_meta), function(key) {
    htmltools::tags$tr(
      htmltools::tags$td(key),
      htmltools::tags$td(as.character(method_meta[[key]]))
    )
  })

  # Append generation time
  rows <- c(rows, list(
    htmltools::tags$tr(
      htmltools::tags$td("Generated"),
      htmltools::tags$td(generation_time)
    )
  ))

  htmltools::tags$table(class = "method-table", rows)
}


# ---- HTML Assembly -------------------------------------------------------------

#' Build the complete DE HTML report
#'
#' @param de_df Normalised DE data.frame
#' @param overview_cards Summary card HTML
#' @param volcano_widget Volcano plotly widget (or NULL)
#' @param de_table_widget DE table DT widget
#' @param dotplot_widget Dot plotly widget (or NULL)
#' @param violin_widget Violin plotly widget (or NULL)
#' @param method_meta Method info list
#' @param generation_time ISO-8601 timestamp string
#' @param output Output file path
#' @param title Report title
#' @param self_contained Logical (v0.1.0: reserved, not yet implemented)
#' @param warnings Character vector of non-fatal warnings
#' @keywords internal
build_html <- function(de_df, overview_cards, volcano_widget,
                       de_table_widget, dotplot_widget, violin_widget,
                       method_meta, generation_time,
                       output, title, self_contained, warnings) {

  # ---- Build sections ----

  # Overview
  sections <- list()
  sections$overview <- de_section(
    "overview", "Overview", visible = TRUE,
    warning_banner(warnings),
    htmltools::tags$div(class = "summary-cards", overview_cards)
  )

  # Volcano
  sections$volcano <- de_section(
    "volcano", "Volcano Plot",
    if (!is.null(volcano_widget)) plot_block("Differential Expression Volcano", volcano_widget)
    else no_data_block("Volcano plot could not be generated. Check that de_df has logFC and p-value columns.")
  )

  # DE Table
  sections$de_table <- de_section(
    "de_table", "DE Table",
    de_table_widget
  )

  # Dot Plot — interactive panel or static fallback
  if (!is.null(dotplot_widget)) {
    # dotplot_widget is an interactive panel (htmltools tag from .build_interactive_marker_dotplot_panel)
    sections$dotplot <- dotplot_widget
  } else {
    sections$dotplot <- de_section(
      "dotplot", "Dot Plot",
      no_data_block("Dot plot requires a Seurat object and marker data for expression visualisation.")
    )
  }

  # Violin Plot — interactive panel or static fallback
  if (!is.null(violin_widget)) {
    sections$violin <- violin_widget
  } else {
    sections$violin <- de_section(
      "violin", "Violin Plot",
      no_data_block("Violin plot requires a Seurat object for expression data.")
    )
  }

  # Method Info
  sections$method <- de_section(
    "method", "Method Info",
    build_method_info(method_meta, generation_time)
  )

  # ---- Build navigation ----
  nav_items <- list(
    htmltools::tags$div(class = "de-nav-label", "Report"),
    nav_item("overview", "Overview"),
    nav_item("volcano", "Volcano Plot"),
    nav_item("de_table", "DE Table"),
    nav_item("dotplot", "Dot Plot"),
    nav_item("violin", "Violin Plot"),
    nav_item("method", "Method Info")
  )

  # Self-contained notice
  if (!self_contained) {
    notice <- paste(
      "[!] This report is NOT self-contained.",
      "To share it with others, either:",
      "(1) send the HTML file together with the generated 'lib' folder, or",
      "(2) set self_contained = TRUE (not yet implemented in v0.1.0)."
    )
    warning_block <- htmltools::tags$div(
      class = "warning-banner",
      style = "margin:12px 24px 0 24px;",
      htmltools::tags$strong("\u26a0\ufe0f Sharing Notice:"),
      htmltools::tags$span(notice)
    )
    sections$overview <- htmltools::tagAppendChild(
      sections$overview, warning_block
    )
  }

  # ---- Build full page ----
  header <- htmltools::tags$header(class = "report-header",
    htmltools::tags$div(class = "report-title", title),
    htmltools::tags$div(class = "report-meta",
      "Differential Expression Report \u2014 scReportDE v0.1.0")
  )

  top_nav <- htmltools::tags$nav(class = "de-nav", nav_items)
  main    <- htmltools::tags$main(class = "de-content", sections)
  footer  <- htmltools::tags$footer(class = "report-footer",
    "Generated by scReportDE v0.1.0  |  scReport Ecosystem")

  script_tag <- htmltools::tags$script(
    htmltools::HTML(report_js())
  )

  page <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "UTF-8"),
      htmltools::tags$meta(
        name = "viewport",
        content = "width=device-width, initial-scale=1.0"
      ),
      htmltools::tags$title(title),
      htmltools::tags$style(htmltools::HTML(report_css()))
    ),
    htmltools::tags$body(
      htmltools::tags$div(class = "container",
        header,
        top_nav,
        htmltools::tags$div(class = "report-body",
          main
        ),
        footer
      ),
      script_tag
    )
  )

  htmltools::save_html(page, file = output)

  message("DE report written to: ", normalizePath(output, mustWork = FALSE))
  invisible(output)
}
