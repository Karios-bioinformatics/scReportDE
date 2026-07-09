# scReportDE

<p align="center">
  <strong>Interactive differential expression reports for single-cell analysis.</strong>
</p>

<p align="center"><img src="https://img.shields.io/badge/Version-v0.1.0-blue" alt="Version"> <img src="https://img.shields.io/badge/Status-Initial%20release-success" alt="Status"> <img src="https://img.shields.io/badge/Layer-scReport%20Module-lightgrey" alt="Layer"> <img src="https://img.shields.io/badge/Focus-Differential%20Expression-purple" alt="Focus"> <img src="https://img.shields.io/badge/License-MIT-yellow" alt="License"></p>

scReportDE accepts either:

- a Seurat object, then runs `FindMarkers()` / `FindAllMarkers()` internally; or
- a precomputed differential-expression `data.frame`, such as marker output from
  Seurat.

It generates an interactive HTML report with summary cards, a WebGL volcano
plot, a DT result table, an interactive marker DotPlot, an interactive Violin
Plot, and method metadata.

## Installation

```r
# From source
install.packages("path/to/scReportDE", repos = NULL, type = "source")

# Or source directly during development
pkg_root <- "path/to/scReportDE"
for (f in list.files(file.path(pkg_root, "R"), full.names = TRUE)) source(f)
```

Core dependencies: `plotly`, `htmltools`, `jsonlite`, `DT`.

Expression-level plots require a Seurat object and the suggested Seurat
packages: `Seurat`, `SeuratObject`.

## Quick Start

```r
library(scReportDE)
library(Seurat)

# Pairwise DE from a Seurat object
obj <- readRDS("my_seurat.rds")
build_screport_de(
  seurat_obj  = obj,
  group_col   = "condition",
  ident_1     = "Treatment",
  ident_2     = "Control",
  mode        = "pairwise",
  output_file = "DE_Treatment_vs_Control.html"
)

# Precomputed marker table
marker_df <- FindAllMarkers(obj, only.pos = TRUE)
build_screport_de(
  seurat_obj  = obj,
  de_df       = marker_df,
  group_col   = "cluster",
  output_file = "DE_precomputed_markers.html"
)

# All-markers mode
build_screport_de(
  seurat_obj  = obj,
  group_col   = "cluster",
  mode        = "all_markers",
  output_file = "DE_AllMarkers.html"
)
```

When both `seurat_obj` and `de_df` are supplied, the precomputed `de_df` is used
for DE results. The Seurat object is still used, when available, to build
expression-level DotPlot and Violin Plot panels.

## Report Sections

| Section | Content |
|---------|---------|
| Overview | Summary cards for total genes, significant genes, up/down counts, comparison, and method |
| Volcano Plot | Plotly WebGL scatter plot of `log2FC` versus `-log10(adjusted p-value)` |
| DE Table | Searchable, sortable, paginated DT table with full DE results |
| Dot Plot | Interactive marker DotPlot with group, direction, top-N, max-gene, bubble-size, and custom-gene controls |
| Violin Plot | Interactive expression distribution view with gene/group switching |
| Method Info | Input source, mode, assay, slot, thresholds, identities, package version, and generation time |

## Volcano Plot Semantics

For bidirectional pairwise DE results, the volcano plot behaves like a standard
volcano plot:

- x-axis: `log2FC`
- y-axis: `-log10(p_val_adj)` when available, otherwise `-log10(p_val)`
- color: significance and direction
- dashed horizontal line: adjusted p-value threshold

For marker-only or positive-only data, all logFC values may be non-negative.
In that case the plot is intentionally one-sided. It should be read as a marker
ranking view rather than a symmetric up/down contrast. Gene labels are hidden by
default for marker-only data to avoid unreadable label piles; hover text remains
available for every point.

The volcano trace uses Plotly `scattergl` for WebGL rendering on larger marker
tables.

## Function Reference

### `build_screport_de()`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seurat_obj` | `NULL` | Seurat object. Required when `de_df` is not supplied. |
| `de_df` | `NULL` | Precomputed DE or marker table. Takes priority over computed DE. |
| `group_col` | `NULL` | Metadata column used as identity/grouping variable. |
| `ident_1` | `NULL` | Identity 1 for pairwise mode. |
| `ident_2` | `NULL` | Identity 2 for pairwise mode. |
| `mode` | `"pairwise"` | `"pairwise"` or `"all_markers"`. Ignored when `de_df` is supplied. |
| `assay` | `NULL` | Assay passed to Seurat DE and expression plotting. |
| `slot` | `"data"` | Expression slot used for DE and expression plotting. |
| `test_use` | `"wilcox"` | Statistical test passed to Seurat. |
| `logfc_threshold` | `0.25` | LogFC threshold for DE filtering and method metadata. |
| `min_pct` | `0.1` | Minimum detection fraction. |
| `only_pos` | `FALSE` | Return only positive markers when computing with Seurat. |
| `top_n` | `20` | Legacy top-gene setting used by volcano defaults. |
| `volcano_label_top_n` | `NULL` | Override maximum volcano labels. `NULL` auto-selects conservative defaults. |
| `volcano_show_labels` | `NULL` | Override volcano label visibility. `NULL` disables labels for marker-only data. |
| `dotplot_identity_layers` | `NULL` | Metadata columns available as DotPlot grouping layers. |
| `dotplot_marker_pool_top_n` | `50` | Top marker genes per identity included in the DotPlot gene pool. |
| `dotplot_pool_max_genes` | `500` | Maximum precomputed DotPlot gene-pool size. |
| `dotplot_top_n` | `10` | Default number of top markers selected per group. |
| `dotplot_max_display_genes` | `50` | Maximum genes displayed in the DotPlot at once. |
| `dotplot_direction` | `"up"` | Default DotPlot marker direction: `"up"`, `"down"`, or `"both"`. |
| `dotplot_extra_genes` | `NULL` | Extra genes forced into the DotPlot gene pool. |
| `dotplot_size_min` | `3` | Minimum DotPlot bubble diameter. |
| `dotplot_size_max` | `14` | Maximum DotPlot bubble diameter. |
| `output_file` | `"scReport_DE.html"` | Output HTML path. |
| `title` | `"Differential Expression Report"` | Report title. |
| `self_contained` | `FALSE` | Whether to attempt self-contained output. Current reports are normally shared with their dependency folder. |

### Return Value

`build_screport_de()` returns invisibly a list with:

- `de_df`: normalized DE table used by the report
- `output_file`: normalized output path
- `metadata`: method and parameter metadata
- `warnings`: non-fatal warnings encountered during report generation

## Expected Input Columns

`scReportDE` normalizes common DE result schemas. It detects gene identifiers
from a gene column or row names, and detects logFC from common names such as
`avg_log2FC`, `avg_logFC`, or `log2FoldChange`.

The most useful optional columns are:

- adjusted p-value: `p_val_adj`, `padj`, `FDR`, or similar
- raw p-value: `p_val`, `pvalue`, or similar
- detection fractions: `pct.1` and `pct.2`
- marker group: `cluster`, `group`, or another identity-like column

## File Structure

```text
scReportDE/
|-- DESCRIPTION
|-- NAMESPACE
|-- LICENSE
|-- README.md
|-- R/
|   |-- utils.R
|   |-- de_compute.R
|   |-- de_plots.R
|   |-- de_html.R
|   |-- dotplot_compute.R
|   |-- dotplot_html.R
|   |-- violin_html.R
|   `-- build_screport_de.R
`-- inst/
    |-- assets/
    |   |-- dotplot_panel.css
    |   `-- dotplot_panel.js
    |-- test_de_basic.R
    `-- test_de_html_contract.R
```

## Testing

```bash
cd path/to/scReportDE
Rscript inst/test_de_basic.R
Rscript inst/test_de_html_contract.R
```

The current test coverage is focused on smoke tests, precomputed `de_df` paths,
missing optional columns, empty DE results, all-markers mode, `group_col`
handling, normalization edge cases, and generated HTML contract checks.

## Design Notes

- **No-crash report generation**: plot sections are guarded so one failed plot
  should not prevent the report from being written.
- **Precomputed DE priority**: supplied `de_df` is treated as the source of
  truth for DE results.
- **Expression panels need Seurat**: DotPlot and Violin Plot need `seurat_obj`
  because they read expression matrices and metadata.
- **Hidden-tab widgets**: the report JavaScript explicitly wakes/resizes
  htmlwidgets after tab switches so DT and Plotly sections render reliably.
- **DotPlot overflow control**: DotPlot is designed around a capped display
  gene count and an internal scrollable plotting area for large marker panels.
- **Sharing**: reports are not self-contained by default. Keep the generated
  HTML file together with its dependency folder when moving or sharing output.

## Version History

- **v0.1.0**: Initial release with pairwise/all-marker DE support,
  precomputed `de_df` input, WebGL volcano rendering, marker-only volcano
  behavior, robust DT result table rendering, interactive DotPlot controls,
  interactive Violin Plot, Method Info, hidden-tab widget rendering fixes, and
  HTML contract checks.
