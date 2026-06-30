# scReportDE v0.1.0

**Differential Expression HTML Report Generator** — part of the scReport ecosystem.

Takes a Seurat object (runs FindMarkers/FindAllMarkers internally) or a pre-computed DE data.frame, and produces an interactive HTML report with volcano plot, DT table, dot plot, and violin plot.

## Installation

```r
# From source
install.packages("path/to/scReportDE", repos = NULL, type = "source")

# Or source directly (for development)
pkg_root <- "path/to/scReportDE"
for (f in list.files(file.path(pkg_root, "R"), full.names = TRUE)) source(f)
```

Requires: `Seurat`, `plotly`, `htmltools`, `jsonlite`, `DT`.

## Quick Start

```r
library(scReportDE)
library(Seurat)

# Option 1: Pairwise DE from Seurat object
obj <- readRDS("my_seurat.rds")
build_screport_de(
  seurat_obj  = obj,
  group_col   = "condition",
  ident_1     = "Treatment",
  ident_2     = "Control",
  mode        = "pairwise",
  output_file = "DE_Treatment_vs_Control.html"
)

# Option 2: Pre-computed DE data.frame
de_results <- FindMarkers(obj, ident.1 = "A", ident.2 = "B")
build_screport_de(de_df = de_results)

# Option 3: AllMarkers mode
build_screport_de(
  seurat_obj  = obj,
  mode        = "all_markers",
  output_file = "DE_AllMarkers.html"
)
```

## Report Sections

| Section | Content |
|---------|---------|
| Overview | Summary cards: n DEGs, up/down/sig counts, comparison, method |
| Volcano Plot | Interactive plotly: logFC vs -log10(p), significant gene highlights |
| DE Table | Searchable, sortable, paginated DT table with full results |
| Dot Plot | Top N genes expression dot plot (requires Seurat object) |
| Violin Plot | Top gene expression distribution across groups (requires Seurat object) |
| Method Info | All parameters and generation metadata |

## Function Reference

### `build_screport_de()`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seurat_obj` | `NULL` | Seurat object (required if no de_df) |
| `de_df` | `NULL` | Pre-computed DE data.frame |
| `group_col` | `NULL` | Metadata column for identity |
| `ident_1` | `NULL` | Identity 1 (pairwise) |
| `ident_2` | `NULL` | Identity 2 (pairwise) |
| `mode` | `"pairwise"` | `"pairwise"` or `"all_markers"` |
| `assay` | `NULL` | Assay name |
| `slot` | `"data"` | Expression slot |
| `test_use` | `"wilcox"` | Statistical test |
| `logfc_threshold` | `0.25` | LogFC threshold |
| `min_pct` | `0.1` | Minimum detection fraction |
| `only_pos` | `FALSE` | Only positive markers |
| `top_n` | `20` | Top genes for labels/dot plot |
| `output_file` | `"scReport_DE.html"` | Output HTML path |
| `title` | `"Differential Expression Report"` | Report title |
| `self_contained` | `FALSE` | Reserved for future use |

### Return value

A list with `de_df`, `output_file`, `metadata`, and `warnings`.

## File Structure

```
scReportDE/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── README.md
├── R/
│   ├── utils.R              # %||%, fmt_num(), fmt_pval(), normalize_de_df(), safe_compute_de()
│   ├── de_compute.R         # compute_de(), build_method_meta()
│   ├── de_plots.R           # plot_volcano(), plot_dotplot(), plot_violin()
│   ├── de_html.R            # CSS, JS, card builders, section builders, build_html()
│   └── build_screport_de.R  # Main API: build_screport_de()
└── inst/
    └── test_de_basic.R      # Full smoke test (7 test cases)
```

## Testing

```bash
cd path/to/scReportDE
Rscript inst/test_de_basic.R
```

7 test cases:
1. Pairwise DE from Seurat
2. Pre-computed de_df only (no Seurat object)
3. Missing optional columns (pct.1/pct.2)
4. Empty DE results (0 genes)
5. AllMarkers mode
6. group_col parameter
7. normalize_de_df edge cases (4 sub-tests)

## Design Notes

- **No crash policy**: Every plot is wrapped in tryCatch. If a plot fails, the report still renders with a placeholder.
- **de_df priority**: If both `seurat_obj` and `de_df` are provided, `de_df` takes precedence.
- **Column auto-detection**: Supports `avg_log2FC`, `avg_logFC`, and `log2FoldChange`; detects gene from column or rownames.
- **Visual style**: Matches scReportLite/scReportComposition (white bg, green accent #00b894, left sidebar, cards).
- **Sharing**: Not self-contained by default. Console message reminds to share the `lib` folder too. `self_contained = TRUE` reserved for future.

## Version History

- **v0.1.0** (2026-06-30): Initial release — pairwise & all_markers DE, volcano plot, DT table, dot plot, violin plot, Method Info.
