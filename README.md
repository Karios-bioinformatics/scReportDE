# scReportDE

<p align="center"> <strong>A differential expression reporting module for the scReport ecosystem.</strong> </p>

<p align="center">
 <img src="https://img.shields.io/badge/Version-v0.0.0--alpha-blue" alt="Version">
 <img src="https://img.shields.io/badge/Status-Early%20Alpha-orange" alt="Status">
 <img src="https://img.shields.io/badge/Layer-scReport%20Module-lightgrey" alt="Layer">
 <img src="https://img.shields.io/badge/Focus-Differential%20Expression-purple" alt="Focus"> 
 <img src="https://img.shields.io/badge/License-MIT-yellow" alt="License">
</p>

## Overview

**scReportDE** is a planned differential expression reporting module in the **scReport** ecosystem.

It is designed to organize, visualize, and communicate differential expression results from single-cell bioinformatics workflows.

The core idea is:

> Convert pre-computed differential expression results into structured, interactive, and shareable HTML reports.

`scReportDE` does not aim to replace upstream differential expression tools. Instead, it focuses on the reporting layer after differential expression analysis has already been performed.

## Current status

The current release is an early alpha project definition.

`v0.0.0-alpha` establishes the repository, project scope, README, and development roadmap. It does not yet provide a functional R package implementation.

The first functional release is planned as `v0.1.0`.

## Position in the scReport ecosystem

`scReportDE` is one module of the broader `scReport` ecosystem.

In the scReport design:

- `scReportLite` focuses on cell-level views such as QC, Feature, PCA, UMAP, marker linkage, and selected gene expression.
- `scReportComposition` focuses on sample-level and group-level cell composition.
- `scReportDE` will focus on differential expression result reporting.
- Future modules may cover enrichment analysis, trajectory analysis, cell communication, regulatory networks, spatial omics, and multi-omics reporting.

`scReportDE` is intended to answer questions such as:

- Which genes are differentially expressed between groups?
- Which genes are upregulated or downregulated in a cluster or cell type?
- Which comparisons show the strongest transcriptional changes?
- Which differential expression results should be linked to UMAP, composition, and enrichment reports?
- How can DE tables, volcano plots, MA plots, and gene-level summaries be organized into a single report?

## Core concept

The central object of `scReportDE` is the differential expression result table.

A typical DE table may contain:

|comparison|cluster|gene|avg_log2FC|p_val|p_val_adj|pct.1|pct.2|
|---|---|---|---|---|---|---|---|
|Disease_vs_Control|Macrophage|IL1B|1.82|1e-10|2e-8|0.71|0.24|
|Disease_vs_Control|T cell|GZMB|1.15|3e-6|1e-4|0.42|0.18|
|Disease_vs_Control|Epithelial|KRT8|-0.93|2e-5|6e-4|0.21|0.55|

The exact required columns may vary depending on the upstream tool, but `gene`, effect size, adjusted p-value, and comparison labels are expected to be central fields.

## Planned v0.1.0 scope

The first functional version of `scReportDE` will focus on **pre-computed differential expression result reporting**.

Planned features for `v0.1.0`:

- Accept pre-computed DE result tables.
- Summarize DE results by comparison, cluster, or cell type.
- Generate a single-page HTML DE report.
- Provide interactive DE result tables.
- Provide volcano plots.
- Provide MA plots.
- Provide top upregulated and downregulated gene summaries.
- Provide comparison-level summary cards.
- Support common Seurat-style marker / DE result columns.

## Out of scope for v0.1.0

The following features are intentionally not included in the first functional version:

- Running differential expression tests internally.
- Replacing Seurat, Scanpy, MAST, edgeR, DESeq2, limma, or other upstream tools.
- Pseudobulk modeling.
- Multi-factor statistical design.
- Pathway enrichment analysis.
- Gene set enrichment analysis.
- Cell-cell communication analysis.
- Trajectory-specific DE analysis.
- Cross-module cell locking.

These features may be supported by future scReport modules or later versions.

## Input data

`scReportDE` is designed to work with differential expression result tables such as:

|   |   |   |   |   |   |   |
|---|---|---|---|---|---|---|
|comparison|group_1|group_2|cluster|gene|avg_log2FC|p_val_adj|
|Disease_vs_Control|Disease|Control|Macrophage|IL1B|1.82|2e-8|
|Disease_vs_Control|Disease|Control|Macrophage|CXCL8|1.44|5e-6|
|Disease_vs_Control|Disease|Control|T cell|GZMB|1.15|1e-4|

Minimum expected columns:

- `gene`
- effect size column, such as `avg_log2FC`, `log2FC`, or `logFC`
- adjusted p-value column, such as `p_val_adj`, `padj`, or `FDR`

Recommended columns:

- `comparison`
- `cluster`
- `cell_type`
- `group_1`
- `group_2`
- `pct.1`
- `pct.2`
- `p_val`

## Planned usage

The planned core workflow is:

```
library(scReportDE)

build_de_report(
  de_df = de_df,
  gene_col = "gene",
  logfc_col = "avg_log2FC",
  padj_col = "p_val_adj",
  comparison_col = "comparison",
  cluster_col = "cluster",
  output = "scReportDE.html"
)
```

Or with a Seurat-style marker table:

```
build_de_report(
  de_df = marker_df,
  gene_col = "gene",
  logfc_col = "avg_log2FC",
  padj_col = "p_val_adj",
  cluster_col = "cluster",
  output = "scReportDE_marker_report.html"
)
```

## Example report panels

A typical `scReportDE` report may include:

1. Overview summary cards
2. Comparison selector
3. Cluster / cell type selector
4. Volcano plot
5. MA plot
6. Top upregulated genes
7. Top downregulated genes
8. Interactive DE result table
9. Optional gene-level detail panel
10. Optional export-ready summary table

## Design principles

`scReportDE` follows the general design principles of the scReport ecosystem:

- Reporting after analysis, not replacing analysis.
- Use pre-computed upstream results.
- Keep statistical modeling separate from report generation.
- Provide clear interactive visualization and result organization.
- Preserve compatibility with future scReport modules.
- Maintain stable identifiers for genes, clusters, cell types, samples, and comparisons.

## Relationship to other scReport modules

`scReportDE` is designed to connect naturally with other modules in the scReport ecosystem.

```
scReportLite
  → Where are the cells?
  → What clusters or cell types are being inspected?

scReportComposition
  → Which cell types change in proportion across samples or groups?

scReportDE
  → Which genes change within a cluster, cell type, or comparison?

scReportEnrichment
  → What biological functions are associated with the changed genes?
```

In the future, DE results may be linked back to cell-level and group-level context through shared fields such as:

- `cell_id`
- `cluster`
- `cell_type`
- `sample`
- `group`
- `comparison`
- `gene`

## Cell-centric design direction

The long-term scReport ecosystem is intended to support **cell-centric global tracking** across modules.

A guiding principle is:

> A cell should not lose its identity when the user moves across analysis modules.

For `scReportDE`, this means that DE results should be traceable back to the relevant cluster, cell type, group comparison, marker table, UMAP context, and downstream enrichment interpretation.

## Roadmap

### v0.0.0-alpha

Repository initialization, project definition, README, and roadmap.

### v0.1.0

Pre-computed DE result reporting.

Potential additions:

- Interactive DE table.
- Volcano plot.
- MA plot.
- Top upregulated / downregulated gene summaries.
- Basic comparison and cluster selectors.

### v0.2.0

Improved DE report structure.

Potential additions:

- Multi-comparison support.
- Cluster / cell type faceting.
- Gene search and detail panels.
- Export-ready filtered result tables.

### v0.3.0

Integration with enrichment reporting.

Potential additions:

- Passing selected DE gene sets to `scReportEnrichment`.
- Upregulated and downregulated gene set export.
- Linked DE-to-enrichment report flow.

### Future

Integration with the main `scReport` ecosystem.

Potential directions:

- Cross-module linking with UMAP, composition, and enrichment modules.
- Cell-centric global tracking.
- Pseudobulk DE result display.
- More flexible support for upstream DE tools.

## Citation

A Zenodo DOI will be added after the first archived release.

Current DOI:

```
To be added after release.
```

## License

This project is released under the MIT License.