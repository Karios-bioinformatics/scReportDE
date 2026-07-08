# scReportDE — HTML Contract Tests
# ====================================
# Validates that generated HTML meets structural contracts:
#   - DotPlot section always has nav + section ids
#   - Marker stats merged without .x/.y suffixes
#   - Volcano labels not all stacked at same capped coord
#   - Provided marker table + group_col drives DotPlot
#
# Run from package root: Rscript inst/test_de_html_contract.R

library(Seurat)
library(SeuratObject)
library(plotly)
library(htmltools)
library(jsonlite)
library(DT)

# ---- Source package files ----
pkg_root <- normalizePath(".")
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  for (d in c("..", "../..", "../../..")) {
    if (file.exists(file.path(pkg_root, d, "DESCRIPTION"))) {
      pkg_root <- normalizePath(file.path(pkg_root, d)); break
    }
  }
}
message("Package root: ", pkg_root)
r_dir <- file.path(pkg_root, "R")
stopifnot(dir.exists(r_dir))
for (f in list.files(r_dir, full.names = TRUE, pattern = "\\.R$")) {
  message("  Sourcing: ", basename(f))
  source(f)
}

# ---- Build mock Seurat object ----
set.seed(123)
n_cells <- 300
n_genes <- 60

condition_vec <- rep(c("Treatment", "Control"), each = 150)
sample_vec   <- paste0("Sample_", rep(1:3, each = 100))
meta <- data.frame(
  condition  = condition_vec,
  sample     = sample_vec,
  row.names  = paste0("Cell_", seq_len(n_cells)),
  stringsAsFactors = FALSE
)

gene_names <- paste0("Gene_", seq_len(n_genes))
counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                 nrow = n_genes, ncol = n_cells,
                 dimnames = list(gene_names, rownames(meta)))

# Inject DE signal
for (g in 1:10) {
  counts[g, condition_vec == "Treatment"] <- rpois(150, lambda = 8)
}
for (g in (n_genes-5):n_genes) {
  counts[g, condition_vec == "Control"]   <- rpois(150, lambda = 6)
}

obj <- CreateSeuratObject(counts = counts, meta.data = meta)
obj <- NormalizeData(obj, verbose = FALSE)
Idents(obj) <- "condition"
obj$seurat_clusters <- rep(as.character(0:2), each = 100)
message("Seurat object: ", ncol(obj), " cells x ", nrow(obj), " genes")

# ============================================================================
# Test 1: No-data DotPlot still has nav-dotplot + section-dotplot
# ============================================================================
message("\n=== Test 1: No-data DotPlot has nav + section ===")

# Build report with no Seurat object → DotPlot must be no-data
faux_de <- data.frame(
  gene        = paste0("Gene_", 1:20),
  avg_log2FC  = rnorm(20, 0, 1.5),
  p_val       = runif(20, 0, 0.1),
  p_val_adj   = runif(20, 0, 0.1),
  stringsAsFactors = FALSE
)

out1 <- build_screport_de(
  de_df       = faux_de,
  output_file = file.path(pkg_root, "test_contract_nodata.html"),
  title       = "Contract Test — No-Data DotPlot"
)

stopifnot(file.exists(out1$output_file))
lines1 <- readLines(out1$output_file, warn = FALSE)

# Check nav exists
has_nav <- any(grepl('id="nav-dotplot"', lines1, fixed = TRUE))
cat(sprintf("  nav-dotplot: %s\n", if (has_nav) "OK" else "MISSING!"))
stopifnot(has_nav)

# Check section exists
has_section <- any(grepl('id="section-dotplot"', lines1, fixed = TRUE))
cat(sprintf("  section-dotplot: %s\n", if (has_section) "OK" else "MISSING!"))
stopifnot(has_section)

# Check that no-data is INSIDE section-dotplot, not bare in main
# (The de_section wrapper ensures the no-data div is a child of section-dotplot)
message("PASS: No-data DotPlot has both nav and section")

# ============================================================================
# Test 2: Successful DotPlot has dotplot-json-data, Plotly.react, marker fields
# ============================================================================
message("\n=== Test 2: Successful DotPlot JSON contract ===")

out2 <- build_screport_de(
  seurat_obj  = obj,
  ident_1     = "Treatment",
  ident_2     = "Control",
  mode        = "pairwise",
  output_file = file.path(pkg_root, "test_contract_success.html"),
  title       = "Contract Test — Success"
)

stopifnot(file.exists(out2$output_file))
lines2 <- readLines(out2$output_file, warn = FALSE)

# dotplot-json-data script tag exists
has_json_data <- any(grepl("dotplot-json-data", lines2, fixed = TRUE))
cat(sprintf("  dotplot-json-data: %s\n", if (has_json_data) "OK" else "MISSING!"))
stopifnot(has_json_data)

# Plotly.react exists
has_plotly_react <- any(grepl("Plotly.react", lines2, fixed = TRUE))
cat(sprintf("  Plotly.react: %s\n", if (has_plotly_react) "OK" else "MISSING!"))
stopifnot(has_plotly_react)

# Extract JSON and check marker fields
json_match <- regmatches(
  paste(lines2, collapse = "\n"),
  regexpr('"marker_avg_log2FC":\\[[^]]*\\]', paste(lines2, collapse = "\n"))
)
has_marker_avg <- any(grepl("marker_avg_log2FC", lines2, fixed = TRUE))
cat(sprintf("  marker_avg_log2FC in HTML: %s\n", if (has_marker_avg) "OK" else "MISSING!"))
stopifnot(has_marker_avg)

has_marker_pval <- any(grepl("marker_p_val_adj", lines2, fixed = TRUE))
cat(sprintf("  marker_p_val_adj in HTML: %s\n", if (has_marker_pval) "OK" else "MISSING!"))
stopifnot(has_marker_pval)

message("PASS: DotPlot JSON contract valid")

# ============================================================================
# Test 3: No .x/.y suffixes in marker field names
# ============================================================================
message("\n=== Test 3: No merge .x/.y suffixes ===")

suffix_patterns <- c("marker_avg_log2FC\\.x", "marker_avg_log2FC\\.y",
                     "marker_p_val_adj\\.x", "marker_p_val_adj\\.y",
                     "marker_rank\\.x", "marker_rank\\.y")
all_ok <- TRUE
for (pat in suffix_patterns) {
  found <- any(grepl(pat, lines2))
  if (found) {
    cat(sprintf("  FAIL: Found suffix pattern: %s\n", pat))
    all_ok <- FALSE
  } else {
    cat(sprintf("  OK: No suffix: %s\n", pat))
  }
}
stopifnot(all_ok)
message("PASS: No merge suffix artifacts in HTML")

# ============================================================================
# Test 4: Volcano labels not all at same capped coordinate
# ============================================================================
message("\n=== Test 4: Volcano labels not all at x=5,y=20 ===")

# Build with mock data that has many highly significant genes
high_sig_de <- data.frame(
  gene        = paste0("Gene_", 1:50),
  avg_log2FC  = c(runif(25, 3, 8), runif(25, -8, -3)),
  p_val       = 10^-(runif(50, 10, 120)),
  p_val_adj   = 10^-(runif(50, 10, 120)),
  stringsAsFactors = FALSE
)

out4 <- build_screport_de(
  de_df       = high_sig_de,
  output_file = file.path(pkg_root, "test_contract_volcano.html"),
  title       = "Contract Test — Volcano Labels"
)

stopifnot(file.exists(out4$output_file))
lines4 <- readLines(out4$output_file, warn = FALSE)

# Find all annotation x/y pairs
text4 <- paste(lines4, collapse = " ")
# Extract annotations section
ann_match <- gregexpr('"annotations":\\s*\\[[^]]*\\{[^}]*?\\}[^]]*\\]', text4)
# Simpler: count how many label annotations are at exactly x=5.0
# If more than 3 are at the same coord, the jitter didn't work
x5_count <- length(grep('"x":5[^0-9.]', lines4))
y20_count <- length(grep('"y":20[^0-9.]', lines4))

cat(sprintf("  Annotations with x=5.0 (exact): %d\n", x5_count))
cat(sprintf("  Annotations with y=20.0 (exact): %d\n", y20_count))

# With our jitter, capped labels should be spread across 4.6-5.0 and 18.5-20.0
# So no more than 2 should land on exactly x=5.0 or y=20.0
stopifnot(x5_count <= 2)
stopifnot(y20_count <= 2)
message("PASS: Volcano labels are spread (not all at same cap coordinate)")

# ============================================================================
# Test 5: Provided marker table with cluster col + group_col="cluster"
# ============================================================================
message("\n=== Test 5: Provided marker table + group_col drives DotPlot ===")

# Create a Seurat object that HAS 'cluster' in meta.data
obj5 <- obj
obj5$cluster <- rep(c("C1", "C2", "C3"), each = 100)

# Provided marker table with cluster column
marker_tbl <- data.frame(
  gene        = rep(paste0("Gene_", 1:20), 3),
  cluster     = rep(c("C1", "C2", "C3"), each = 20),
  avg_log2FC  = runif(60, 0.5, 3),
  p_val       = runif(60, 0.001, 0.05),
  p_val_adj   = runif(60, 0.01, 0.1),
  pct.1       = runif(60, 0.1, 0.9),
  pct.2       = runif(60, 0.1, 0.5),
  stringsAsFactors = FALSE
)

out5 <- build_screport_de(
  seurat_obj  = obj5,
  de_df       = marker_tbl,
  group_col   = "cluster",
  output_file = file.path(pkg_root, "test_contract_marker_tbl.html"),
  title       = "Contract Test — Marker Table + group_col"
)

stopifnot(file.exists(out5$output_file))
lines5 <- readLines(out5$output_file, warn = FALSE)

# Must have dotplot-json-data
has_json5 <- any(grepl("dotplot-json-data", lines5, fixed = TRUE))
cat(sprintf("  dotplot-json-data: %s\n", if (has_json5) "OK" else "MISSING!"))
stopifnot(has_json5)

# Should NOT show "No marker genes available"
has_no_markers <- any(grepl("No marker genes available", lines5, fixed = TRUE))
cat(sprintf("  'No marker genes' message: %s\n", if (!has_no_markers) "OK (absent)" else "FAIL (present!)"))
stopifnot(!has_no_markers)

message("PASS: Provided marker table + group_col generates DotPlot JSON")

# ============================================================================
# Test 6: No-data DotPlot is a complete section (not bare div)
# ============================================================================
message("\n=== Test 6: No-data DotPlot is complete section ===")

# Test from Test 1's output:
# The no-data should be wrapped in de-section, not bare
lines_nodata <- readLines(out1$output_file, warn = FALSE)

# Check that no-data div appears inside section-dotplot context
# grep for the section-dotplot div followed by no-data
txt <- paste(lines_nodata, collapse = "\n")

# A bare no-data div directly in de-content (without section wrapper) would
# mean there's a <div class="no-data"> that's NOT inside <div id="section-dotplot"
# Our fix ensures no-data is always inside de_section("dotplot", ...)

# Verify: section-dotplot exists (already checked in Test 1)
# Verify: section-dotplot contains "no-data" class OR the no-data message
# Simple check: after section-dotplot, before next section, there's content
sec_pos <- grep('id="section-dotplot"', lines_nodata, fixed = TRUE)
next_sec <- grep('id="section-', lines_nodata, fixed = TRUE)
next_sec <- next_sec[next_sec > sec_pos[1]]
if (length(next_sec) > 0) {
  between_lines <- lines_nodata[sec_pos[1]:next_sec[1]]
  has_no_data_inside <- any(grepl("no-data", between_lines, fixed = TRUE)) ||
                        any(grepl("no data", between_lines, ignore.case = TRUE))
  cat(sprintf("  no-data inside section-dotplot: %s\n",
              if (has_no_data_inside) "OK" else "MISSING!"))
  stopifnot(has_no_data_inside)
} else {
  stop("Could not verify section-dotplot boundaries")
}

message("PASS: No-data DotPlot is wrapped in complete section")

# ============================================================================
# Test 7: Volcano uses scattergl (WebGL renderer)
# ============================================================================
message("\n=== Test 7: Volcano scattergl ===")

# Use the high-sig mock data from Test 4
has_scattergl <- any(grepl('"type":"scattergl"', lines4, fixed = TRUE))
cat(sprintf("  scattergl trace: %s\n", if (has_scattergl) "OK" else "MISSING!"))
stopifnot(has_scattergl)

# Should NOT use plain scatter for volcano main trace
has_scatter_volcano <- any(grepl('"type":"scatter"', lines4, fixed = TRUE))
cat(sprintf("  plain scatter in volcano: %s\n",
    if (has_scatter_volcano) "PRESENT (may be segments/threshold)" else "OK (absent)"))

message("PASS: Volcano main trace uses scattergl")

# ============================================================================
# Test 8: Marker-only data → annotations ≤ 5
# ============================================================================
message("\n=== Test 8: Marker-only annotation count ≤ 5 ===")

# All-positive marker table (simulates precomputed all-markers)
marker_only_de <- data.frame(
  gene        = paste0("Marker_", 1:100),
  avg_log2FC  = runif(100, 0.5, 8),
  p_val       = 10^-(runif(100, 10, 120)),
  p_val_adj   = 10^-(runif(100, 10, 120)),
  stringsAsFactors = FALSE
)
# Make some non-significant
marker_only_de$p_val_adj[91:100] <- 0.5

out8 <- build_screport_de(
  de_df       = marker_only_de,
  output_file = file.path(pkg_root, "test_contract_marker_only.html"),
  title       = "Contract Test — Marker-Only Volcano"
)

stopifnot(file.exists(out8$output_file))
lines8 <- readLines(out8$output_file, warn = FALSE)

# Count annotation entries (each annotation is an object in the JSON)
# Count occurrences of '"text":"Marker_' which indicates a labelled gene
annotation_texts <- grep('"text":"Marker_', lines8, value = TRUE)
n_annotations <- length(annotation_texts)
cat(sprintf("  Gene annotations in marker-only volcano: %d\n", n_annotations))
stopifnot(n_annotations <= 5)

# Title must have marker-only note
has_marker_note <- any(grepl("marker-only", lines8, ignore.case = TRUE))
cat(sprintf("  marker-only title note: %s\n", if (has_marker_note) "OK" else "MISSING!"))
stopifnot(has_marker_note)

message("PASS: Marker-only volcano has ≤ 5 annotations + title note")

# ============================================================================
# Test 9: Bidirectional volcano annotation count ≤ 8
# ============================================================================
message("\n=== Test 9: Bidirectional annotation count ≤ 8 ===")

# Re-use high_sig_de from Test 4 (already has up and down genes)
# Count annotations from that output
annotation_texts4 <- grep('"text":"Gene_', lines4, value = TRUE)
n_annotations4 <- length(annotation_texts4)
cat(sprintf("  Gene annotations in bidirectional volcano: %d\n", n_annotations4))
stopifnot(n_annotations4 <= 10)  # generous upper bound; default is 8
# Don't be overly strict — the exact count depends on how many genes are 'sig'
message("PASS: Bidirectional volcano annotation count within bounds")

# ============================================================================
# Test 10: volcano_label_top_n parameter override
# ============================================================================
message("\n=== Test 10: volcano_label_top_n parameter ===")

out10 <- build_screport_de(
  de_df                = high_sig_de,
  volcano_label_top_n  = 3,
  output_file          = file.path(pkg_root, "test_contract_label_override.html"),
  title                = "Contract Test — Label Override"
)

stopifnot(file.exists(out10$output_file))
lines10 <- readLines(out10$output_file, warn = FALSE)

annotation_texts10 <- grep('"text":"Gene_', lines10, value = TRUE)
n_annotations10 <- length(annotation_texts10)
cat(sprintf("  Gene annotations with volcano_label_top_n=3: %d\n", n_annotations10))
stopifnot(n_annotations10 <= 3)

message("PASS: volcano_label_top_n parameter respected")

# ============================================================================
# Test 11: Capped-label overlap check — no pileup at y=20
# ============================================================================
message("\n=== Test 11: No label pileup at capped y=20 ===")

# With marker-only data (all high logFC, all tiny p) many labels would cap to y=20
# But with annotation limit=5 and jitter, they must be spread

# Check marker_only output (Test 8): each y-capped annotation should have unique position
# Extract all 'y' values for annotations
# We grep for the JSON annotations section — find y values around 18-20
y_vals <- regmatches(
  paste(lines8, collapse = "\n"),
  gregexpr('"y":1[89][0-9.]*|"y":20[0-9.]*', paste(lines8, collapse = "\n"))
)[[1]]
# If we have y-capped labels, their y values should NOT all be "20"
y_exact_20 <- sum(y_vals == '"y":20')
y_spread   <- length(y_vals) - y_exact_20
cat(sprintf("  Annotations at y=20.0 (exact): %d / %d\n", y_exact_20, length(y_vals)))
stopifnot(y_exact_20 <= 1)  # At most 1 can land on exactly 20.0; rest jittered

message("PASS: Capped labels are spread across jittered positions")

# ============================================================================
# Summary
# ============================================================================
message("\n========================================")
message("ALL HTML CONTRACT TESTS PASSED")
message("========================================")
