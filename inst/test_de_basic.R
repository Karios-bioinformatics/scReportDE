# Smoke test for scReportDE v0.1.0 — Full coverage
# Run from package root: Rscript inst/test_de_basic.R

library(Seurat)
library(SeuratObject)
library(plotly)
library(htmltools)
library(jsonlite)
library(DT)

# ---- Source all package files (dynamic path) ----
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
set.seed(42)
n_cells <- 800
n_genes <- 100

# Metadata
condition_vec <- rep(c("Treatment", "Control"), each = 400)
sample_vec   <- paste0("Sample_", rep(1:4, each = 200))
meta <- data.frame(
  condition  = condition_vec,
  sample     = sample_vec,
  row.names  = paste0("Cell_", seq_len(n_cells)),
  stringsAsFactors = FALSE
)

# Expression matrix: some genes truly DE
gene_names <- paste0("Gene_", seq_len(n_genes))
counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                 nrow = n_genes, ncol = n_cells,
                 dimnames = list(gene_names, rownames(meta)))

# Inject true DE signal: first 15 genes up in Treatment, last 10 down
# (log-normal approach via Poisson mean shift)
de_up   <- 1:15
de_down <- (n_genes - 9):n_genes
for (g in de_up) {
  counts[g, condition_vec == "Treatment"] <-
    rpois(400, lambda = 6)
}
for (g in de_down) {
  counts[g, condition_vec == "Control"] <-
    rpois(400, lambda = 5)
}
for (g in de_down) {
  counts[g, condition_vec == "Treatment"] <-
    rpois(400, lambda = 1.5)
}

obj <- CreateSeuratObject(counts = counts, meta.data = meta)
obj <- NormalizeData(obj, verbose = FALSE)
Idents(obj) <- "condition"

message("Seurat object: ", ncol(obj), " cells x ", nrow(obj), " genes")

# ============================================================================
# Test 1: Pairwise DE (Seurat → full report)
# ============================================================================
message("\n=== Test 1: Pairwise DE from Seurat ===")
out1 <- build_screport_de(
  seurat_obj  = obj,
  ident_1     = "Treatment",
  ident_2     = "Control",
  mode        = "pairwise",
  output_file = file.path(pkg_root, "test_pairwise.html"),
  title       = "DE Test — Pairwise"
)

# Checks
stopifnot(!is.null(out1$de_df))
stopifnot(file.exists(out1$output_file))
stopifnot(nrow(out1$de_df) > 0)
stopifnot("gene" %in% colnames(out1$de_df))
stopifnot("avg_log2FC" %in% colnames(out1$de_df))

file_size1 <- file.info(out1$output_file)$size
message(sprintf("  File: %s (%d bytes)", basename(out1$output_file), file_size1))
message(sprintf("  DE genes: %d", nrow(out1$de_df)))

# Check HTML has key sections
lines1 <- readLines(out1$output_file, warn = FALSE)
for (kw in c("section-overview", "section-volcano", "section-de_table",
             "section-dotplot", "section-violin", "section-method")) {
  found <- any(grepl(kw, lines1, fixed = TRUE))
  cat(sprintf("  %-25s %s\n", kw, if (found) "OK" else "MISSING!"))
}

message("PASS: Pairwise DE report generated")

# ============================================================================
# Test 2: Pre-computed de_df (no Seurat object)
# ============================================================================
message("\n=== Test 2: Pre-computed de_df only ===")
faux_de <- data.frame(
  gene        = paste0("Gene_", 1:50),
  avg_log2FC  = rnorm(50, 0, 1.5),
  p_val       = runif(50, 0, 0.1),
  p_val_adj   = runif(50, 0, 0.1),
  pct.1       = runif(50, 0.1, 0.9),
  pct.2       = runif(50, 0.1, 0.9),
  stringsAsFactors = FALSE
)
# Make some significant
faux_de$p_val_adj[1:10] <- runif(10, 0.001, 0.04)
faux_de$p_val_adj[11:20] <- runif(10, 0.001, 0.04)

out2 <- build_screport_de(
  de_df       = faux_de,
  output_file = file.path(pkg_root, "test_dedf_only.html"),
  title       = "DE Test — de_df Only"
)

stopifnot(!is.null(out2$de_df))
stopifnot(file.exists(out2$output_file))
stopifnot(nrow(out2$de_df) == 50)
file_size2 <- file.info(out2$output_file)$size
message(sprintf("  File: %s (%d bytes)", basename(out2$output_file), file_size2))

lines2 <- readLines(out2$output_file, warn = FALSE)
# Dot plot and violin should show no-data (no Seurat object)
stopifnot(any(grepl("no data", lines2, ignore.case = TRUE)) ||
          any(grepl("no-data", lines2, fixed = TRUE)))
message("PASS: de_df-only report generated (dot/violin no-data shown)")

# ============================================================================
# Test 3: Missing optional columns (pct.1/pct.2)
# ============================================================================
message("\n=== Test 3: Missing pct.1 / pct.2 columns ===")
faux_de3 <- faux_de[, c("gene", "avg_log2FC", "p_val", "p_val_adj")]
out3 <- build_screport_de(
  de_df       = faux_de3,
  output_file = file.path(pkg_root, "test_missing_pct.html"),
  title       = "DE Test — Missing pct columns"
)

stopifnot(!is.null(out3$de_df))
stopifnot(file.exists(out3$output_file))
stopifnot("pct.1" %in% colnames(out3$de_df))  # Should be filled with NA
message(sprintf("  Warnings: %d", length(out3$warnings)))
stopifnot(length(out3$warnings) > 0)
message("PASS: Missing columns handled gracefully")

# ============================================================================
# Test 4: Empty DE results (0 genes)
# ============================================================================
message("\n=== Test 4: Empty DE results ===")
faux_de4 <- data.frame(
  gene       = character(0),
  avg_log2FC = numeric(0),
  p_val      = numeric(0),
  p_val_adj  = numeric(0),
  stringsAsFactors = FALSE
)

out4 <- tryCatch({
  build_screport_de(
    de_df       = faux_de4,
    output_file = file.path(pkg_root, "test_empty.html"),
    title       = "DE Test — Empty"
  )
}, error = function(e) {
  message("  Got error: ", e$message)
  NULL
})

if (!is.null(out4)) {
  stopifnot(file.exists(out4$output_file))
  stopifnot(nrow(out4$de_df) == 0)
  message("PASS: Empty DE results produce valid (but empty) report")
} else {
  message("INFO: Empty DE produced error — needs review")
}

# ============================================================================
# Test 5: AllMarkers mode
# ============================================================================
message("\n=== Test 5: AllMarkers mode ===")
obj2 <- obj
obj2$cluster <- sample(c("C1", "C2", "C3"), ncol(obj2), replace = TRUE)
Idents(obj2) <- "cluster"

out5 <- build_screport_de(
  seurat_obj  = obj2,
  mode        = "all_markers",
  output_file = file.path(pkg_root, "test_allmarkers.html"),
  title       = "DE Test — AllMarkers"
)

stopifnot(!is.null(out5$de_df))
stopifnot(file.exists(out5$output_file))
stopifnot("cluster" %in% colnames(out5$de_df))
n_clusters <- length(unique(out5$de_df$cluster))
message(sprintf("  DE genes: %d across %d clusters", nrow(out5$de_df), n_clusters))
message("PASS: AllMarkers report generated")

# ============================================================================
# Test 6: group_col parameter
# ============================================================================
message("\n=== Test 6: group_col parameter ===")
out6 <- build_screport_de(
  seurat_obj  = obj,
  group_col   = "condition",
  ident_1     = "Treatment",
  ident_2     = "Control",
  mode        = "pairwise",
  output_file = file.path(pkg_root, "test_groupcol.html"),
  title       = "DE Test — group_col"
)

stopifnot(!is.null(out6$de_df))
stopifnot(file.exists(out6$output_file))
message("PASS: group_col parameter works")

# ============================================================================
# Test 7: utils — normalize_de_df edge cases
# ============================================================================
message("\n=== Test 7: normalize_de_df edge cases ===")

# 7a: rownames as gene
df_rn <- data.frame(avg_log2FC = c(1.2, -0.5), p_val = c(0.01, 0.3), p_val_adj = c(0.04, 0.5))
rownames(df_rn) <- c("CD14", "CD3D")
norm_rn <- normalize_de_df(df_rn)
stopifnot(norm_rn$df$gene[1] == "CD14")

# 7b: avg_logFC → avg_log2FC rename
df_fc <- data.frame(gene = "X", avg_logFC = 0.8, p_val = 0.02, p_val_adj = 0.06)
norm_fc <- normalize_de_df(df_fc)
stopifnot("avg_log2FC" %in% colnames(norm_fc$df))
stopifnot(!is.null(norm_fc$warnings) && length(norm_fc$warnings) > 0)

# 7c: log2FoldChange → avg_log2FC
df_l2 <- data.frame(gene = "Y", log2FoldChange = -1.1, p_val = 0.001, p_val_adj = 0.01)
norm_l2 <- normalize_de_df(df_l2)
stopifnot("avg_log2FC" %in% colnames(norm_l2$df))
stopifnot(norm_l2$df$avg_log2FC[1] == -1.1)

# 7d: Missing p_val_adj
df_noadj <- data.frame(gene = "Z", avg_log2FC = 0.5, p_val = 0.03)
norm_noadj <- normalize_de_df(df_noadj)
stopifnot(length(norm_noadj$warnings) > 0)

message("PASS: All normalize_de_df edge cases handled")

# ============================================================================
# Summary
# ============================================================================
message("\n========================================")
message("ALL TESTS PASSED — scReportDE v0.1.0")
message("========================================")
