input_compatibility_object <- function() {
  counts <- matrix(
    c(1, 2, 3, 4, 5, 6,
      6, 5, 4, 3, 2, 1,
      2, 3, 1, 5, 4, 6),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("tx3", "tx1", "tx2"), paste0("cell", 1:6))
  )
  CreateSCE(
    counts,
    data.frame(group = rep(c("A", "B"), each = 3),
               batch = rep(c("one", "two"), 3), row.names = colnames(counts)),
    data.frame(gene_id = c("g2", "g2", "g1"),
               transcript_id = c("iso3", "iso1", "iso2"),
               row.names = rownames(counts)),
    active.group.id = "group", quiet = TRUE
  )
}

test_that("normalization preserves zeros for dense and sparse count assays", {
  raw <- matrix(c(3, 1, 0, 0, 2, 2, 0, 0, 0), nrow = 3,
                dimnames = list(paste0("tx", 1:3), paste0("cell", 1:3)))
  object <- CreateSCE(
    raw, data.frame(group = c("A", "B", "B"), row.names = colnames(raw)),
    data.frame(gene_id = c("g1", "g1", "g2"), row.names = rownames(raw)),
    quiet = TRUE
  )
  stored_zeros <- Matrix::sparseMatrix(
    i = rep(1:3, 3), j = rep(1:3, each = 3), x = as.vector(raw),
    dims = dim(raw), dimnames = dimnames(raw)
  )
  expected_ft <- matrix(
    c(sqrt(3) + 2, 1 + sqrt(2), 0, 0, sqrt(2) + sqrt(3), sqrt(2) + sqrt(3), 0, 0, 0),
    nrow = 3, dimnames = dimnames(raw)
  )
  for (input in list(raw, Matrix::Matrix(raw, sparse = FALSE), counts(object), stored_zeros)) {
    counts(object) <- input
    for (method in c("LogNormalize", "FT")) {
      result <- NormalizeCounts(object, method.use = method, scale.factor = 4, quiet = TRUE)
      normalized <- assay(result, if (method == "FT") "ftcounts" else "logcounts")
      expect_equal(as.matrix(normalized), if (method == "FT") expected_ft else log1p(raw))
      expect_true(all(is.finite(normalized)))
      expect_identical(counts(result), input)
      expect_identical(dimnames(normalized), dimnames(raw))
      if (inherits(input, "sparseMatrix")) expect_s4_class(normalized, "sparseMatrix")
    }
  }
})

test_that("transcript subsetting recalculates QC for dense and sparse assays", {
  raw <- rbind(tx1 = c(2, 0, 0), tx2 = c(3, 4, 0), tx3 = c(1, 0, 0), tx4 = c(5, 5, 5))
  colnames(raw) <- paste0("cell", 1:3)
  object <- CreateSCE(
    raw, data.frame(group = c("A", "B", "B"), row.names = colnames(raw)),
    data.frame(gene_id = c("g1", "g1", "g2", "g3"), row.names = rownames(raw)),
    quiet = TRUE
  )
  stored_zeros <- Matrix::sparseMatrix(
    i = rep(1:4, 3), j = rep(1:3, each = 4), x = as.vector(raw),
    dims = dim(raw), dimnames = dimnames(raw)
  )
  for (input in list(raw, Matrix::Matrix(raw, sparse = FALSE), counts(object), stored_zeros)) {
    counts(object) <- input
    result <- SubsetTranscripts(object, transcripts = c("tx3", "tx1", "tx2"), quiet = TRUE)
    expect_identical(counts(result), input[1:3, , drop = FALSE])
    expect_equal(unname(result$nCount), c(6, 4, 0))
    expect_identical(result$nTranscript, c(3L, 1L, 0L))
    expect_identical(result$nGene, c(2L, 1L, 0L))
    empty <- SubsetTranscripts(object, transcripts = character(), quiet = TRUE)
    expect_equal(dim(empty), c(0L, 3L))
    expect_identical(empty$nGene, c(0L, 0L, 0L))
  }

  # Dense square assays may be stored as symmetric Matrix objects.
  symmetric_counts <- matrix(c(1, 2, 2, 0), 2,
                             dimnames = list(c("tx1", "tx2"), c("cell1", "cell2")))
  square_object <- object[1:2, 1:2, drop = FALSE]
  for (input in list(symmetric_counts, Matrix::Matrix(symmetric_counts, sparse = FALSE),
                     Matrix::Matrix(symmetric_counts, sparse = TRUE))) {
    counts(square_object) <- input
    result <- SubsetTranscripts(square_object, transcripts = c("tx1", "tx2"), quiet = TRUE)
    expect_equal(unname(result$nCount), c(3, 2))
    expect_identical(result$nTranscript, c(2L, 1L))
    expect_identical(result$nGene, c(1L, 1L))
  }
})

test_that("grouped functions reject missing labels before combining columns", {
  object <- input_compatibility_object()
  grouped_calls <- list(
    function(x, by) GetUsage(x, genes = "g2", group.by = by, quiet = TRUE),
    function(x, by) GetDiversity(x, genes = "g2", group.by = by, quiet = TRUE),
    function(x, by) GetExpression(x, genes = "g2", group.by = by,
                                 assay.use = "counts", quiet = TRUE),
    function(x, by) RunDIU(x, group.by = by, bootstrap = FALSE, quiet = TRUE),
    function(x, by) RunDIV(x, group.by = by, boot.iter = 2, quiet = TRUE),
    function(x, by) RunDEI(x, group.by = by, assay.use = "counts", quiet = TRUE),
    function(x, by) PlotUsage(x, gene = "g2", group.by = by, quiet = TRUE),
    function(x, by) PlotDiversity(x, genes = "g2", group.by = by, quiet = TRUE),
    function(x, by) PlotExpression(x, transcripts = "tx3", group.by = by,
                                  assay.use = "counts", quiet = TRUE),
    function(x, by) PlotCellQC(x, group.by = by, combine = FALSE)
  )
  for (as_factor in c(FALSE, TRUE)) {
    labels <- c("A", NA, "NA", "B", "B", "B")
    object$group <- if (as_factor) factor(labels) else labels
    for (run in grouped_calls) {
      expect_error(run(object, "group"), "[Mm]issing|anyMissing")
      expect_error(run(object, c("batch", "group")), "[Mm]issing|anyMissing")
    }
    expect_error(GetUsage(object, genes = "g2", quiet = TRUE), "[Mm]issing|anyMissing")
  }

  # A literal "NA" label remains a valid group when no value is missing.
  object$group <- c("A", "NA", "NA", "B", "B", "B")
  usage <- GetUsage(object, genes = "g2", group.by = "group", quiet = TRUE)
  expect_setequal(usage$group, c("A", "B", "NA"))
  expect_equal(usage$cts[usage$group == "NA" & usage$transcript == "tx3"], 5)
})

test_that("active gene column names with spaces preserve IDs and row alignment", {
  object <- input_compatibility_object()
  rowData(object)[["gene label"]] <- c("zeta", "zeta", "alpha")
  rowData(object)$gene_label <- rowData(object)[["gene label"]]
  reference <- SetGenes(object, "gene_label", quiet = TRUE)
  object <- SetGenes(object, "gene label", quiet = TRUE)

  for (id in c("", "transcript_id")) {
    actual <- SetTranscripts(object, id)
    expected <- SetTranscripts(reference, id)
    query <- if (id == "") c("tx2", "tx3", "tx1") else c("iso2", "iso3", "iso1")
    expect_equal(
      GetUsage(actual, genes = c("alpha", "zeta"), quiet = TRUE),
      GetUsage(expected, genes = c("alpha", "zeta"), quiet = TRUE)
    )
    expression <- GetExpression(actual, transcripts = query, assay.use = "counts", quiet = TRUE)
    expect_equal(expression,
                 GetExpression(expected, transcripts = query, assay.use = "counts", quiet = TRUE))
    expect_equal(expression$avgExpr[expression$group == "A" & expression$transcript == query[2]], 2)
    expect_true(all(expression$gene[expression$transcript == query[2]] == "zeta"))
    expect_equal(
      RunDEI(actual, transcripts = query, assay.use = "counts", quiet = TRUE),
      RunDEI(expected, transcripts = query, assay.use = "counts", quiet = TRUE)
    )
    plot <- PlotUsage(actual, gene = "zeta", quiet = TRUE)
    expect_equal(plot$data, PlotUsage(expected, gene = "zeta", quiet = TRUE)$data)
    expect_no_error(ggplot2::ggplot_build(plot))
  }
})

test_that("RunDEI tests a single requested transcript", {
  object <- input_compatibility_object()
  for (id in c("", "transcript_id")) {
    object <- SetTranscripts(object, id)
    transcript <- if (id == "") "tx3" else "iso3"
    result <- RunDEI(object, transcripts = transcript, assay.use = "counts", quiet = TRUE)
    expect_identical(result$stats$transcript, transcript)
    expect_equal(result$data$avgExpr.1, 2)
    expect_equal(result$data$avgExpr.2, 5)
    expect_equal(result$stats$log2FC, log2(2 / 5))
    expect_equal(result$stats$pval, 0.1)
    expect_equal(result$stats$padj, 0.1)
    expect_equal(result,
                 RunDEI(object, transcripts = c(transcript, "missing"),
                        assay.use = "counts", quiet = TRUE))
    expect_error(RunDEI(object, transcripts = "missing", assay.use = "counts", quiet = TRUE),
                 "None of the transcripts")
  }
})

test_that("PlotExpression supplies names for unnamed embedding dimensions", {
  object <- input_compatibility_object()
  coords <- matrix(c(6, 1, 5, 2, 4, 3, 9, 7, 8, 1, 3, 2), ncol = 2,
                   dimnames = list(colnames(object), NULL))
  for (dimension_names in list(NULL, c("UMAP1", "UMAP2"), c("", NA_character_))) {
    colnames(coords) <- dimension_names
    reducedDim(object, "UMAP") <- coords
    for (label in c(FALSE, TRUE)) {
      plot <- PlotExpression(object, transcripts = "tx3", assay.use = "counts",
                             plot.type = "reducedDim", dim.use = "UMAP",
                             group.subset = "B", label = label, quiet = TRUE)
      built <- ggplot2::ggplot_build(plot)
      expect_equal(built$data[[1]]$x, unname(coords[4:6, 1]))
      expect_equal(built$data[[1]]$y, unname(coords[4:6, 2]))
      if (label) {
        expect_equal(built$data[[2]]$x, mean(coords[4:6, 1]))
        expect_equal(built$data[[2]]$y, mean(coords[4:6, 2]))
      }
      expect_identical(reducedDim(object, "UMAP"), coords)
    }
  }
})
