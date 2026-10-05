make_div_bootstrap_sce <- function(sparse = TRUE) {
  counts <- matrix(
    c(5, 5, 10, 0, 5, 5, 8, 2, 10, 0, 5, 5, 5, 5, 10, 0),
    nrow = 8, ncol = 2,
    dimnames = list(paste0("tx", 1:8), c("cell1", "cell2"))
  )
  counts <- counts[, rep(1:2, each = 3)]
  colnames(counts) <- paste0("cell", 1:6)
  if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
  CreateSCE(
    counts,
    data.frame(group = rep(c("Group A", "Group B"), each = 3),
               row.names = colnames(counts)),
    data.frame(gene_alias = rep(c("up", "down", "equal", "below"), each = 2),
               tx_alias = paste0("isoform", 1:8), row.names = rownames(counts)),
    active.gene.id = "gene_alias", active.transcript.id = "tx_alias",
    active.group.id = "group", quiet = TRUE
  )
}

test_that("bootstrap support uses finite differences from matching iterations", {
  summarize <- Hypatia:::.BootstrapDiversitySupport
  result <- summarize(
    c(0.4, NA, 0.2, 0.1, Inf), c(0.1, 0.2, NA, 0.4, 0), 0.2, 0.975
  )
  expect_identical(result$boot.valid.iter, 2L)
  expect_equal(result$div.diff, 0, tolerance = 1e-12)
  expect_equal(result$div.diff.lower, -0.285, tolerance = 1e-12)
  expect_equal(result$div.diff.upper, 0.285, tolerance = 1e-12)
  expect_equal(result$support.positive, 0.5)
  expect_equal(result$support.negative, 0.5)
  expect_false(result$supported)

  none <- summarize(c(NA, 0.2), c(0.3, NA), 0.1, 0.975)
  expect_identical(none$boot.valid.iter, 0L)
  expect_true(all(is.na(none[setdiff(names(none), "boot.valid.iter")])))
  expect_type(none$supported, "logical")
  expect_type(none$div.diff.lower, "double")
  one <- summarize(c(0.4, NA), c(0.1, 0.2), 0.2, 0.975)
  expect_identical(one$boot.valid.iter, 1L)
  expect_equal(one$support.positive, 1)
  expect_true(is.na(one$div.diff.lower))
  expect_true(is.na(one$supported))
})

test_that("effect exceedance is strict and the support cutoff is inclusive", {
  summarize <- Hypatia:::.BootstrapDiversitySupport
  boundary <- summarize(c(0.125, -0.125, 0), c(0, 0, 0), 0.125, 0.975)
  expect_equal(boundary$support.positive, 0)
  expect_equal(boundary$support.negative, 0)
  expect_false(boundary$supported)
  at_cutoff <- summarize(c(rep(0.25, 39), 0), rep(0, 40), 0.125, 0.975)
  expect_equal(at_cutoff$support.positive, 0.975)
  expect_true(at_cutoff$supported)
  below_cutoff <- summarize(c(rep(0.25, 38), 0, 0), rep(0, 40), 0.125, 0.975)
  expect_false(below_cutoff$supported)
  stricter <- summarize(c(rep(0.25, 39), 0), rep(0, 40), 0.125, 0.99)
  expect_false(stricter$supported)
})

test_that("RunDIV runs the full default bootstrap and honors combined grouping", {
  object <- make_div_bootstrap_sce()
  colData(object)$batch <- rep("Batch 1", ncol(object))
  set.seed(42)
  result <- RunDIV(object, group.by = c("group", "batch"), quiet = TRUE)
  expect_identical(result$stats$boot.valid.iter, rep(5000L, 4))
  expect_true(all(lengths(result$data$div.1) == 5000L))
  expect_true(all(result$stats$group.1 == "Group A_Batch 1"))
  expect_true(all(result$stats$group.2 == "Group B_Batch 1"))
  expect_identical(result$stats$supported, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(result$stats$div.diff, c(0.24, -0.375, 0, 0.375), tolerance = 1e-12)
})

test_that("RunDIV defaults to directional support on the retained gene diversities", {
  for (sparse in c(FALSE, TRUE)) {
    object <- make_div_bootstrap_sce(sparse)
    set.seed(1024)
    result <- RunDIV(object, boot.iter = 20, div.diff.thresh = 0.25, quiet = TRUE)
    stats <- result$stats[match(c("up", "down", "equal", "below"), result$stats$gene), ]
    expect_equal(stats$div.diff, c(0.375, -0.375, 0, 0.24), tolerance = 1e-12)
    expect_equal(stats$div.diff.lower, stats$div.diff, tolerance = 1e-12)
    expect_equal(stats$div.diff.upper, stats$div.diff, tolerance = 1e-12)
    expect_equal(stats$support.positive, c(1, 0, 0, 0))
    expect_equal(stats$support.negative, c(0, 1, 0, 0))
    expect_identical(stats$supported, c(TRUE, TRUE, FALSE, FALSE))
    expect_identical(stats$boot.valid.iter, rep(20L, 4))
    expect_false(any(c("pval", "padj") %in% names(stats)))
    expect_equal(result$data$n.transcripts, rep(2L, 4))
    expect_true(all(result$data$group.1 == "Group A"))
    expect_true(all(result$data$group.2 == "Group B"))
    expect_true(inherits(assay(object, "counts"), "sparseMatrix"))
    set.seed(1024)
    expect_identical(result, RunDIV(object, boot.iter = 20,
                                  div.diff.thresh = 0.25, quiet = TRUE))
  }
})

test_that("default bootstrap draws preserve each original group size", {
  counts <- matrix(c(9, 1, 3, 7, 5, 5, 8, 2, 2, 8, 6, 4, 4, 6), nrow = 2,
                   dimnames = list(c("tx1", "tx2"), paste0("cell", 1:7)))
  object <- CreateSCE(
    counts,
    data.frame(group = c(rep("A", 3), rep("B", 4)), row.names = colnames(counts)),
    data.frame(gene_id = c("gene", "gene"), row.names = rownames(counts)),
    active.group.id = "group", quiet = TRUE
  )
  entropy <- function(x) (1 - sum((x / sum(x))^3)) / 2
  set.seed(302)
  expected <- replicate(40, {
    a <- sample(1:3, size = 3, replace = TRUE)
    b <- sample(1:4, size = 4, replace = TRUE) + 3
    c(entropy(rowSums(counts[, a, drop = FALSE])),
      entropy(rowSums(counts[, b, drop = FALSE])))
  })
  set.seed(302)
  result <- RunDIV(object, boot.iter = 40, quiet = TRUE)
  expect_equal(result$data$div.1[[1]], expected[1, ], tolerance = 1e-12)
  expect_equal(result$data$div.2[[1]], expected[2, ], tolerance = 1e-12)
  delta <- expected[1, ] - expected[2, ]
  expect_equal(result$stats$div.diff, mean(delta), tolerance = 1e-12)
  expect_equal(result$stats$div.diff.lower, unname(quantile(delta, 0.025)))
  expect_equal(result$stats$div.diff.upper, unname(quantile(delta, 0.975)))
  expect_equal(result$stats$support.positive, mean(delta > 0.1))
  expect_equal(result$stats$support.negative, mean(delta < -0.1))
})

test_that("invalid bootstrap genes remain visible with unavailable support", {
  object <- make_div_bootstrap_sce()
  set.seed(202)
  result <- RunDIV(object, entropy.use = "NormalizedShannon", boot.iter = 4, quiet = TRUE)
  missing <- result$stats[result$stats$gene %in% c("up", "down", "below"), ]
  expect_equal(nrow(result$stats), 4L)
  expect_identical(missing$boot.valid.iter, rep(0L, 3))
  expect_true(all(is.na(missing$supported)))
  expect_true(all(is.na(missing$support.positive)))
  expect_true(all(is.na(missing$div.diff)))
  expect_true(all(is.na(missing$div.diff.lower)))
  expect_false(any(is.nan(result$stats$avgDiv.1)))
  expect_true(all(result$stats$boot.valid.iter[result$stats$gene == "equal"] == 4L))
  expect_error(RunDIV(object, min.gene.cts = 1e6, boot.iter = 4, quiet = TRUE),
               "0 genes passed")
})

test_that("bootstrap method arguments are validated", {
  object <- make_div_bootstrap_sce()
  expect_false(any(c("test.use", "p.adj") %in% names(formals(RunDIV))))
  expect_equal(formals(RunDIV)$boot.iter, 5000)
  expect_equal(formals(RunDIV)$boot.fraction, 1)
  for (threshold in list(-1, NA_real_, Inf)) {
    expect_error(RunDIV(object, div.diff.thresh = threshold, quiet = TRUE))
  }
  for (support in list(0.5, 0, 1.01, NA_real_, Inf)) {
    expect_error(RunDIV(object, support.thresh = support, quiet = TRUE))
  }
  expect_error(RunDIV(object, boot.iter = 1.5, quiet = TRUE), "count")
  expect_error(RunDIV(object, boot.fraction = 0.5, boot.ncells = 2, quiet = TRUE),
               "only one")
})
