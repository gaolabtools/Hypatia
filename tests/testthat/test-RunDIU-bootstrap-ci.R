make_diu_ci_sce <- function(counts, gene.ids, group, sparse = TRUE) {
  if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
  CreateSCE(
    counts,
    data.frame(group = group, row.names = colnames(counts)),
    data.frame(gene_alias = gene.ids, tx_alias = paste0("isoform", seq_len(nrow(counts))),
               row.names = rownames(counts)),
    active.gene.id = "gene_alias", active.transcript.id = "tx_alias",
    active.group.id = "group", quiet = TRUE
  )
}

test_that("RunDIU uses full-size resampling and signed percentile intervals", {
  counts <- rbind(
    c(9, 3, 5, 8, 2, 6, 4), c(2, 4, 6, 1, 3, 5, 7),
    c(1, 7, 5, 2, 8, 4, 6), c(8, 6, 4, 9, 7, 5, 3)
  )
  dimnames(counts) <- list(paste0("tx", 1:4), paste0("cell", 1:7))
  # Interleaved genes and unequal group sizes exercise pooling and row alignment.
  set.seed(302)
  expected <- replicate(40, {
    a <- rowSums(counts[, sample(1:3, size = 3, replace = TRUE), drop = FALSE])
    b <- rowSums(counts[, sample(1:4, size = 4, replace = TRUE) + 3, drop = FALSE])
    denominator.a <- c(a[1] + a[3], a[2] + a[4], a[1] + a[3], a[2] + a[4])
    denominator.b <- c(b[1] + b[3], b[2] + b[4], b[1] + b[3], b[2] + b[4])
    c(a / denominator.a, b / denominator.b)
  })
  for (sparse in c(FALSE, TRUE)) {
    object <- make_diu_ci_sce(counts, c("z", "a", "z", "a"),
                              c(rep("Group A", 3), rep("Group B", 4)), sparse)
    set.seed(302)
    result <- RunDIU(object, boot.iter = 40, quiet = TRUE)
    indices <- match(result$data$transcript, rowData(object)$tx_alias)
    expect_equal(result$data$boot.prop.1,
                 lapply(indices, function(i) expected[i, ]), tolerance = 1e-12)
    expect_equal(result$data$boot.prop.2,
                 lapply(indices, function(i) expected[i + 4, ]), tolerance = 1e-12)
    expect_equal(result$data$boot.prop.diff,
                 lapply(indices, function(i) expected[i, ] - expected[i + 4, ]),
                 tolerance = 1e-12)
    for (i in seq_len(nrow(result$stats))) {
      selected <- match(result$stats$transcript[i], rowData(object)$tx_alias)
      delta <- expected[selected, ] - expected[selected + 4, ]
      expect_equal(result$stats$boot.prop.diff.mean[i], mean(delta), tolerance = 1e-12)
      expect_equal(result$stats$boot.prop.diff.lower[i], unname(quantile(delta, 0.025)))
      expect_equal(result$stats$boot.prop.diff.upper[i], unname(quantile(delta, 0.975)))
    }
    expect_identical(result$stats$boot.valid.iter, rep(40L, 2))
    set.seed(302)
    expect_identical(result, RunDIU(object, boot.iter = 40, quiet = TRUE))

    set.seed(302)
    explicit <- RunDIU(object, boot.iter = 40, boot.fraction = 1, quiet = TRUE)
    expect_identical(result, explicit)
  }
})

test_that("RunDIU runs 5000 bootstrap iterations by default", {
  counts <- matrix(c(8, 2, 8, 2, 2, 8, 2, 8), nrow = 2,
                   dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4)))
  object <- make_diu_ci_sce(counts, rep("gene", 2), c("A", "A", "B", "B"))
  expect_equal(formals(RunDIU)$boot.iter, 5000)
  expect_equal(formals(RunDIU)$boot.fraction, 1)
  expect_equal(formals(RunDIU)$boot.conf, 0.95)
  set.seed(1024)
  result <- RunDIU(object, quiet = TRUE)
  expect_true(all(lengths(result$data$boot.prop.diff) == 5000L))
  expect_identical(result$stats$boot.valid.iter, 5000L)
  expect_equal(abs(result$stats$boot.prop.diff.mean), 0.6, tolerance = 1e-12)
  expect_equal(result$stats$boot.prop.diff.lower, result$stats$boot.prop.diff.mean,
               tolerance = 1e-12)
  expect_equal(result$stats$boot.prop.diff.upper, result$stats$boot.prop.diff.mean,
               tolerance = 1e-12)

  single <- RunDIU(object, boot.iter = 1, quiet = TRUE)
  expect_identical(single$stats$boot.valid.iter, 1L)
  expect_equal(abs(single$stats$boot.prop.diff.mean), 0.6, tolerance = 1e-12)
  expect_true(is.na(single$stats$boot.prop.diff.lower))
  expect_true(is.na(single$stats$boot.prop.diff.upper))
})

test_that("RunDIU requires two finite differences for interval bounds", {
  counts <- matrix(c(0, 0, 8, 2, 2, 8, 2, 8), nrow = 2,
                   dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4)))
  object <- make_diu_ci_sce(counts, rep("gene", 2), c("A", "A", "B", "B"))
  # Seed 1 samples the zero-count A cell in iterations 1, 2, and 4.
  set.seed(1)
  one <- RunDIU(object, min.gene.cts = 0, boot.iter = 4, boot.ncells = 1, quiet = TRUE)
  expect_identical(one$stats$boot.valid.iter, 1L)
  expect_equal(abs(one$stats$boot.prop.diff.mean), 0.6, tolerance = 1e-12)
  expect_true(is.na(one$stats$boot.prop.diff.lower))
  expect_true(is.na(one$stats$boot.prop.diff.upper))
  selected <- match(one$stats$transcript, one$data$transcript)
  expect_identical(which(is.finite(one$data$boot.prop.diff[[selected]])), 3L)

  set.seed(1)
  none <- RunDIU(object, min.gene.cts = 0, boot.iter = 2, boot.ncells = 1, quiet = TRUE)
  expect_identical(none$stats$boot.valid.iter, 0L)
  expect_true(is.na(none$stats$boot.prop.diff.mean))
  expect_true(is.na(none$stats$boot.prop.diff.lower))
  expect_true(is.na(none$stats$boot.prop.diff.upper))
})

test_that("DIU bootstrap settings do not change theoretical or permutation results", {
  counts <- matrix(c(8, 2, 6, 4, 2, 8, 3, 7), nrow = 2,
                   dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4)))
  object <- make_diu_ci_sce(counts, rep("gene", 2), c("A", "A", "B", "B"))
  set.seed(202)
  without <- RunDIU(object, bootstrap = FALSE, permutation = TRUE, perm.iter = 10,
                    quiet = TRUE)
  set.seed(202)
  with <- RunDIU(object, boot.iter = 10, permutation = TRUE, perm.iter = 10,
                 quiet = TRUE)
  columns <- c("max.prop.diff", "transcript", "pval", "padj", "cramers.v",
               "pval.perm", "padj.perm", "pval.perm.delta", "padj.perm.delta")
  expect_identical(with$stats[columns], without$stats[columns])
})
