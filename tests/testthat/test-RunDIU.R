test_that("RunDIU returns appropriate results", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  res <- RunDIU(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    genes = "ENSG00000135945",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 1,
    permutation = FALSE,
    bootstrap = FALSE
  )

  expect_class(res, "list")
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")
  expect_true("prop.diff" %in% names(res$data))
  expect_false("delta" %in% names(res$data))
  expect_true("max.prop.diff" %in% names(res$stats))
  expect_false("max.delta" %in% names(res$stats))
  expect_true("cramers.v" %in% names(res$stats))
  expect_false("effect.size" %in% names(res$stats))

})

test_that("RunDIU requires cells on both sides of a comparison", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  expect_error(
    RunDIU(gbm, group.1 = unique(colData(gbm)$cell_type), quiet = TRUE),
    "at least one cell in both groups"
  )
})

test_that("RunDIU reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  gbm <- SetTranscripts(gbm, id = "transcript_name")

  res <- RunDIU(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 1,
    genes = "ENSG00000135945",
    quiet = TRUE
  )

  expect_true(all(res$data$transcript %in% rowData(gbm)$transcript_name))
  expect_true(all(res$stats$transcript %in% rowData(gbm)$transcript_name))
})

test_that("RunDIU supports stats::p.adjust methods", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  res <- RunDIU(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 1,
    genes = "ENSG00000135945",
    p.adj = "none",
    quiet = TRUE
  )

  expect_equal(res$stats$padj, res$stats$pval)
  expect_error(
    RunDIU(gbm, p.adj = "Bonferroni", quiet = TRUE),
    "element of set"
  )
})

test_that("RunDIU uses the uncorrected Pearson statistic for two-transcript genes", {
  countData <- matrix(
    c(
      5, 5, 10, 10,
      10, 10, 5, 5
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4))
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(countData)),
    data.frame(gene_id = rep("gene1", 2), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  result <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = FALSE,
    bootstrap = FALSE,
    p.adj = "none",
    quiet = TRUE
  )
  count_table <- matrix(c(10, 20, 20, 10), ncol = 2)
  pearson <- suppressWarnings(chisq.test(count_table, correct = FALSE))
  yates <- suppressWarnings(chisq.test(count_table, correct = TRUE))

  expect_equal(result$stats$pval, unname(pearson$p.value))
  expect_equal(
    result$stats$cramers.v,
    sqrt(unname(pearson$statistic) / sum(count_table))
  )
  expect_false(isTRUE(all.equal(result$stats$pval, unname(yates$p.value))))
})

test_that("RunDIU reports reproducible cell-label permutation p-values", {
  countData <- Matrix::Matrix(
    matrix(
      c(
        9, 3, 6, 8, 2, 4,
        1, 7, 4, 2, 8, 6,
        2, 1, 3, 5, 1, 2
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = rep("gene1", 3), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(123)
  first <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = TRUE,
    perm.iter = 20,
    bootstrap = FALSE,
    p.adj = "none",
    quiet = TRUE
  )
  set.seed(123)
  second <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = TRUE,
    perm.iter = 20,
    bootstrap = FALSE,
    p.adj = "none",
    quiet = TRUE
  )

  permutation_columns <- c(
    "pval.perm", "padj.perm", "pval.perm.delta", "padj.perm.delta"
  )
  expect_true(all(permutation_columns %in% names(first$stats)))
  expect_equal(first$stats[permutation_columns], second$stats[permutation_columns])
  expect_true(all(first$stats$pval.perm > 0 & first$stats$pval.perm <= 1))
  expect_true(all(first$stats$pval.perm.delta > 0 & first$stats$pval.perm.delta <= 1))
  expect_equal(first$stats$padj.perm, first$stats$pval.perm)
  expect_equal(first$stats$padj.perm.delta, first$stats$pval.perm.delta)

  observed_grp1 <- rowSums(countData[, 1:3, drop = FALSE])
  observed_grp2 <- rowSums(countData[, 4:6, drop = FALSE])
  observed_max_delta <- max(abs(
    observed_grp1 / sum(observed_grp1) - observed_grp2 / sum(observed_grp2)
  ))
  total_counts <- rowSums(countData)
  set.seed(123)
  null_max_delta <- replicate(20, {
    grp1_idx <- sample.int(ncol(countData))[1:3]
    perm_grp1 <- rowSums(countData[, grp1_idx, drop = FALSE])
    perm_grp2 <- total_counts - perm_grp1
    max(abs(
      perm_grp1 / sum(perm_grp1) - perm_grp2 / sum(perm_grp2)
    ))
  })
  expected_delta_pval <-
    (1 + sum(null_max_delta >= observed_max_delta)) / (1 + length(null_max_delta))
  expect_equal(first$stats$pval.perm.delta, expected_delta_pval)

  disabled <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = FALSE,
    bootstrap = FALSE,
    quiet = TRUE
  )
  expect_type(disabled$stats$pval.perm, "double")
  expect_true(all(is.na(disabled$stats$pval.perm)))
  expect_true(all(is.na(disabled$stats$padj.perm)))
  expect_type(disabled$stats$pval.perm.delta, "double")
  expect_true(all(is.na(disabled$stats$pval.perm.delta)))
  expect_true(all(is.na(disabled$stats$padj.perm.delta)))

  set.seed(123)
  delta_only <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = TRUE,
    perm.iter = 20,
    perm.statistics = "max_delta",
    bootstrap = FALSE,
    quiet = TRUE
  )
  expect_true(all(is.na(delta_only$stats$pval.perm)))
  expect_equal(delta_only$stats$pval.perm.delta, first$stats$pval.perm.delta)

  set.seed(123)
  chisq_only <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = TRUE,
    perm.iter = 20,
    perm.statistics = "chisq",
    bootstrap = FALSE,
    quiet = TRUE
  )
  expect_equal(chisq_only$stats$pval.perm, first$stats$pval.perm)
  expect_true(all(is.na(chisq_only$stats$pval.perm.delta)))
})

test_that("RunDIU supports permutation p-values alongside Fisher tests", {
  countData <- matrix(
    c(9, 3, 6, 8, 2, 4, 1, 7, 4, 2, 8, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("tx1", "tx2"), paste0("cell", 1:6))
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = c("gene1", "gene1"), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(321)
  result <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    method.use = "Fisher",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    permutation = TRUE,
    perm.iter = 10,
    bootstrap = FALSE,
    p.adj = "none",
    quiet = TRUE
  )

  expect_true(is.finite(result$stats$pval))
  expect_true(is.finite(result$stats$pval.perm))
  expect_true(is.finite(result$stats$pval.perm.delta))
  expect_equal(result$stats$padj, result$stats$pval)
  expect_equal(result$stats$padj.perm, result$stats$pval.perm)
  expect_equal(result$stats$padj.perm.delta, result$stats$pval.perm.delta)
})

test_that("RunDIU validates permutation iterations", {
  object <- CreateSCE(
    matrix(c(5, 4, 3, 2), nrow = 2,
           dimnames = list(c("tx1", "tx2"), c("cell1", "cell2"))),
    data.frame(group = c("A", "B"), row.names = c("cell1", "cell2")),
    data.frame(gene_id = c("gene1", "gene1"), row.names = c("tx1", "tx2")),
    active.group.id = "group",
    quiet = TRUE
  )

  expect_error(
    RunDIU(object, permutation = TRUE, perm.iter = 0, quiet = TRUE),
    ">= 1"
  )
  expect_error(
    RunDIU(object, perm.statistics = "not_a_statistic", quiet = TRUE),
    "subset"
  )
  expect_error(
    RunDIU(object, perm.statistics = character(), quiet = TRUE),
    "length"
  )
})

test_that("RunDIU returns typed empty bootstrap results when disabled", {

  countData <- Matrix::Matrix(
    matrix(
      c(
        9, 3, 6, 8, 2, 4,
        1, 7, 4, 2, 8, 6,
        2, 1, 3, 5, 1, 2
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  colData <- data.frame(
    group = rep(c("A", "B"), each = 3),
    row.names = colnames(countData)
  )
  rowData <- data.frame(
    gene_id = rep("gene1", 3),
    row.names = rownames(countData)
  )
  object <- CreateSCE(
    countData,
    colData,
    rowData,
    active.group.id = "group",
    quiet = TRUE
  )

  res <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    bootstrap = FALSE,
    cell.dispersion = FALSE,
    p.adj = "none",
    quiet = TRUE
  )

  expect_true(all(lengths(res$data$boot.prop.1) == 0))
  expect_true(all(lengths(res$data$boot.prop.2) == 0))
  expect_true(all(lengths(res$data$boot.prop.diff) == 0))
  expect_type(res$data$cell.n.1, "integer")
  expect_type(res$data$cell.n.2, "integer")
  expect_true(all(is.na(res$data$cell.n.1)))
  expect_true(all(is.na(res$data$cell.prop.sd.2)))
  expect_true(all(is.na(res$stats$boot.prop.diff.mean)))
  expect_true(all(is.na(res$stats$boot.prop.diff.lower)))
  expect_true(all(is.na(res$stats$boot.prop.diff.upper)))
  expect_true(all(res$stats$boot.valid.iter == 0L))
})

test_that("RunDIU reports cell-to-cell transcript proportion dispersion", {
  countData <- Matrix::Matrix(
    matrix(
      c(
        8, 0, 4, 0, 3, 1,
        2, 0, 0, 0, 1, 3,
        0, 0, 1, 0, 0, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = rep("gene1", 3), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  without_dispersion <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 2,
    bootstrap = FALSE,
    cell.dispersion = FALSE,
    p.adj = "none",
    quiet = TRUE
  )
  res_messages <- capture.output(
    res <- RunDIU(
      object,
      group.1 = "A",
      group.2 = "B",
      min.gene.cts = 0,
      min.gene.pct = 0,
      min.tx.cts = 2,
      bootstrap = FALSE,
      cell.dispersion = TRUE,
      p.adj = "none",
      quiet = FALSE
    ),
    type = "message"
  )
  expect_true(any(grepl("Calculating cell-to-cell dispersion...", res_messages, fixed = TRUE)))

  expect_equal(res$data$transcript, c("tx1", "tx2"))
  expect_equal(res$data$cell.n.1, c(2L, 2L))
  expect_equal(res$data$cell.n.2, c(2L, 2L))
  expect_equal(res$data$cell.prop.mean.1, c(0.9, 0.1), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.median.1, c(0.9, 0.1), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.sd.1, rep(stats::sd(c(0.8, 1)), 2), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.iqr.1, c(0.1, 0.1), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.mean.2, c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.zero.frac.1, c(0, 0.5), tolerance = 1e-12)
  expect_equal(res$data$cell.prop.zero.frac.2, c(0, 0), tolerance = 1e-12)
  legacy_cols <- setdiff(names(without_dispersion$data), grep("^cell\\.", names(without_dispersion$data), value = TRUE))
  expect_equal(res$data[legacy_cols], without_dispersion$data[legacy_cols], tolerance = 1e-12)
  expect_equal(res$stats, without_dispersion$stats, tolerance = 1e-12)

  empty <- Hypatia:::.CellUsageDispersion(
    Matrix::Matrix(matrix(0, nrow = 2, ncol = 2,
                          dimnames = list(c("tx1", "tx2"), c("c1", "c2"))), sparse = TRUE),
    c("gene1", "gene1")
  )
  expect_equal(empty$cell.n, c(0L, 0L))
  expect_true(all(is.na(empty$cell.prop.mean)))

  single <- Hypatia:::.CellUsageDispersion(
    Matrix::Matrix(matrix(c(1, 0), nrow = 2,
                          dimnames = list(c("tx1", "tx2"), "c1")), sparse = TRUE),
    c("gene1", "gene1")
  )
  expect_equal(single$cell.n, c(1L, 1L))
  expect_true(all(is.na(single$cell.prop.sd)))
})

test_that("RunDIU bootstraps transcript proportions without changing test statistics", {

  countData <- Matrix::Matrix(
    matrix(
      c(
        9, 3, 6, 8, 2, 4,
        1, 7, 4, 2, 8, 6,
        2, 1, 3, 5, 1, 2
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  colData <- data.frame(
    group = rep(c("A", "B"), each = 3),
    row.names = colnames(countData)
  )
  rowData <- data.frame(
    gene_id = rep("gene1", 3),
    row.names = rownames(countData)
  )
  object <- CreateSCE(
    countData,
    colData,
    rowData,
    active.group.id = "group",
    quiet = TRUE
  )

  without_bootstrap <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    p.adj = "none",
    quiet = TRUE
  )

  set.seed(1024)
  expected_prop_1 <- matrix(numeric(3 * 4), nrow = 3)
  expected_prop_2 <- matrix(numeric(3 * 4), nrow = 3)
  for (i in seq_len(4)) {
    grp1_idx <- sample(seq_len(3), size = 1, replace = TRUE)
    grp2_idx <- sample(seq_len(3), size = 1, replace = TRUE)
    expected_prop_1[, i] <- countData[, grp1_idx] / sum(countData[, grp1_idx])
    expected_prop_2[, i] <- countData[, grp2_idx + 3] / sum(countData[, grp2_idx + 3])
  }
  expected_prop_diff <- expected_prop_1 - expected_prop_2

  set.seed(1024)
  with_bootstrap <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    bootstrap = TRUE,
    boot.iter = 4,
    boot.ncells = 1,
    boot.conf = 0.8,
    p.adj = "none",
    quiet = TRUE
  )

  legacy_data_cols <- setdiff(names(without_bootstrap$data), grep("^boot\\.", names(without_bootstrap$data), value = TRUE))
  legacy_stats_cols <- c(
    "group.1", "group.2", "gene", "max.prop.diff", "transcript",
    "pval", "padj", "cramers.v", "approx"
  )
  expect_equal(
    with_bootstrap$data[legacy_data_cols],
    without_bootstrap$data[legacy_data_cols],
    tolerance = 1e-12
  )
  expect_equal(
    with_bootstrap$stats[legacy_stats_cols],
    without_bootstrap$stats[legacy_stats_cols],
    tolerance = 1e-12
  )

  expect_equal(with_bootstrap$data$boot.prop.1, lapply(seq_len(3), function(i) expected_prop_1[i, ]), tolerance = 1e-12)
  expect_equal(with_bootstrap$data$boot.prop.2, lapply(seq_len(3), function(i) expected_prop_2[i, ]), tolerance = 1e-12)
  expect_equal(with_bootstrap$data$boot.prop.diff, lapply(seq_len(3), function(i) expected_prop_diff[i, ]), tolerance = 1e-12)

  selected_idx <- match(with_bootstrap$stats$transcript, rownames(countData))
  selected_diff <- expected_prop_diff[selected_idx, ]
  expect_equal(with_bootstrap$stats$boot.prop.diff.mean, mean(selected_diff), tolerance = 1e-12)
  expect_equal(with_bootstrap$stats$boot.prop.diff.lower, unname(quantile(selected_diff, 0.1)), tolerance = 1e-12)
  expect_equal(with_bootstrap$stats$boot.prop.diff.upper, unname(quantile(selected_diff, 0.9)), tolerance = 1e-12)
  expect_equal(with_bootstrap$stats$boot.valid.iter, 4L)

  set.seed(2048)
  with_fraction <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    bootstrap = TRUE,
    boot.iter = 2,
    boot.fraction = 1,
    p.adj = "none",
    quiet = TRUE
  )
  set.seed(2048)
  with_fixed_count <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    bootstrap = TRUE,
    boot.iter = 2,
    boot.ncells = 3,
    p.adj = "none",
    quiet = TRUE
  )
  expect_equal(with_fraction$data$boot.prop.1, with_fixed_count$data$boot.prop.1)
  expect_equal(with_fraction$data$boot.prop.2, with_fixed_count$data$boot.prop.2)
  expect_equal(with_fraction$data$boot.prop.diff, with_fixed_count$data$boot.prop.diff)

  set.seed(4096)
  simulated_without_bootstrap <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    simulate.p = TRUE,
    p.adj = "none",
    quiet = TRUE
  )
  set.seed(4096)
  simulated_with_bootstrap <- RunDIU(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    simulate.p = TRUE,
    bootstrap = TRUE,
    boot.iter = 2,
    boot.ncells = 1,
    p.adj = "none",
    quiet = TRUE
  )
  expect_equal(
    simulated_with_bootstrap$stats[c("pval", "padj", "cramers.v")],
    simulated_without_bootstrap$stats[c("pval", "padj", "cramers.v")],
    tolerance = 1e-12
  )

  expect_error(
    RunDIU(object, bootstrap = TRUE, boot.fraction = 0.5, boot.ncells = 1, quiet = TRUE),
    "only one"
  )
  expect_error(
    RunDIU(object, bootstrap = TRUE, boot.conf = 1, quiet = TRUE),
    "strictly between"
  )
})
