make_div_permutation_sce <- function(sparse = TRUE, active.transcripts = FALSE) {
  counts <- rbind(
    c(5, 5, 5, 10, 10, 10), c(5, 5, 5, 0, 0, 0),
    rep(5, 6), rep(5, 6),
    c(10, 10, 10, 5, 5, 5), c(0, 0, 0, 5, 5, 5)
  )
  dimnames(counts) <- list(paste0("tx", 1:6), paste0("cell", 1:6))
  if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
  CreateSCE(
    counts,
    data.frame(group = rep(c("Group A", "Group B"), each = 3),
               batch = "Batch 1", row.names = colnames(counts)),
    data.frame(gene_alias = rep(c("up", "equal", "down"), each = 2),
               tx_alias = paste0("isoform", 1:6), row.names = rownames(counts)),
    active.gene.id = "gene_alias",
    active.transcript.id = if (active.transcripts) "tx_alias" else NULL,
    active.group.id = "group", quiet = TRUE
  )
}

test_that("DIV permutations match the analytically known null and include ties", {
  # Three balanced cells versus three monoform cells. Tsallis(q=3) diversity
  # is 3/8 versus 0. Only 2 of the choose(6,3)=20 partitions attain |3/8|:
  # the original grouping and its reversal. The equal gene always ties zero.
  set.seed(321)
  k <- replicate(99, sum(sample.int(6)[1:3] <= 3))
  expected <- (1 + sum(k %in% c(0, 3))) / 100
  for (sparse in c(FALSE, TRUE)) {
    object <- make_div_permutation_sce(sparse)
    set.seed(321)
    result <- Hypatia:::.PermutationDiversityStats(
      assay(object)[, 1:3, drop = FALSE], assay(object)[, 4:6, drop = FALSE],
      rowData(object)$gene_alias, Hypatia:::.DiversityFunction("Tsallis"), 99
    )
    expect_equal(result$div.observed.1, c(3/8, 3/8, 0))
    expect_equal(result$div.observed.2, c(0, 3/8, 3/8))
    expect_equal(result$div.diff.observed, c(3/8, 0, -3/8))
    expect_equal(result$pval.perm, c(expected, 1, expected))
    expect_identical(result$perm.valid.iter, rep(99L, 3))
  }
})

test_that("optional DIV permutations preserve bootstrap results and identifiers", {
  for (sparse in c(FALSE, TRUE)) {
    for (active.transcripts in c(FALSE, TRUE)) {
      object <- make_div_permutation_sce(sparse, active.transcripts)
      set.seed(71)
      disabled <- RunDIV(object, boot.iter = 5, quiet = TRUE)
      expect_true(all(is.na(disabled$stats$pval.perm)))
      expect_true(all(is.na(disabled$stats$padj.perm)))
      expect_identical(disabled$stats$perm.valid.iter, rep(NA_integer_, 3))
      set.seed(71)
      enabled <- RunDIV(object, boot.iter = 5, permutation = TRUE,
                        perm.iter = 99, quiet = TRUE)
      added <- c("pval.perm", "padj.perm", "perm.valid.iter")
      expect_identical(enabled$data, disabled$data)
      expect_identical(enabled$stats[setdiff(names(enabled$stats), added)],
                       disabled$stats[setdiff(names(disabled$stats), added)])
      expect_equal(enabled$stats$padj.perm,
                   p.adjust(enabled$stats$pval.perm, "BH"))
      expect_true(all(enabled$stats$pval.perm >= 1/100))
      expect_identical(enabled$stats$perm.valid.iter, rep(99L, 3))
      set.seed(71)
      expect_identical(enabled, RunDIV(object, boot.iter = 5, permutation = TRUE,
                                      perm.iter = 99, quiet = TRUE))
    }
  }
})

test_that("permutations use original unequal group sizes, not bootstrap sizes", {
  counts <- rbind(c(5, 5, 10, 10, 10), c(5, 5, 0, 0, 0))
  # Of choose(5,2)=10 partitions, only the original two balanced cells
  # versus three monoform cells attains the observed absolute difference.
  set.seed(98)
  k <- replicate(50, sum(sample.int(5)[1:2] <= 2))
  set.seed(98)
  result <- Hypatia:::.PermutationDiversityStats(
    counts[, 1:2, drop = FALSE], counts[, 3:5, drop = FALSE],
    c("gene", "gene"), Hypatia:::.DiversityFunction("Tsallis"), 50
  )
  expect_equal(result$div.diff.observed, 3/8)
  expect_equal(result$pval.perm, (1 + sum(k == 2)) / 51)

  object <- make_div_permutation_sce()
  set.seed(23)
  a <- RunDIV(object, boot.iter = 3, permutation = TRUE, perm.iter = 30,
              boot.ncells = 1, div.diff.thresh = 0, quiet = TRUE)
  set.seed(23)
  b <- RunDIV(object, boot.iter = 3, permutation = TRUE, perm.iter = 30,
              boot.ncells = 1, div.diff.thresh = 1, quiet = TRUE)
  expect_identical(a$stats$pval.perm, b$stats$pval.perm)
  expect_equal(a$stats$div.diff.observed, c(-3/8, 0, 3/8))
})

test_that("DIV permutation adjustment is scoped to each comparison", {
  object <- make_div_permutation_sce()
  colData(object)$group <- rep(c("Group A", "Group B", "Group C"), each = 2)
  for (groups in list(list(NULL, NULL), list("Group A", NULL),
                      list(c("Group A", "Group B"), "Group C"),
                      list("Group A", "Group C"))) {
    set.seed(16)
    result <- RunDIV(object, group.by = c("group", "batch"),
      group.1 = if (is.null(groups[[1]])) NULL else paste0(groups[[1]], "_Batch 1"),
      group.2 = if (is.null(groups[[2]])) NULL else paste0(groups[[2]], "_Batch 1"),
      permutation = TRUE, perm.iter = 20, p.adj = "bonferroni",
      boot.iter = 2, quiet = TRUE)
    expected_rows <- if (is.null(groups[[1]])) 9L else 3L
    expect_equal(nrow(result$stats), expected_rows)
    for (comparison in split(result$stats, result$stats$group.1)) {
      expect_equal(comparison$padj.perm,
                   p.adjust(comparison$pval.perm, "bonferroni"))
    }
  }
})

test_that("undefined DIV permutation draws are counted without conditioning the null", {
  counts <- rbind(c(5, 0, 0, 10, 0, 0), c(5, 0, 0, 0, 0, 0))
  set.seed(123)
  draws <- replicate(50, seq_len(6) %in% sample.int(6)[1:3])
  set.seed(123)
  result <- Hypatia:::.PermutationDiversityStats(
    counts[, 1:3, drop = FALSE], counts[, 4:6, drop = FALSE],
    rep("gene", 2), Hypatia:::.DiversityFunction("Tsallis"), 50
  )
  expect_equal(result$div.diff.observed, 3/8)
  expect_identical(result$perm.valid.iter, as.integer(sum(draws[1, ] != draws[4, ])))
  expect_true(is.na(result$pval.perm))

  object <- make_div_permutation_sce()
  result <- RunDIV(object, entropy.use = "NormalizedShannon", permutation = TRUE,
                   perm.iter = 20, boot.iter = 2, quiet = TRUE)
  expect_true(all(is.na(result$stats$pval.perm[result$stats$gene != "equal"])))
  expect_true(all(is.na(result$stats$div.diff.observed[result$stats$gene != "equal"])))
  expect_equal(result$stats$pval.perm[result$stats$gene == "equal"], 1)
})

test_that("DIV observed and null diversities honor transcript filtering and entropy options", {
  object <- make_div_permutation_sce()
  # tx2 and tx6 have 15 total counts in one group: retained exactly at 15,
  # excluded at 16. A single retained transcript has zero Tsallis diversity.
  boundary <- RunDIV(object, min.tx.cts = 15, boot.iter = 2,
                     permutation = TRUE, perm.iter = 10, quiet = TRUE)
  removed <- RunDIV(object, min.tx.cts = 16, boot.iter = 2,
                    permutation = TRUE, perm.iter = 10, quiet = TRUE)
  expect_equal(boundary$stats$div.diff.observed, c(-3/8, 0, 3/8))
  expect_equal(removed$stats$div.diff.observed, c(0, 0))
  expect_equal(removed$data$n.transcripts, c(1L, 1L))
  expect_equal(removed$stats$pval.perm, c(1, 1))

  counts <- rbind(rep(6, 6), rep(3, 6), rep(1, 6))
  for (renormalize in c(FALSE, TRUE)) {
    result <- Hypatia:::.PermutationDiversityStats(
      counts[, 1:3], counts[, 4:6], rep("gene", 3),
      Hypatia:::.DiversityFunction("Shannon", top.n = 2, renormalize = renormalize), 5
    )
    proportions <- if (renormalize) c(2/3, 1/3) else c(0.6, 0.3)
    expect_equal(result$div.observed.1, -sum(proportions * log(proportions)))
    expect_equal(result$pval.perm, 1)
  }
})

test_that("DIV permutation inputs are validated", {
  object <- make_div_permutation_sce()
  expect_error(RunDIV(object, permutation = NA, quiet = TRUE))
  expect_error(RunDIV(object, permutation = 1, quiet = TRUE))
  for (value in c(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(RunDIV(object, perm.iter = value, quiet = TRUE))
  }
  expect_error(RunDIV(object, p.adj = "invalid", quiet = TRUE))
  expect_error(RunDIV(object, genes = "missing", permutation = TRUE, quiet = TRUE))
  expect_error(RunDIV(object, min.gene.cts = 1e6, permutation = TRUE, quiet = TRUE),
               "0 genes passed")
})
