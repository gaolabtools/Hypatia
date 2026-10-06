make_diu_permutation_sce <- function(sparse = TRUE) {
  countData <- rbind(
    c(10, 0, 0, 10, 0, 0),
    c(10, 0, 0, 10, 0, 0),
    c(56, 32, 0, 8, 0, 0),
    c(8, 32, 0, 56, 0, 0),
    c(56, 56, 56, 8, 8, 8),
    c(8, 8, 8, 56, 56, 56)
  )
  dimnames(countData) <- list(paste0("tx", 1:6), paste0("cell", 1:6))
  if (sparse) countData <- Matrix::Matrix(countData, sparse = TRUE)
  CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = rep(c("flat", "shifted", "dense"), each = 2),
               row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )
}

test_that("RunDIU permutation p-values use finite draws separately for each gene", {
  set.seed(123)
  draws <- replicate(50, seq_len(6) %in% sample.int(6)[1:3])

  # The flat gene is expressed only in cells 1 and 4, with identical usage.
  # The shifted gene also has expression in cell 2, with intermediate usage.
  # Its observed statistics are 54 (Pearson) and 9/16 (maximum difference).
  # A valid shuffle ties these when cells 1 and 4 are in opposite groups;
  # isolating cell 2 instead gives zero for both statistics.
  flat_valid <- draws[1, ] != draws[4, ]
  shifted_valid <- colSums(draws[c(1, 2, 4), , drop = FALSE]) %in% 1:2
  # The dense gene has finite statistics in every draw. Only the original
  # partition and its reversal tie its observed maximum usage difference.
  dense_extreme <- colSums(draws[1:3, , drop = FALSE]) %in% c(0, 3)
  expect_equal(sum(flat_valid), 31)
  expect_equal(sum(shifted_valid), 45)
  expect_equal(sum(dense_extreme), 4)
  expected <- c(flat = 1, shifted = 32 / 46, dense = 5 / 51)

  for (sparse in c(FALSE, TRUE)) {
    object <- make_diu_permutation_sce(sparse)
    set.seed(123)
    result <- RunDIU(
      object,
      permutation = TRUE,
      perm.iter = 50,
      bootstrap = FALSE,
      only.valid = TRUE,
      quiet = TRUE
    )

    expect_setequal(result$stats$gene, names(expected))
    expect_true(all(result$stats$approx == "valid"))
    expected_p <- unname(expected[result$stats$gene])
    expect_equal(result$stats$pval.perm, expected_p)
    expect_equal(result$stats$pval.perm.delta, expected_p)
    expect_equal(result$stats$padj.perm, p.adjust(expected_p, "BH"))
    expect_equal(result$stats$padj.perm.delta, p.adjust(expected_p, "BH"))
  }
})

test_that("RunDIU returns NA permutation p-values when no draws are finite", {
  object <- make_diu_permutation_sce()
  set.seed(1)
  expect_identical(sample.int(6)[1:3], c(1L, 4L, 3L))
  # Both cells expressing the flat gene enter group A in the only shuffle.
  set.seed(1)
  result <- RunDIU(
    object,
    permutation = TRUE,
    perm.iter = 1,
    bootstrap = FALSE,
    quiet = TRUE
  )

  expected <- c(flat = NA_real_, shifted = 0.5, dense = 0.5)
  expect_setequal(result$stats$gene, names(expected))
  expected_p <- unname(expected[result$stats$gene])
  expect_equal(result$stats$pval.perm, expected_p)
  expect_equal(result$stats$pval.perm.delta, expected_p)
  expect_equal(result$stats$padj.perm, p.adjust(expected_p, "BH"))
  expect_equal(result$stats$padj.perm.delta, p.adjust(expected_p, "BH"))
  expect_true(all(is.finite(result$stats$pval)))
  expect_true(all(is.finite(result$stats$padj)))
})
