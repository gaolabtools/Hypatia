test_that("RunDIV works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    res <- RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10)
    RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "Shannon", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, boot.iter = 10, entropy.use = "NormalizedShannon", entropy.thresh = 0.6, quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "Renyi", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "NormalizedRenyi", entropy.thresh = 0.6, quiet = TRUE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "GiniSimpson", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "InverseSimpson", quiet = FALSE)
  })

  expect_class(res, "list")
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")
  expect_true("div.diff" %in% names(res$stats))
  expect_false("delta.div" %in% names(res$stats))

})

test_that("RunDIV requires cells on both sides of a comparison", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  expect_error(
    RunDIV(gbm, group.1 = unique(colData(gbm)$cell_type), quiet = TRUE),
    "at least one cell in both groups"
  )
})

test_that("RunDIV works with active transcript IDs", {

  countData <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 2, 3, 4),
    x = c(5, 3, 2, 6),
    dims = c(2, 4),
    dimnames = list(paste0("tx", 1:2), paste0("cell", 1:4))
  )
  colData <- data.frame(group = c("A", "A", "B", "B"), row.names = colnames(countData))
  rowData <- data.frame(
    gene_id = rep("gene1", 2),
    tx_name = c("tx_a", "tx_b"),
    row.names = rownames(countData)
  )
  object <- CreateSCE(countData, colData, rowData, active.group.id = "group", quiet = TRUE)
  object <- SetTranscripts(object, id = "tx_name")

  expect_no_error({
    res <- RunDIV(
      object,
      group.1 = "A",
      group.2 = "B",
      min.gene.cts = 0,
      min.gene.pct = 0,
      min.tx.cts = 0,
      boot.iter = 3,
      boot.size = 1,
      genes = "gene1",
      p.adj = "none",
      quiet = TRUE
    )
  })
  expect_equal(res$stats$gene, "gene1")
  expect_equal(res$stats$padj, res$stats$pval)
})

test_that("RunDIV rejects p.adjust method names not used by stats::p.adjust", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  expect_error(
    RunDIV(gbm, p.adj = "Bonferroni", quiet = TRUE),
    "element of set"
  )
})
