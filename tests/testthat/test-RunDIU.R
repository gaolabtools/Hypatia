test_that("RunDIU returns appropriate results", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  res <- RunDIU(gbm, min.gene.cts = 0, min.gene.pct = 0, min.tx.cts = 1)

  expect_class(res, "list")
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")
  expect_true("prop.diff" %in% names(res$data))
  expect_false("delta" %in% names(res$data))
  expect_true("max.prop.diff" %in% names(res$stats))
  expect_false("max.delta" %in% names(res$stats))

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
