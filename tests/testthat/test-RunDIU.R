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

})
