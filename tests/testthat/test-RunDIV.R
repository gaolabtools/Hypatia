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
    res <- RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0)
    RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "Shannon", quiet = TRUE)
    #RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "NormalizedShannon", diversity.cutoff = 0.6, quiet = TRUE)
    RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "Renyi", quiet = TRUE)
    #RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "NormalizedRenyi", diversity.cutoff = 0.6, quiet = TRUE)
    RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "GiniSimpson", quiet = TRUE)
    RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, method.use = "InverseSimpson", quiet = TRUE)
  })

  expect_class(res, "list")
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")

})
