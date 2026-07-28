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

})
