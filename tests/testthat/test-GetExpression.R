test_that("GetExpression outputs a data frame", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    res <- GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10))
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = "Tumor")
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = c("Tumor", "Astrocyte"))
  })

  expect_class(res, "data.frame")

})
