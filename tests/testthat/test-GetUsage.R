test_that("GetUsage outputs a data frame", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error(res <- GetUsage(gbm, genes = c("ENSG00000135945", "ENSG00000048052", "ENSG00000049618")))

  expect_class(res, "data.frame")
  expect_true(nrow(res) >= 3)

})
