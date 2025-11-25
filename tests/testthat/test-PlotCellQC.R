test_that("PlotCellQC outputs QC metrics", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error({
    p <- PlotCellQC(gbm)
    p_list <- PlotCellQC(gbm, combine = FALSE)

  })

  expect_class(p, "patchwork")
  expect_class(p_list, "list")

})
