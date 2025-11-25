test_that("PlotUsage outputs a ggplot object", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(p <- PlotUsage(gbm, gene = "ENSG00000135945", group.by = "cell_type"))

  expect_class(p, "ggplot")

})
