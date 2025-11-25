test_that("RunDEI outputs a data frame", {

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
    res <- RunDEI(gbm)
    RunDEI(gbm, group.1 = "Tumor")
    RunDEI(gbm, group.1 = "Tumor", group.2 = "Oligodendrocyte")
    RunDEI(gbm, group.1 = "Tumor", group.2 = c("Oligodendrocyte", "Astrocyte"))
    })


  expect_class(res, "data.frame")
  expect_true(nrow(res) >= 3)

})
