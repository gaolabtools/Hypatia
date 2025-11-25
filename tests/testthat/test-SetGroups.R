test_that("SetGroups updates the active.group.id object metadata", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- SetGroups(gbm, id = "cell_type"))

  expect_true(metadata(gbm)$active.group.id == "cell_type")

})
