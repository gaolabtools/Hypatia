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

test_that("GetUsage reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  gbm <- SetTranscripts(gbm, id = "transcript_name")

  res <- GetUsage(
    gbm,
    genes = "ENSG00000135945",
    group.subset = "Tumor",
    quiet = TRUE
  )

  expect_true(all(res$transcript %in% rowData(gbm)$transcript_name))
})
