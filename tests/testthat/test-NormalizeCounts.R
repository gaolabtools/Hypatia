test_that("NormalizeCounts returns object with normalized assay", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
      quiet = TRUE
    )

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    gbm <- NormalizeCounts(gbm, method.use = "TPM")
    gbm <- NormalizeCounts(gbm, method.use = "FT")
  })

  expect_true(all(assayNames(gbm) %in% c("counts", "logcounts", "tpmcounts", "ftcounts")))

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- NormalizeCounts(gbm, method.use = "TPM", gtf = gbm_gtf, gtf.transcript.id = "transcript_id"))

  expect_true(all(assayNames(gbm) %in% c("counts", "tpmcounts")))

})
