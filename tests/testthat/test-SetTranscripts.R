test_that("SetTranscripts updates the active.transcript.id object metadata", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error({
    gbm_txname <- SetTranscripts(gbm, id = "transcript_name")
    gbm_reset <- SetTranscripts(gbm_txname, id = "")
  })

  expect_true({
    metadata(gbm)$active.transript.id == ""
    metadata(gbm_txname)$active.transcript.id == "transcript_name"
    metadata(gbm_reset)$active.transcript.id == ""
  })


})
