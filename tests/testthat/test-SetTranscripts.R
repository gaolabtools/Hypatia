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

  expect_equal(metadata(gbm)$active.transcript.id, "")
  expect_equal(metadata(gbm_txname)$active.transcript.id, "transcript_name")
  expect_equal(metadata(gbm_reset)$active.transcript.id, "")
})

test_that("SetTranscripts validates transcript ID columns", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)$tx_missing <- rowData(gbm)$transcript_name
  rowData(gbm)$tx_missing[1] <- NA_character_
  rowData(gbm)$tx_duplicate <- rowData(gbm)$transcript_name
  rowData(gbm)$tx_duplicate[2] <- rowData(gbm)$tx_duplicate[1]
  rowData(gbm)$tx_factor <- factor(rowData(gbm)$transcript_name)

  expect_error(SetTranscripts(gbm, id = "not_a_column"))
  expect_error(SetTranscripts(gbm, id = "tx_missing"))
  expect_error(SetTranscripts(gbm, id = "tx_duplicate"))
  expect_error(SetTranscripts(gbm, id = "tx_factor"))

})
