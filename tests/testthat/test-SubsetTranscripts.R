test_that("SubsetTranscripts works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error({
    gbm_nCell <- SubsetTranscripts(gbm, subset = nCell >= 3)
    gbm_nCell_inv <- SubsetTranscripts(gbm, subset = nCell >= 3, invert = TRUE)
    gbm_sq <- SubsetTranscripts(gbm, subset = structural_category == "full-splice_match")
    gbm_feat <- SubsetTranscripts(gbm, transcripts = rownames(gbm)[1:5])
  })

  expect_true({
    all(rowData(gbm_nCell)$nCell >= 3)
    all(rowData(gbm_nCell_inv)$nCell < 3)
    all(rowData(gbm_sq)$structural_category == "full-splice_match")
    all(rownames(gbm_feat) == rownames(gbm)[1:5])
  })

})
