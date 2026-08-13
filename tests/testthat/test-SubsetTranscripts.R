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

test_that("SubsetTranscripts requires exactly one selector", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_error(SubsetTranscripts(gbm, quiet = TRUE), "exactly one")
  expect_error(
    SubsetTranscripts(gbm, subset = nCell >= 3, transcripts = rownames(gbm)[1:5], quiet = TRUE),
    "exactly one"
  )
})

test_that("SubsetTranscripts supports non-syntactic rowData names", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)[["structural category"]] <- rowData(gbm)$structural_category

  expect_no_error({
    gbm_sq <- SubsetTranscripts(
      gbm,
      subset = `structural category` == "full-splice_match",
      quiet = TRUE
    )
  })
  expect_true(all(rowData(gbm_sq)[["structural category"]] == "full-splice_match"))
})

test_that("SubsetTranscripts preserves rowData columns when active transcript IDs are used", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)$rn_temp_holder <- "keep_me"
  gbm <- SetTranscripts(gbm, id = "transcript_name")
  selected_transcripts <- rowData(gbm)$transcript_name[1:5]

  gbm_subset <- SubsetTranscripts(
    gbm,
    transcripts = selected_transcripts,
    quiet = TRUE
  )

  expect_equal(rownames(gbm_subset), rownames(gbm)[1:5])
  expect_true("rn_temp_holder" %in% names(rowData(gbm_subset)))
  expect_true(all(rowData(gbm_subset)$rn_temp_holder == "keep_me"))
})

test_that("SubsetTranscripts recomputes cell QC metrics after transcript filtering", {

  countData <- Matrix::sparseMatrix(
    i = c(1, 2, 3),
    j = c(1, 3, 3),
    x = c(1, 2, 3),
    dims = c(3, 3),
    dimnames = list(paste0("tx", 1:3), paste0("cell", 1:3))
  )
  colData <- data.frame(group = c("A", "B", "A"), row.names = colnames(countData))
  rowData <- data.frame(gene_id = c("gene1", "gene1", "gene2"), row.names = rownames(countData))

  object <- CreateSCE(countData, colData, rowData, quiet = TRUE)
  object <- SubsetTranscripts(object, transcripts = c("tx2", "tx3"), quiet = TRUE)

  expect_equal(unname(colData(object)$nCount), c(0, 0, 5))
  expect_equal(colData(object)$nTranscript, c(0L, 0L, 2L))
  expect_equal(colData(object)$nGene, c(0L, 0L, 2L))
})

test_that("SubsetTranscripts reports invalid subset expressions", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_error(
    SubsetTranscripts(gbm, subset = missing_column == "x", quiet = TRUE),
    "subset term was not valid"
  )
})
