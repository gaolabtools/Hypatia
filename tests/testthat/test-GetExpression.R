test_that("GetExpression outputs a data frame", {

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
    res <- GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10))
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = "Tumor")
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = c("Tumor", "Astrocyte"))
  })

  expect_class(res, "data.frame")
  expect_true("pct" %in% names(res))
  expect_false("transcript.pct" %in% names(res))

})

test_that("GetExpression reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  gbm <- SetTranscripts(gbm, id = "transcript_name")

  transcripts <- rowData(gbm)$transcript_name[1:5]
  res <- GetExpression(gbm, transcripts = transcripts, group.subset = "Tumor", quiet = TRUE)

  expect_setequal(res$transcript, transcripts)
  expect_true(all(res$pct >= 0 & res$pct <= 1))
})

test_that("GetExpression works with gene filters and grouped names", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  colData(gbm)[["cell type"]] <- colData(gbm)$cell_type

  res <- GetExpression(
    gbm,
    genes = "ENSG00000135945",
    group.by = "cell type",
    group.subset = "Tumor",
    quiet = TRUE
  )

  expect_true(all(res$gene == "ENSG00000135945"))
  expect_true("pct" %in% names(res))
})
