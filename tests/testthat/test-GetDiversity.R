test_that("GetDiversity works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    res <- GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Shannon", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "NormalizedShannon", entropy.thresh = 0.6, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Renyi", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "NormalizedRenyi", entropy.thresh = 0.6, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "GiniSimpson", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "InverseSimpson", quiet = TRUE)
  })

  expect_class(res, "data.frame")

})

test_that("GetDiversity reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  gbm <- SetTranscripts(gbm, id = "transcript_name")

  res <- GetDiversity(
    gbm,
    genes = "ENSG00000135945",
    group.subset = "Tumor",
    min.tx.cts = 0,
    quiet = TRUE
  )

  expect_true(all(res$transcript %in% rowData(gbm)$transcript_name))
})
