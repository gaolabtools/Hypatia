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
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "Shannon", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "NormalizedShannon", diversity.cutoff = 0.6, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "Renyi", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "NormalizedRenyi", diversity.cutoff = 0.6, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "GiniSimpson", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, method.use = "InverseSimpson", quiet = TRUE)
  })

  expect_class(res, "data.frame")

})
