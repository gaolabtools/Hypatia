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

test_that("NormalizedShannon equals Shannon / log(2)", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  genes <- unique(gbm_rowData$gene_id)[1:5]

  res_sh <- GetDiversity(gbm, genes = genes, min.tx.cts = 0,
                         entropy.use = "Shannon")
  res_norm <- GetDiversity(gbm, genes = genes, min.tx.cts = 0,
                           entropy.use = "NormalizedShannon")

  expect_equal(res_norm$div, res_sh$div / log(2))
  expect_equal(res_norm$div.class, res_sh$div.class)

})

test_that("NormalizedRenyi equals Renyi / log(2)", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  genes <- unique(gbm_rowData$gene_id)[1:5]

  res_ren <- GetDiversity(gbm, genes = genes, min.tx.cts = 0,
                          entropy.use = "Renyi")
  res_norm <- GetDiversity(gbm, genes = genes, min.tx.cts = 0,
                           entropy.use = "NormalizedRenyi")

  expect_equal(res_norm$div, res_ren$div / log(2))
  expect_equal(res_norm$div.class, res_ren$div.class)

})
