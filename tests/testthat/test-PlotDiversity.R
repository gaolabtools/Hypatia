test_that("PlotDiversity works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    p_lp <- PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0)
    p_den <- PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id), plot.type = "density", min.tx.cts = 0)
    p_pc <- PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], plot.type = "pcoord", min.tx.cts = 0)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "Shannon")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "NormalizedShannon")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "Renyi")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "NormalizedRenyi")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "GiniSimpson")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, method.use = "InverseSimpson")
  })

  expect_class(p_lp, "ggplot")
  expect_class(p_den, "ggplot")
  expect_class(p_pc, "ggplot")

})
