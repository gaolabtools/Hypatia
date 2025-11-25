test_that("PlotExpression outputs a ggplot object", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  gbm.umap <- data.frame(row.names = colnames(gbm),
                         "UMAP1" = runif(length(colnames(gbm))),
                         "UMAP2" = runif(length(colnames(gbm))))

  reducedDim(gbm, "UMAP") <- gbm.umap
  reducedDim(gbm, "only_UMAP1") <- gbm.umap[, "UMAP1", drop = FALSE]

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    p_vln <- PlotExpression(gbm, transcripts = sample(rownames(gbm_countData), 10))
    p_rd <- PlotExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), plot.type = "reducedDim", dim.use = "UMAP")
    p_heat <- PlotExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), plot.type = "heatmap")
  })

  expect_class(p_vln, "ggplot")
  expect_class(p_rd, "ggplot")
  expect_class(p_heat, "ggplot")

  expect_error(PlotExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), plot.type = "reducedDim", dim.use = "only_UMAP1", verbose = FALSE))

})
