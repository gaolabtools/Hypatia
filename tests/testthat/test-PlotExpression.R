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

  expect_error(
    PlotExpression(
      gbm,
      transcripts = sample(rownames(gbm_countData), 10),
      plot.type = "reducedDim",
      dim.use = "only_UMAP1",
      quiet = TRUE
    ),
    "At least two dimensions"
  )

})

test_that("PlotExpression builds with more groups than default colors", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  colData(gbm)$many_groups <- rep(paste0("group", seq_len(11)), length.out = ncol(gbm))

  p_vln <- PlotExpression(
    gbm,
    transcripts = rownames(gbm)[1],
    group.by = "many_groups",
    quiet = TRUE
  )
  p_heat <- PlotExpression(
    gbm,
    transcripts = rownames(gbm)[1],
    group.by = "many_groups",
    plot.type = "heatmap",
    quiet = TRUE
  )

  expect_no_error(ggplot2::ggplot_build(p_vln))
  expect_no_error(ggplot2::ggplot_build(p_heat))
})

test_that("PlotExpression labels combined group variables consistently", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  p <- PlotExpression(
    gbm,
    transcripts = rownames(gbm)[1],
    group.by = c("cell_type", "project"),
    group.subset = "Tumor_Project",
    quiet = TRUE
  )

  expect_equal(p$labels$x, "cell_type_project")
  expect_equal(p$labels$fill, "cell_type_project")
})
