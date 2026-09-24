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
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Shannon")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Shannon", top.n = 2)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "NormalizedShannon")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Renyi")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Tsallis", order = 1)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Renyi", order = 1)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "NormalizedRenyi", order = 1)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "Renyi", top.n = 2)
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "NormalizedRenyi")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "GiniSimpson")
    PlotDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:3], min.tx.cts = 0, entropy.use = "InverseSimpson")
  })

  expect_class(p_lp, "ggplot")
  expect_class(p_den, "ggplot")
  expect_class(p_pc, "ggplot")
  expect_s3_class(p_den$facet, "FacetNull")
  expect_false("nrow" %in% names(formals(PlotDiversity)))
  expect_equal(p_den$labels$colour, "cell_type")
  expect_equal(all.vars(p_den$layers[[1]]$mapping$colour), "group_var")
  expect_null(p_den$layers[[1]]$mapping$fill)
  expect_equal(p_den$layers[[1]]$aes_params$linewidth, 1)
  p_lp_data <- ggplot_build(p_lp)$data[[2]]
  expect_equal(p_lp$labels$shape, "class")
  expect_setequal(unique(stats::na.omit(p_lp_data$shape)), c(19, 21))
  expect_true(all(p_lp_data$fill == "white"))
  expect_equal(p_lp$guides$guides$colour$params$order, 1)
  expect_equal(p_lp$guides$guides$shape$params$order, 2)
  expect_equal(p_pc$labels$shape, "class")
  expect_equal(p_pc$guides$guides$colour$params$order, 1)
  expect_equal(p_pc$guides$guides$shape$params$order, 2)

})

test_that("PlotDiversity labels combined group variables consistently", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  p <- PlotDiversity(
    gbm,
    genes = unique(gbm_rowData$gene_id)[1:3],
    group.by = c("cell_type", "project"),
    group.subset = "Tumor_Project",
    min.tx.cts = 0,
    quiet = TRUE
  )

  expect_equal(p$labels$colour, "cell_type_project")
})
