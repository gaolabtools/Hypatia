test_that("PlotUsage outputs a ggplot object", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(p <- PlotUsage(gbm, gene = "ENSG00000135945", group.by = "cell_type"))

  expect_class(p, "ggplot")

})

test_that("PlotUsage filters transcripts by counts in at least one group", {

  countData <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2),
    j = c(1, 3, 1, 3),
    x = c(5, 7, 1, 1),
    dims = c(3, 4),
    dimnames = list(paste0("tx", 1:3), paste0("cell", 1:4))
  )
  colData <- data.frame(group = c("A", "A", "B", "B"), row.names = colnames(countData))
  rowData <- data.frame(gene_id = rep("gene1", 3), row.names = rownames(countData))
  object <- CreateSCE(countData, colData, rowData, quiet = TRUE)

  p <- PlotUsage(
    object,
    gene = "gene1",
    group.by = "group",
    min.tx.cts = 5,
    quiet = TRUE
  )

  expect_setequal(as.character(p$data$transcripts_query), "tx1")
})

test_that("PlotUsage labels combined group variables consistently", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  p <- PlotUsage(
    gbm,
    gene = "ENSG00000135945",
    group.by = c("cell_type", "project"),
    group.subset = "Tumor_Project",
    quiet = TRUE
  )

  expect_equal(p$labels$x, "cell_type_project")
})
