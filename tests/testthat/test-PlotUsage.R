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

test_that("PlotUsage shares categorical colors without changing the heatmap gradient", {
  for (n in c(3, 12)) {
    counts <- matrix(rep(seq_len(n), 4), nrow = n,
                     dimnames = list(paste0("tx", seq_len(n)), paste0("cell", 1:4)))
    object <- CreateSCE(
      counts,
      data.frame(group = c("A", "A", "B", "B"), row.names = colnames(counts)),
      data.frame(gene_id = rep("gene1", n), row.names = rownames(counts)),
      active.group.id = "group", quiet = TRUE
    )
    custom <- setNames(grDevices::hcl.colors(n, "Dark 3"), rownames(counts))
    for (type in c("stackedbar", "bar", "pie")) {
      for (colors in list(NULL, custom)) {
        plot <- PlotUsage(object, gene = "gene1", plot.type = type,
                          colors = colors, quiet = TRUE)
        scale <- ggplot2::ggplot_build(plot)$plot$scales$get_scales("fill")
        labels <- scale$get_limits()
        expected <- if (is.null(colors)) Hypatia:::.DefaultGroupColors(n) else unname(colors[labels])
        expect_equal(unname(scale$map(labels)), expected)
      }
      plot <- PlotUsage(object, gene = "gene1", plot.type = type,
                        min.tx.prop = 2 / sum(seq_len(n)), quiet = TRUE)
      scale <- ggplot2::ggplot_build(plot)$plot$scales$get_scales("fill")
      expect_true("Other" %in% scale$get_limits())
      expect_equal(scale$map("Other"), "#9E9E9E")
    }
    heatmap <- PlotUsage(object, gene = "gene1", plot.type = "heatmap", quiet = TRUE)
    scale <- ggplot2::ggplot_build(heatmap)$plot$scales$get_scales("fill")
    expect_equal(scale$map(c(0, 0.5, 1)), c("#FFFFFF", "#9C9BE9", "#FFEB3B"))
    expect_equal(scale$get_limits(), c(0, 1))
  }
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
