test_that("the group palette provides distinct deterministic categorical colors", {
  palette <- Hypatia:::.DefaultGroupColors
  expect_identical(palette(0), character())
  expect_identical(palette(3), c("#4E79A7", "#E15759", "#59A14F"))
  expect_identical(palette(6), palette(10)[1:6])
  for (n in c(1, 6, 10, 11, 20, 50)) {
    colors <- palette(n)
    expect_length(colors, n)
    expect_equal(anyDuplicated(colors), 0L)
    expect_identical(colors, palette(n))
    expect_no_error(grDevices::col2rgb(colors))
  }
  expect_identical(palette(12), grDevices::hcl.colors(12, "Dark 3"))
})

test_that("cell-group plots share defaults and respect named custom colors", {
  for (n in c(2, 12)) {
    groups <- sprintf("Group %02d", seq_len(n))
    counts <- matrix(rep(seq_len(17), length.out = 8 * 4 * n), nrow = 8)
    counts[1, seq(1, ncol(counts), by = 2)] <- 0
    rownames(counts) <- paste0("tx", seq_len(nrow(counts)))
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
    sce <- CreateSCE(
      countData = counts,
      colData = data.frame(cell_type = rep(groups, each = 4),
                           row.names = colnames(counts)),
      rowData = data.frame(gene_id = rep(paste0("gene", 1:4), each = 2),
                           row.names = rownames(counts)),
      active.group.id = "cell_type", quiet = TRUE
    )
    custom <- setNames(rev(grDevices::hcl.colors(n, "Dark 3")), rev(groups))
    for (colors in list(NULL, custom)) {
      expected <- if (is.null(colors)) Hypatia:::.DefaultGroupColors(n) else unname(colors[groups])
      qc <- PlotCellQC(sce, group.by = "cell_type", colors = colors, combine = FALSE)
      plots <- c(qc, list(
        lollipop = PlotDiversity(sce, genes = paste0("gene", 1:4),
                                 min.tx.cts = 0, colors = colors, quiet = TRUE),
        density = PlotDiversity(sce, genes = paste0("gene", 1:4),
                                min.tx.cts = 0, plot.type = "density",
                                colors = colors, quiet = TRUE),
        violin = PlotExpression(sce, transcripts = "tx1", assay.use = "counts",
                                colors = colors, quiet = TRUE),
        heatmap = PlotExpression(sce, transcripts = "tx1", assay.use = "counts",
                                 plot.type = "heatmap", colors = colors, quiet = TRUE)
      ))
      for (name in names(plots)) {
        built <- ggplot2::ggplot_build(plots[[name]])
        aesthetic <- if (name %in% c("scatter", "lollipop", "density")) "colour" else "fill"
        scale <- built$plot$scales$get_scales(aesthetic)
        expect_equal(unname(scale$map(groups)), expected, info = paste(n, name))
      }
    }
  }
})

test_that("parallel coordinates use the shared palette for gene colors", {
  for (n in c(2, 12)) {
    genes <- paste0("gene", seq_len(n))
    counts <- matrix(rep(c(9, 1, 6, 4, 2, 8, 3, 7), length.out = n * 2 * 4),
                     nrow = n * 2,
                     dimnames = list(paste0("tx", seq_len(n * 2)), paste0("cell", 1:4)))
    object <- CreateSCE(
      counts,
      data.frame(group = c("A", "A", "B", "B"), row.names = colnames(counts)),
      data.frame(gene_id = rep(genes, each = 2), row.names = rownames(counts)),
      active.group.id = "group", quiet = TRUE
    )
    gene_order <- rev(genes)
    custom <- setNames(grDevices::hcl.colors(n, "Dark 3"), genes)
    for (colors in list(NULL, custom)) {
      plot <- PlotDiversity(object, genes = gene_order, plot.type = "pcoord",
                            group.order = c("B", "A"), colors = colors, quiet = TRUE)
      built <- ggplot2::ggplot_build(plot)
      scale <- built$plot$scales$get_scales("colour")
      expected <- if (is.null(colors)) {
        Hypatia:::.DefaultGroupColors(n)
      } else {
        unname(colors[gene_order])
      }
      expect_identical(unname(scale$map(gene_order)), expected)
      expect_setequal(built$data[[1]]$colour, expected)
      expect_setequal(built$data[[2]]$colour, expected)
    }
  }
})
