test_that("RunDEI outputs split data and statistics", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    res <- RunDEI(gbm)
    RunDEI(gbm, group.1 = "Tumor")
    RunDEI(gbm, group.1 = "Tumor", group.2 = "Oligodendrocyte")
    RunDEI(gbm, group.1 = "Tumor", group.2 = c("Oligodendrocyte", "Astrocyte"))
    })


  expect_named(res, c("data", "stats"))
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")
  expect_true(nrow(res$stats) >= 3)

})

test_that("RunDEI respects transcript filters", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  transcripts <- rownames(gbm)[1:2]
  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = transcripts,
    min.pct = 0,
    quiet = TRUE
  )

  expect_setequal(res$stats$transcript, transcripts)
})

test_that("RunDEI requires cells on both sides of a comparison", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  expect_error(
    RunDEI(gbm, group.1 = unique(colData(gbm)$cell_type), quiet = TRUE),
    "at least one cell in both groups"
  )
})

test_that("RunDEI works with active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  gbm <- SetTranscripts(gbm, id = "transcript_name")

  transcripts <- rowData(gbm)$transcript_name[1:2]
  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = transcripts,
    min.pct = 0,
    quiet = TRUE
  )

  expect_setequal(res$stats$transcript, transcripts)
})

test_that("RunDEI reports when no transcripts pass detection thresholds", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  expect_error(
    RunDEI(
      gbm,
      group.1 = "Tumor",
      group.2 = "Oligodendrocyte",
      min.pct = 1,
      quiet = TRUE
    ),
    "0 transcripts passed filtering"
  )
})

test_that("RunDEI supports stats::p.adjust methods", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = rownames(gbm)[1:2],
    min.pct = 0,
    p.adj = "none",
    quiet = TRUE
  )

  expect_equal(res$stats$padj, res$stats$pval)
  expect_error(
    RunDEI(gbm, p.adj = "Bonferroni", quiet = TRUE),
    "element of set"
  )
})

test_that("RunDEI separates summaries and statistics and reports expression dispersion", {
  count_data <- matrix(
    c(0, 2, 4, 4,
      1, 1, 2, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4))
  )
  cell_data <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(count_data)
  )
  transcript_data <- data.frame(
    gene_id = c("gene1", "gene1"),
    row.names = rownames(count_data)
  )
  object <- CreateSCE(
    count_data, cell_data, transcript_data,
    active.group.id = "group", quiet = TRUE
  )
  assay(object, "logcounts") <- count_data

  result_messages <- capture.output(
    result <- suppressWarnings(RunDEI(
      object, group.1 = "A", group.2 = "B", min.pct = 0,
      quiet = FALSE, cell.dispersion = TRUE
    )),
    type = "message"
  )
  expect_true(any(grepl("Calculating cell-to-cell dispersion...", result_messages, fixed = TRUE)))
  disabled <- suppressWarnings(RunDEI(
    object, group.1 = "A", group.2 = "B", min.pct = 0,
    quiet = TRUE, cell.dispersion = FALSE
  ))

  tx1 <- result$data[result$data$transcript == "tx1", ]
  expect_named(
    result,
    c("data", "stats")
  )
  expect_equal(tx1$cell.n.1, 2L)
  expect_equal(tx1$cell.expr.median.1, 1)
  expect_equal(tx1$cell.expr.sd.1, stats::sd(c(0, 2)))
  expect_equal(tx1$cell.expr.iqr.1, stats::IQR(c(0, 2)))
  expect_equal(tx1$cell.expr.median.2, 4)
  expect_equal(tx1$cell.expr.sd.2, 0)
  expect_equal(tx1$cell.expr.iqr.2, 0)
  expect_equal(result$stats, disabled$stats, tolerance = 1e-12)
  legacy_cols <- c(
    "group.1", "group.2", "gene", "transcript", "pct.1", "pct.2",
    "avgExpr.1", "avgExpr.2"
  )
  expect_equal(result$data[legacy_cols], disabled$data[legacy_cols])
  expect_true(all(is.na(disabled$data$cell.n.1)))
  expect_true(all(is.na(disabled$data$cell.expr.sd.2)))
})
