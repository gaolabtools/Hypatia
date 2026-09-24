test_that("GetExpression outputs a data frame", {

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
    res <- GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10))
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = "Tumor")
    GetExpression(gbm, transcripts = sample(rownames(gbm_countData), 10), group.subset = c("Tumor", "Astrocyte"))
  })

  expect_class(res, "data.frame")
  expect_true("pct" %in% names(res))
  expect_false("transcript.pct" %in% names(res))

})

test_that("GetExpression reports cell-to-cell expression dispersion", {
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
    result <- GetExpression(
      object, transcripts = "tx1", quiet = FALSE, cell.dispersion = TRUE
    ),
    type = "message"
  )
  expect_true(any(grepl("Calculating cell-to-cell dispersion...", result_messages, fixed = TRUE)))
  disabled <- GetExpression(
    object, transcripts = "tx1", quiet = TRUE, cell.dispersion = FALSE
  )

  group_a <- result[result$group == "A", ]
  group_b <- result[result$group == "B", ]
  expect_equal(group_a$cell.n, 2L)
  expect_equal(group_a$cell.expr.median, 1)
  expect_equal(group_a$cell.expr.sd, stats::sd(c(0, 2)))
  expect_equal(group_a$cell.expr.iqr, stats::IQR(c(0, 2)))
  expect_equal(group_b$cell.expr.median, 4)
  expect_equal(group_b$cell.expr.sd, 0)
  expect_equal(group_b$cell.expr.iqr, 0)
  expect_equal(result[c("group", "gene", "transcript", "pct", "avgExpr")],
               disabled[c("group", "gene", "transcript", "pct", "avgExpr")])
  expect_true(all(is.na(disabled$cell.n)))
  expect_true(all(is.na(disabled$cell.expr.sd)))
})

test_that("GetExpression reports active transcript IDs", {

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

  transcripts <- rowData(gbm)$transcript_name[1:5]
  res <- GetExpression(gbm, transcripts = transcripts, group.subset = "Tumor", quiet = TRUE)

  expect_setequal(res$transcript, transcripts)
  expect_true(all(res$pct >= 0 & res$pct <= 1))
})

test_that("GetExpression works with gene filters and grouped names", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  colData(gbm)[["cell type"]] <- colData(gbm)$cell_type

  res <- GetExpression(
    gbm,
    genes = "ENSG00000135945",
    group.by = "cell type",
    group.subset = "Tumor",
    quiet = TRUE
  )

  expect_true(all(res$gene == "ENSG00000135945"))
  expect_true("pct" %in% names(res))
})
