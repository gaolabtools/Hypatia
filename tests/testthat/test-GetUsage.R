test_that("GetUsage outputs a data frame", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error(res <- GetUsage(gbm, genes = c("ENSG00000135945", "ENSG00000048052", "ENSG00000049618")))

  expect_class(res, "data.frame")
  expect_true(nrow(res) >= 3)

})

test_that("GetUsage reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  gbm <- SetTranscripts(gbm, id = "transcript_name")

  res <- GetUsage(
    gbm,
    genes = "ENSG00000135945",
    group.subset = "Tumor",
    cell.dispersion = TRUE,
    quiet = TRUE
  )

  expect_true(all(res$transcript %in% rowData(gbm)$transcript_name))
  expect_true(all(!is.na(res$cell.n)))
})

test_that("GetUsage reports cell-to-cell transcript proportion dispersion", {
  countData <- Matrix::Matrix(
    matrix(
      c(
        8, 0, 4, 0, 3, 1,
        2, 0, 0, 0, 1, 3,
        0, 0, 1, 0, 0, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = rep("gene1", 3), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  disabled <- GetUsage(
    object,
    genes = "gene1",
    min.tx.cts = 2,
    cell.dispersion = FALSE,
    quiet = TRUE
  )
  enabled_messages <- capture.output(
    enabled <- GetUsage(
      object,
      genes = "gene1",
      min.tx.cts = 2,
      cell.dispersion = TRUE,
      quiet = FALSE
    ),
    type = "message"
  )
  expect_true(any(grepl("Calculating cell-to-cell dispersion...", enabled_messages, fixed = TRUE)))

  expect_type(disabled$cell.n, "integer")
  expect_true(all(is.na(disabled$cell.n)))
  expect_true(all(is.na(disabled$cell.prop.mean)))
  legacy_cols <- setdiff(names(disabled), grep("^cell\\.", names(disabled), value = TRUE))
  expect_equal(enabled[legacy_cols], disabled[legacy_cols], tolerance = 1e-12)

  group_a <- enabled[enabled$group == "A", ]
  group_b <- enabled[enabled$group == "B", ]
  expect_equal(group_a$transcript, c("tx1", "tx2"))
  expect_equal(group_a$cell.n, c(2L, 2L))
  expect_equal(group_a$cell.prop.mean, c(0.9, 0.1), tolerance = 1e-12)
  expect_equal(group_a$cell.prop.median, c(0.9, 0.1), tolerance = 1e-12)
  expect_equal(group_a$cell.prop.sd, rep(sd(c(0.8, 1)), 2), tolerance = 1e-12)
  expect_equal(group_a$cell.prop.iqr, c(0.1, 0.1), tolerance = 1e-12)
  expect_equal(group_a$cell.prop.zero.frac, c(0, 0.5), tolerance = 1e-12)
  expect_equal(group_b$cell.n, c(2L, 2L))
  expect_equal(group_b$cell.prop.mean, c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(group_b$cell.prop.sd, rep(sd(c(0.75, 0.25)), 2), tolerance = 1e-12)
  expect_equal(group_b$cell.prop.iqr, c(0.25, 0.25), tolerance = 1e-12)
  expect_equal(group_b$cell.prop.zero.frac, c(0, 0), tolerance = 1e-12)
})
