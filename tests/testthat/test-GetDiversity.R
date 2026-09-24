test_that("GetDiversity works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    res <- GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Shannon", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Shannon", top.n = 2, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "NormalizedShannon", prop.thresh = 0.2, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Renyi", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Tsallis", order = 1, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Renyi", order = 1, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "NormalizedRenyi", order = 1, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "Renyi", top.n = 2, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "NormalizedRenyi", prop.thresh = 0.2, quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "GiniSimpson", quiet = TRUE)
    GetDiversity(gbm, genes = unique(gbm_rowData$gene_id)[1:5], min.tx.cts = 0, entropy.use = "InverseSimpson", quiet = TRUE)
  })

  expect_class(res, "data.frame")

})

test_that("GetDiversity reports active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  gbm <- SetTranscripts(gbm, id = "transcript_name")

  res <- GetDiversity(
    gbm,
    genes = "ENSG00000135945",
    group.subset = "Tumor",
    min.tx.cts = 0,
    cell.dispersion = TRUE,
    quiet = TRUE
  )

  expect_true(all(res$transcript %in% rowData(gbm)$transcript_name))
  expect_true(all(!is.na(res$cell.n)))
})

test_that("GetDiversity reports biological cell-to-cell diversity dispersion", {
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

  disabled <- GetDiversity(
    object,
    genes = "gene1",
    min.tx.cts = 2,
    cell.dispersion = FALSE,
    quiet = TRUE
  )
  enabled_messages <- capture.output(
    enabled <- GetDiversity(
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
  expect_true(all(is.na(disabled$cell.div.mean)))
  legacy_cols <- setdiff(names(disabled), grep("^cell\\.", names(disabled), value = TRUE))
  expect_equal(enabled[legacy_cols], disabled[legacy_cols], tolerance = 1e-12)

  group_a <- enabled[enabled$group == "A", ]
  group_b <- enabled[enabled$group == "B", ]
  div_a <- c(0.24, 0)
  div_b <- c(0.28125, 0.28125)
  expect_equal(group_a$cell.n, c(2L, 2L))
  expect_equal(group_a$cell.div.mean, rep(mean(div_a), 2), tolerance = 1e-12)
  expect_equal(group_a$cell.div.median, rep(median(div_a), 2), tolerance = 1e-12)
  expect_equal(group_a$cell.div.sd, rep(sd(div_a), 2), tolerance = 1e-12)
  expect_equal(group_a$cell.div.iqr, rep(IQR(div_a), 2), tolerance = 1e-12)
  expect_equal(group_b$cell.n, c(2L, 2L))
  expect_equal(group_b$cell.div.mean, rep(mean(div_b), 2), tolerance = 1e-12)
  expect_equal(group_b$cell.div.median, rep(median(div_b), 2), tolerance = 1e-12)
  expect_equal(group_b$cell.div.sd, c(0, 0), tolerance = 1e-12)
  expect_equal(group_b$cell.div.iqr, c(0, 0), tolerance = 1e-12)
})
