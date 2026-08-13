test_that("CreateSCE returns a SingleCellExperiment object", {

  expect_class(
    CreateSCE(
    countData = gbm_countData,
    colData = gbm_colData,
    rowData = gbm_rowData,
    ),
  "SingleCellExperiment")

  expect_class(
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.gene.id = "gene_name",
      active.group.id = "cell_type",
      active.transcript.id = "transcript_name",
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
    ),
    "SingleCellExperiment")

})

test_that("CreateSCE quiet suppresses gtf messages", {

  expect_message(
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
      quiet = TRUE
    ),
    NA
  )
})

test_that("CreateSCE computes QC metrics from detected counts", {

  gbm <- CreateSCE(
    countData = gbm_countData,
    colData = gbm_colData,
    rowData = gbm_rowData,
    quiet = TRUE
  )

  expect_equal(unname(colData(gbm)$nCount), unname(Matrix::colSums(gbm_countData)))
  expect_equal(colData(gbm)$nTranscript, as.integer(Matrix::colSums(gbm_countData > 0)))
  expect_equal(unname(rowData(gbm)$nCell), unname(Matrix::rowSums(gbm_countData > 0)))

  expected_nGene <- vapply(seq_len(25), function(j) {
    detected_tx <- which(gbm_countData[, j] > 0)
    length(unique(gbm_rowData$gene_id[detected_tx]))
  }, integer(1))
  expect_equal(colData(gbm)$nGene[seq_len(25)], expected_nGene)
})

test_that("CreateSCE handles zero-count cells", {

  countData <- Matrix::sparseMatrix(
    i = c(1, 2, 3),
    j = c(1, 3, 3),
    x = c(1, 2, 3),
    dims = c(3, 3),
    dimnames = list(paste0("tx", 1:3), paste0("cell", 1:3))
  )
  colData <- data.frame(group = c("A", "B", "A"), row.names = colnames(countData))
  rowData <- data.frame(gene_id = c("gene1", "gene1", "gene2"), row.names = rownames(countData))

  object <- CreateSCE(countData, colData, rowData, quiet = TRUE)

  expect_equal(colData(object)$nCount, c(1, 0, 5))
  expect_equal(colData(object)$nTranscript, c(1L, 0L, 2L))
  expect_equal(colData(object)$nGene, c(1L, 0L, 2L))
})

test_that("CreateSCE validates count and GTF inputs", {

  countData <- matrix(
    c(1, -1, 0, 2),
    nrow = 2,
    dimnames = list(c("tx1", "tx2"), c("cell1", "cell2"))
  )
  colData <- data.frame(group = c("A", "B"), row.names = colnames(countData))
  rowData <- data.frame(gene_id = c("gene1", "gene2"), row.names = rownames(countData))

  expect_error(CreateSCE(countData, colData, rowData, quiet = TRUE))

  nonnumeric_counts <- as.data.frame(abs(countData))
  nonnumeric_counts[[1]] <- as.character(nonnumeric_counts[[1]])
  expect_error(
    CreateSCE(nonnumeric_counts, colData, rowData, quiet = TRUE),
    "numeric columns"
  )

  expect_error(
    CreateSCE(gbm_countData, gbm_colData, gbm_rowData, gtf = gbm_gtf, quiet = TRUE),
    "gtf.transcript.id"
  )
})
