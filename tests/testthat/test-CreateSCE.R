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

make_metadata_alignment_input <- function() {
  counts <- matrix(c(2, 0, 5, 0, 3, 1, 4, 0, 0), nrow = 3, byrow = TRUE,
    dimnames = list(c("tx3", "tx1", "tx2"), c("cell3", "cell1", "cell2")))
  list(
    counts = counts,
    cells = data.frame(group = c("A", "B", "A"), score = c(30L, 10L, 20L),
                       row.names = colnames(counts)),
    transcripts = data.frame(gene_alias = c("g2", "g1", "g2"),
      tx_alias = c("iso3", "iso1", "iso2"), row.names = rownames(counts))
  )
}

test_that("CreateSCE aligns shuffled metadata without reordering counts", {
  input <- make_metadata_alignment_input()
  shuffled_cells <- input$cells[c(2, 3, 1), , drop = FALSE]
  shuffled_transcripts <- input$transcripts[c(3, 1, 2), , drop = FALSE]
  for (counts in list(input$counts, as.data.frame(input$counts),
                      Matrix::Matrix(input$counts, sparse = TRUE))) {
    for (active_transcript in list(NULL, "tx_alias")) {
      expected <- CreateSCE(counts, input$cells, input$transcripts,
        active.gene.id = "gene_alias", active.group.id = "group",
        active.transcript.id = active_transcript, quiet = TRUE)
      expect_message(
        object <- CreateSCE(counts, shuffled_cells, shuffled_transcripts,
          active.gene.id = "gene_alias", active.group.id = "group",
          active.transcript.id = active_transcript, quiet = TRUE), NA
      )
      expect_identical(object, expected)
      expect_identical(rownames(rowData(object)), rownames(input$counts))
      expect_identical(rownames(colData(object)), colnames(input$counts))
      expect_equal(as.matrix(assay(object, "counts")), input$counts)
      expect_s4_class(assay(object, "counts"), "dgCMatrix")
      expect_equal(colData(object)$nCount, c(6, 3, 6))
      expect_identical(colData(object)$nTranscript, c(2L, 1L, 2L))
      expect_identical(colData(object)$nGene, c(1L, 1L, 2L))
      expect_equal(rowData(object)$nCell, c(2, 2, 1))
      expect_identical(rowData(object)$gene_alias, c("g2", "g1", "g2"))
      expect_identical(colData(object)$score, c(30L, 10L, 20L))
    }
  }
  expect_identical(rownames(shuffled_cells), c("cell1", "cell2", "cell3"))
  expect_identical(rownames(shuffled_transcripts), c("tx2", "tx3", "tx1"))
})

test_that("CreateSCE ignores extra metadata rows and preserves single columns", {
  input <- make_metadata_alignment_input()
  cells <- rbind(input$cells, data.frame(group = NA, score = NA_integer_,
                                       row.names = "extra_cell"))
  transcripts <- rbind(input$transcripts,
    data.frame(gene_alias = NA_character_, tx_alias = NA_character_,
               row.names = "extra_tx"))
  object <- CreateSCE(input$counts, cells[c(4, 2, 3, 1), "group", drop = FALSE],
    transcripts[c(3, 4, 1, 2), "gene_alias", drop = FALSE],
    active.gene.id = "gene_alias", active.group.id = "group", quiet = TRUE)
  expect_identical(rownames(object), rownames(input$counts))
  expect_identical(colnames(object), colnames(input$counts))
  expect_identical(rowData(object)$gene_alias, c("g2", "g1", "g2"))
  expect_identical(colData(object)$group, c("A", "B", "A"))
  expect_identical(colData(object)$nGene, c(1L, 1L, 2L))
})

test_that("CreateSCE reports missing metadata IDs even with extra rows", {
  input <- make_metadata_alignment_input()
  cells <- input$cells[1, , drop = FALSE]
  cells <- rbind(cells, data.frame(group = "B", score = 40L, row.names = "extra_cell"))
  expect_error(CreateSCE(input$counts, cells, input$transcripts,
    active.gene.id = "gene_alias", quiet = TRUE),
    "`colData` is missing cell IDs from `countData`: cell1, cell2", fixed = TRUE)

  transcripts <- input$transcripts[2, , drop = FALSE]
  transcripts <- rbind(transcripts,
    data.frame(gene_alias = "g3", tx_alias = "iso4", row.names = "extra_tx"))
  expect_error(CreateSCE(input$counts, input$cells, transcripts,
    active.gene.id = "gene_alias", quiet = TRUE),
    "`rowData` is missing transcript IDs from `countData`: tx3, tx2", fixed = TRUE)
})

test_that("CreateSCE aligns metadata and GTF ranges to the same transcript order", {
  input <- make_metadata_alignment_input()
  gtf <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 101, 201), width = 10),
                               strand = "+")
  S4Vectors::mcols(gtf)$type <- "exon"
  S4Vectors::mcols(gtf)$tx_key <- c("tx1", "tx2", "tx3")
  object <- CreateSCE(input$counts, input$cells[c(3, 1, 2), , drop = FALSE],
    input$transcripts[c(2, 3, 1), , drop = FALSE], active.gene.id = "gene_alias",
    active.transcript.id = "tx_alias", gtf = gtf, gtf.transcript.id = "tx_key",
    quiet = TRUE)
  expect_identical(names(rowRanges(object)), c("tx3", "tx1", "tx2"))
  expect_equal(as.integer(GenomicRanges::start(unlist(rowRanges(object)))), c(201L, 1L, 101L))
  expect_identical(rowData(object)$tx_alias, c("iso3", "iso1", "iso2"))
  expect_identical(rowData(object)$gene_alias, c("g2", "g1", "g2"))
  expect_identical(metadata(object)$GTF, gtf)
})

test_that("CreateSCE rejects nonfinite counts in every supported input format", {
  input <- make_metadata_alignment_input()
  for (value in c(Inf, -Inf, NA_real_, NaN)) {
    counts <- input$counts
    counts[1, 1] <- value
    for (invalid_counts in list(counts, as.data.frame(counts),
                                 Matrix::Matrix(counts, sparse = TRUE))) {
      expect_error(CreateSCE(invalid_counts, input$cells, input$transcripts,
        active.gene.id = "gene_alias", quiet = TRUE))
    }
  }
})

test_that("CreateSCE validates quiet before emitting messages", {
  input <- make_metadata_alignment_input()
  for (quiet in list(NA, 1, c(TRUE, FALSE))) {
    expect_message(
      expect_error(CreateSCE(input$counts, input$cells, input$transcripts,
        active.gene.id = "gene_alias", quiet = quiet), "quiet"), NA
    )
  }
})
