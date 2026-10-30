test_that("SubsetTranscripts works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error({
    gbm_nCell <- SubsetTranscripts(gbm, subset = nCell >= 3)
    gbm_nCell_inv <- SubsetTranscripts(gbm, subset = nCell >= 3, invert = TRUE)
    gbm_sq <- SubsetTranscripts(gbm, subset = structural_category == "full-splice_match")
    gbm_feat <- SubsetTranscripts(gbm, transcripts = rownames(gbm)[1:5])
  })

  expect_true({
    all(rowData(gbm_nCell)$nCell >= 3)
    all(rowData(gbm_nCell_inv)$nCell < 3)
    all(rowData(gbm_sq)$structural_category == "full-splice_match")
    all(rownames(gbm_feat) == rownames(gbm)[1:5])
  })

})

test_that("SubsetTranscripts requires exactly one selector", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_error(SubsetTranscripts(gbm, quiet = TRUE), "exactly one")
  expect_error(
    SubsetTranscripts(gbm, subset = nCell >= 3, transcripts = rownames(gbm)[1:5], quiet = TRUE),
    "exactly one"
  )
})

test_that("SubsetTranscripts supports non-syntactic rowData names", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)[["structural category"]] <- rowData(gbm)$structural_category

  expect_no_error({
    gbm_sq <- SubsetTranscripts(
      gbm,
      subset = `structural category` == "full-splice_match",
      quiet = TRUE
    )
  })
  expect_true(all(rowData(gbm_sq)[["structural category"]] == "full-splice_match"))
})

test_that("SubsetTranscripts preserves rowData columns when active transcript IDs are used", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)$rn_temp_holder <- "keep_me"
  gbm <- SetTranscripts(gbm, id = "transcript_name")
  selected_transcripts <- rowData(gbm)$transcript_name[1:5]

  gbm_subset <- SubsetTranscripts(
    gbm,
    transcripts = selected_transcripts,
    quiet = TRUE
  )

  expect_equal(rownames(gbm_subset), rownames(gbm)[1:5])
  expect_true("rn_temp_holder" %in% names(rowData(gbm_subset)))
  expect_true(all(rowData(gbm_subset)$rn_temp_holder == "keep_me"))
})

test_that("SubsetTranscripts recomputes cell QC metrics after transcript filtering", {

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
  object <- SubsetTranscripts(object, transcripts = c("tx2", "tx3"), quiet = TRUE)

  expect_equal(unname(colData(object)$nCount), c(0, 0, 5))
  expect_equal(colData(object)$nTranscript, c(0L, 0L, 2L))
  expect_equal(colData(object)$nGene, c(0L, 0L, 2L))
})

test_that("SubsetTranscripts reports invalid subset expressions", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  error <- expect_error(
    SubsetTranscripts(gbm, subset = missing_column == "x", quiet = TRUE),
    "Failed to subset transcripts using `subset`:", fixed = TRUE
  )
  expect_match(conditionMessage(error), "missing_column", fixed = TRUE)
})

make_subset_metadata_sce <- function(active_ids = FALSE) {
  counts <- matrix(c(2, 0, 5, 0, 3, 1, 4, 0, 0), nrow = 3, byrow = TRUE,
    dimnames = list(c("tx3", "tx1", "tx2"), c("cell3", "cell1", "cell2")))
  CreateSCE(counts,
    data.frame(group = c("A", "B", "A"), row.names = colnames(counts)),
    data.frame(gene_alias = c("g2", "g1", "g2"),
      tx_alias = c("iso3", "iso1", "iso2"), transcript_strand = c("+", "-", "+"),
      row.names = rownames(counts)),
    active.gene.id = "gene_alias",
    active.transcript.id = if (active_ids) "tx_alias" else NULL, quiet = TRUE)
}

for (reserved_name in c("seqnames", "ranges", "strand", "start", "end", "width", "element")) {
  test_that(paste("SubsetTranscripts retains errors for legacy column", reserved_name), {
    object <- make_subset_metadata_sce()
    # Simulate objects created before CreateSCE rejected reserved names.
    rowData(object)[[reserved_name]] <- c("+", "-", "+")
    original_error <- tryCatch(object[1:2, , drop = FALSE], error = identity)
    expect_s3_class(original_error, "error")
    expect_error(SubsetTranscripts(object, transcripts = "tx1", quiet = TRUE),
                 conditionMessage(original_error), fixed = TRUE)
    expect_error(SubsetTranscripts(object, subset = nCell >= 1, quiet = TRUE),
                 conditionMessage(original_error), fixed = TRUE)
  })
}

test_that("valid transcript subsetting preserves row IDs, annotations, and cell order", {
  for (active_ids in c(FALSE, TRUE)) {
    object <- make_subset_metadata_sce(active_ids)
    selected <- if (active_ids) c("iso2", "iso3") else c("tx2", "tx3")
    for (invert in c(FALSE, TRUE)) {
      expected_rows <- if (invert) 2L else c(1L, 3L)
      by_ids <- SubsetTranscripts(object, transcripts = selected, invert = invert,
                                 quiet = TRUE)
      by_expression <- SubsetTranscripts(object, subset = transcript_strand == "+",
                                        invert = invert, quiet = TRUE)
      expect_identical(by_ids, by_expression)
      expect_identical(rownames(by_ids), rownames(object)[expected_rows])
      expect_identical(colnames(by_ids), colnames(object))
      expect_identical(rowData(by_ids), rowData(object)[expected_rows, , drop = FALSE])
      expect_identical(assay(by_ids, "counts"),
                       assay(object, "counts")[expected_rows, , drop = FALSE])
      expect_true(methods::validObject(by_ids))
    }
  }
})
