test_that("NormalizeCounts returns object with normalized assay", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
      quiet = TRUE
    )

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    gbm <- NormalizeCounts(gbm, method.use = "TPM")
    gbm <- NormalizeCounts(gbm, method.use = "FT")
  })

  expect_true(all(assayNames(gbm) %in% c("counts", "logcounts", "tpmcounts", "ftcounts")))

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- NormalizeCounts(gbm, method.use = "TPM", gtf = gbm_gtf, gtf.transcript.id = "transcript_id"))

  expect_true(all(assayNames(gbm) %in% c("counts", "tpmcounts")))

})

test_that("NormalizeCounts respects quiet for TPM overwrite messages", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
      quiet = TRUE
    )
  gbm <- NormalizeCounts(gbm, method.use = "TPM", quiet = TRUE)

  expect_message(
    NormalizeCounts(gbm, method.use = "TPM", quiet = TRUE),
    NA
  )
  expect_message(
    NormalizeCounts(gbm, method.use = "TPM", quiet = FALSE),
    "overwritten"
  )
})

test_that("TPM uses length-adjusted totals and keeps empty cells at zero", {
  raw_counts <- matrix(
    c(100, 100, 0, 100, 200, 0, 0, 0, 0),
    nrow = 3,
    dimnames = list(c("tx1", "tx2", "tx3"), c("cell1", "cell2", "empty"))
  )
  gtf <- GenomicRanges::GRanges(c("chr1:1-1000", "chr1:1-2000", "chr1:1-500"))
  mcols(gtf)$transcript_id <- rownames(raw_counts)
  # Annotation order must not determine which length is used for a transcript.
  gtf <- gtf[c(3, 1, 2)]
  expected <- matrix(
    c(2e6 / 3, 1e6 / 3, 0, 5e5, 5e5, 0, 0, 0, 0),
    nrow = 3, dimnames = dimnames(raw_counts)
  )

  for (sparse in c(FALSE, TRUE)) {
    input_counts <- if (sparse) Matrix::Matrix(raw_counts, sparse = TRUE) else raw_counts
    object <- CreateSCE(
      input_counts,
      colData = data.frame(group = c("A", "B", "B"), row.names = colnames(raw_counts)),
      rowData = data.frame(gene = c("g1", "g1", "g2"), row.names = rownames(raw_counts)),
      active.gene.id = "gene", quiet = TRUE
    )
    # CreateSCE stores sparse counts; exercise an explicitly dense assay too.
    assay(object, "counts") <- input_counts

    result <- NormalizeCounts(
      object, method.use = "TPM", gtf = gtf,
      gtf.transcript.id = "transcript_id", quiet = TRUE
    )
    tpm <- assay(result, "tpmcounts")
    expect_equal(as.matrix(tpm), expected)
    expect_equal(unname(Matrix::colSums(tpm)), c(1e6, 1e6, 0))
    expect_true(all(is.finite(tpm)))
    expect_identical(assay(result, "counts"), input_counts)
    expect_identical(dimnames(tpm), dimnames(raw_counts))
    if (sparse) expect_s4_class(tpm, "sparseMatrix")

    # The same lengths already stored in rowRanges give the same result.
    repeated <- NormalizeCounts(result, method.use = "TPM", quiet = TRUE)
    expect_equal(as.matrix(assay(repeated, "tpmcounts")), expected)

    # The stored-GTF fallback must use the same normalization too.
    metadata(object)$GTF <- gtf
    from_metadata <- NormalizeCounts(
      object, method.use = "TPM", gtf.transcript.id = "transcript_id", quiet = TRUE
    )
    expect_equal(as.matrix(assay(from_metadata, "tpmcounts")), expected)
  }
})
