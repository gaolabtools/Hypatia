.gtf_example <- function() {
  gtf <- GenomicRanges::GRanges(c(
    "chr1:1-2400:+", "chr1:1-1100:+", "chr1:1-100:+",
    "chr1:1001-1100:+", "chr1:2001-2400:-", "chr1:2001-2400:-"
  ))
  mcols(gtf)$type <- c("gene", "transcript", "exon", "exon", "transcript", "exon")
  mcols(gtf)$tx_key <- c(NA, "txA", "txA", "txA", "txB", "txB")
  counts <- matrix(
    c(100, 100, 200, 100, 0, 0), nrow = 2,
    dimnames = list(c("txB", "txA"), c("cell1", "cell2", "empty"))
  )
  list(
    gtf = gtf, counts = counts,
    cells = data.frame(group = c("A", "B", "B"), row.names = colnames(counts)),
    transcripts = data.frame(
      gene = c("gB", "gA"), label = c("B", "A"), row.names = rownames(counts)
    )
  )
}

.gtf_object <- function(example, with_gtf = TRUE) {
  CreateSCE(
    example$counts, example$cells, example$transcripts,
    active.gene.id = "gene", active.transcript.id = "label",
    gtf = if (with_gtf) example$gtf else NULL,
    gtf.transcript.id = if (with_gtf) "tx_key" else NULL,
    quiet = TRUE
  )
}

test_that("GTF exon ranges match assay rows and preserve annotations", {
  example <- .gtf_example()
  object <- .gtf_object(example)
  expect_s4_class(rowRanges(object), "GRangesList")
  expect_identical(names(rowRanges(object)), rownames(example$counts))
  expect_equal(unname(lengths(rowRanges(object))), c(1L, 2L))
  expect_equal(unname(sum(GenomicRanges::width(rowRanges(object)))), c(400L, 200L))
  expect_identical(metadata(object)$GTF, example$gtf)
  expect_identical(metadata(object)$gtf.transcript.id, "tx_key")
  expect_identical(as.data.frame(rowData(object))[, c("gene", "label")], example$transcripts)
  expect_equal(as.matrix(assay(object, "counts")), example$counts)

  reordered <- example
  reordered$gtf <- example$gtf[c(6, 4, 1, 3, 5, 2)]
  expect_identical(rowRanges(.gtf_object(reordered)), rowRanges(object))
})

test_that("TPM uses spliced exon lengths for dense and sparse counts", {
  example <- .gtf_example()
  expected <- matrix(
    c(1e6 / 3, 2e6 / 3, 5e5, 5e5, 0, 0), nrow = 2,
    dimnames = dimnames(example$counts)
  )
  for (sparse in c(FALSE, TRUE)) {
    object <- .gtf_object(example)
    input <- if (sparse) Matrix::Matrix(example$counts, sparse = TRUE) else example$counts
    assay(object, "counts") <- input
    result <- NormalizeCounts(object, method.use = "TPM", quiet = TRUE)
    expect_equal(as.matrix(assay(result, "tpmcounts")), expected)
    expect_equal(unname(Matrix::colSums(assay(result, "tpmcounts"))), c(1e6, 1e6, 0))
    expect_identical(assay(result, "counts"), input)
    expect_identical(rowData(result), rowData(object))
    if (sparse) expect_s4_class(assay(result, "tpmcounts"), "sparseMatrix")
  }
})

test_that("GTF replacement updates stored annotations and exon ranges together", {
  example <- .gtf_example()
  object <- .gtf_object(example)
  replacement <- example$gtf
  replacement[6] <- GenomicRanges::resize(replacement[6], width = 200, fix = "start")
  mcols(replacement)$replacement_id <- mcols(replacement)$tx_key
  mcols(replacement)$tx_key <- NULL

  result <- NormalizeCounts(
    object, method.use = "TPM", gtf = replacement,
    gtf.transcript.id = "replacement_id", quiet = TRUE
  )
  expect_identical(metadata(result)$GTF, replacement)
  expect_identical(metadata(result)$gtf.transcript.id, "replacement_id")
  expect_equal(unname(sum(GenomicRanges::width(rowRanges(result)))), c(200L, 200L))
  expect_equal(unname(as.matrix(assay(result, "tpmcounts"))[, 1]), c(5e5, 5e5))
  expect_identical(rowData(result), rowData(object))

  # Reusing the remembered ID column also works for an explicitly supplied GTF.
  repeated <- NormalizeCounts(result, method.use = "TPM", gtf = replacement, quiet = TRUE)
  expect_equal(assay(repeated, "tpmcounts"), assay(result, "tpmcounts"))
})

test_that("stored GTFs can rebuild exon ranges using the remembered ID column", {
  example <- .gtf_example()
  object <- .gtf_object(example, with_gtf = FALSE)
  metadata(object)$GTF <- example$gtf
  metadata(object)$gtf.transcript.id <- "tx_key"
  result <- NormalizeCounts(object, method.use = "TPM", quiet = TRUE)
  expect_identical(rowRanges(result), rowRanges(.gtf_object(example)))

  # A single genomic span per transcript is reconstructed from the stored GTF.
  annotations <- rowData(object)
  spans <- example$gtf[c(5, 2)]
  names(spans) <- rownames(object)
  rowRanges(object) <- spans
  rowData(object) <- annotations
  rebuilt <- NormalizeCounts(object, method.use = "TPM", quiet = TRUE)
  expect_equal(assay(rebuilt, "tpmcounts"), assay(result, "tpmcounts"))
  expect_identical(rowData(rebuilt), annotations)

  metadata(object)$gtf.transcript.id <- NULL
  expect_error(NormalizeCounts(object, method.use = "TPM", quiet = TRUE), "gtf.transcript.id")
  expect_no_error(NormalizeCounts(
    object, method.use = "TPM", gtf.transcript.id = "tx_key", quiet = TRUE
  ))
})

test_that("exon ranges and TPM stay aligned after active-ID transcript subsetting", {
  example <- .gtf_example()
  object <- .gtf_object(example)
  selected <- SubsetTranscripts(object, transcripts = "A", quiet = TRUE)
  expect_identical(rownames(selected), "txA")
  expect_identical(names(rowRanges(selected)), "txA")
  expect_equal(unname(lengths(rowRanges(selected))), 2L)
  expect_equal(unname(sum(GenomicRanges::width(rowRanges(selected)))), 200L)
  expect_identical(metadata(selected)$GTF, example$gtf)
  result <- NormalizeCounts(selected, method.use = "TPM", quiet = TRUE)
  expect_equal(unname(as.matrix(assay(result, "tpmcounts"))), matrix(c(1e6, 1e6, 0), nrow = 1))

  selected <- SubsetTranscripts(object, transcripts = character(), quiet = TRUE)
  expect_no_error(NormalizeCounts(selected, method.use = "TPM", quiet = TRUE))
})

test_that("overlapping or repeated exon records do not double-count bases", {
  example <- .gtf_example()
  expected <- .gtf_object(example)
  example$gtf <- c(example$gtf, example$gtf[3], GenomicRanges::resize(example$gtf[3], width = 50))
  expect_identical(rowRanges(.gtf_object(example)), rowRanges(expected))
})

test_that("GTF validation requires exons for each requested transcript", {
  example <- .gtf_example()
  missing_exons <- example
  missing_exons$gtf <- example$gtf[-6]
  expect_error(.gtf_object(missing_exons), "exon.*txB")

  transcript_only <- example
  transcript_only$gtf <- example$gtf[c(2, 5)]
  expect_error(.gtf_object(transcript_only), "exon")

  missing_ids <- example
  mcols(missing_ids$gtf)$tx_key[3] <- NA_character_
  expect_error(.gtf_object(missing_ids), "transcript IDs")

  zero_width <- example
  zero_width$gtf[3] <- GenomicRanges::resize(zero_width$gtf[3], width = 0)
  expect_error(.gtf_object(zero_width), "positive")

  inconsistent <- example
  GenomicRanges::strand(inconsistent$gtf)[3] <- "-"
  expect_error(.gtf_object(inconsistent), "chromosome and strand")
})

test_that("TPM validates external exon lists and refuses unannotated genomic spans", {
  example <- .gtf_example()
  object <- .gtf_object(example)
  metadata(object)$GTF <- NULL
  expect_no_error(NormalizeCounts(object, method.use = "TPM", quiet = TRUE))

  annotations <- rowData(object)
  spans <- example$gtf[c(5, 2)]
  names(spans) <- rownames(object)
  rowRanges(object) <- spans
  rowData(object) <- annotations
  expect_error(NormalizeCounts(object, method.use = "TPM", quiet = TRUE), "exon")

  rowRanges(object) <- GenomicRanges::GRangesList(txB = example$gtf[6], txA = example$gtf[0])
  rowData(object) <- annotations
  expect_error(NormalizeCounts(object, method.use = "TPM", quiet = TRUE), "exon")
})

test_that("TPM includes retained introns represented in transcript exon coordinates", {
  example <- .gtf_example()
  example$gtf <- GenomicRanges::GRanges(c(
    "chr1:1-300:+", "chr1:1-100:+", "chr1:201-300:+",
    "chr1:1-300:+", "chr1:1-300:+"
  ))
  mcols(example$gtf)$type <- c("transcript", "exon", "exon", "transcript", "exon")
  mcols(example$gtf)$tx_key <- c("txA", "txA", "txA", "txB", "txB")
  object <- .gtf_object(example)
  expect_equal(unname(sum(GenomicRanges::width(rowRanges(object)))), c(300L, 200L))
  result <- NormalizeCounts(object, method.use = "TPM", quiet = TRUE)
  # Equal counts: 100 / 0.3 kb versus 100 / 0.2 kb gives a 40:60 TPM split.
  expect_equal(unname(as.matrix(assay(result, "tpmcounts"))[, 1]), c(4e5, 6e5))
})
