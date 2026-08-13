#' Normalize raw isoform counts
#'
#' Adds a normalized assay to a `SingleCellExperiment` object.
#'
#' @param object A `SingleCellExperiment` object.
#' @param method.use Normalization method: `"LogNormalize"`, `"TPM"`, or `"FT"`. Creates or overwrites `"logcounts"`, `"tpmcounts"`, or `"ftcounts"`, respectively.
#' @param scale.factor Scale factor for `"LogNormalize"` and `"FT"`.
#' @param gtf Optional `GRanges` object with transcript ranges for TPM normalization.
#' @param gtf.transcript.id Metadata column in `gtf` containing transcript IDs matching `rownames(object)`.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns The input object with the selected normalized assay added.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import Matrix
#' @importFrom S4Vectors mcols mcols<-

NormalizeCounts <- function(
    object,
    method.use = "LogNormalize",
    scale.factor = 10000,
    gtf = NULL,
    gtf.transcript.id = NULL,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertChoice(method.use, c("LogNormalize", "TPM", "FT"))
  assertNumber(scale.factor, lower = 1, finite = TRUE)
  if (!is.null(gtf)) {
    assertClass(gtf, "GRanges")
    assertString(gtf.transcript.id, null.ok = FALSE)
    assertTRUE(gtf.transcript.id %in% names(mcols(gtf)))
    assertTRUE(all(rownames(object) %in% mcols(gtf)[[gtf.transcript.id]]))
  }
  assertFlag(quiet)

  # Check if selected assay exists
  current_assays <- assayNames(object)
  if (method.use == "LogNormalize" && "logcounts" %in% current_assays) {
    if (!quiet) message("\u2139 Warning: The assay 'logcounts' already exists and will be overwritten.")
  }
  if (method.use == "TPM" && "tpmcounts" %in% current_assays) {
    if (!quiet) message("\u2139 Warning: The assay 'tpmcounts' already exists and will be overwritten.")
  }

  # Normalization
  raw_counts <- counts(object)
  col_sum <- colSums(raw_counts)

  ## LogNormalize
  if (method.use == "LogNormalize") {
    if (!quiet) message("Performing log normalization...")

    norm_counts <- raw_counts %*% Diagonal(x = scale.factor / col_sum)
    norm_counts <- log1p(norm_counts)
    dimnames(norm_counts) <- dimnames(raw_counts)

    assay(object, "logcounts") <- norm_counts
    if (!quiet) message("Done.")
  }

  ## TPM
  if (method.use == "TPM") {
    if (!quiet) message("Performing TPM normalization...")

    # if gtf is supplied, add rowRanges
    if (!is.null(gtf)) {

      if (!quiet) message("\u2139 Using user supplied gtf.")

      # add gtf to metadata
      if (is.null(object@metadata$GTF)) {
        object@metadata$GTF <- gtf
      }

      # add rowRanges
      rowData <- rowData(object)
      names(gtf) <- mcols(gtf)[[gtf.transcript.id]]
      gtf <- gtf[rownames(object)]
      mcols(gtf) <- NULL
      rowRanges(object) <- gtf
      rowData(object) <- rowData
    }

    # if rowRanges exists, skip
    else if (all(lengths(rowRanges(object)) != 0)) {
    }

    # if GTF metadata exists and rowRanges does not, add rowRanges
    else if (!is.null(object@metadata$GTF)) {
      if (!quiet) message("\u2139 GTF metadata detected. Using transcript widths.")
      assertClass(object@metadata$GTF, "GRanges")
      assertTRUE(gtf.transcript.id %in% names(mcols(object@metadata$GTF)))
      assertTRUE(all(rownames(object) %in% mcols(object@metadata$GTF)[[gtf.transcript.id]]))

      # add rowRanges
      rowData <- rowData(object)
      gtf <- object@metadata$GTF
      names(gtf) <- mcols(gtf)[[gtf.transcript.id]]
      gtf <- gtf[rownames(object)]
      mcols(gtf) <- NULL
      rowRanges(object) <- gtf
      rowData(object) <- rowData

    } else {
      stop("Transcript widths not found in rowRanges or GTF metadata. Please provide the gtf.")
    }

    kb_widths <- width(rowRanges(object)) / 1000
    rpk_counts <- (raw_counts / kb_widths)
    scale.factor <- col_sum / 1e6
    norm_counts <- rpk_counts %*% Diagonal(x = 1 / scale.factor)
    dimnames(norm_counts) <- dimnames(raw_counts)

    assay(object, "tpmcounts") <- norm_counts
    if (!quiet) message("Done.")
  }

  if (method.use == "FT") {
    if (!quiet) message("Performing Freeman-Tukey normalization...")

    norm_counts <- raw_counts %*% Diagonal(x = scale.factor / col_sum)
    norm_counts@x <- sqrt(norm_counts@x) + sqrt(norm_counts@x + 1)
    dimnames(norm_counts) <- dimnames(raw_counts)

    assay(object, "ftcounts") <- norm_counts
    if (!quiet) message("Done.")
  }

  return(object)
}
