#' Normalize raw isoform counts
#'
#' Adds a normalized assay to a `SingleCellExperiment` object.
#'
#' @param object A `SingleCellExperiment` object.
#' @param method.use Normalization method: `"LogNormalize"`, `"TPM"`, or `"FT"`. Creates or overwrites `"logcounts"`, `"tpmcounts"`, or `"ftcounts"`, respectively.
#' @param scale.factor Scale factor for `"LogNormalize"` and `"FT"`.
#' @param gtf Optional `GRanges` object with exon annotations for TPM normalization, as described in [CreateSCE()]. Replaces the stored GTF and transcript exon ranges when supplied for TPM.
#' @param gtf.transcript.id Metadata column in `gtf` containing transcript IDs matching `rownames(object)`. If `NULL`, uses `metadata(object)$gtf.transcript.id` when available.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns The input object with the selected normalized assay added.
#' @details
#' - `"LogNormalize"`: Divides counts by the total counts in each cell,
#'   multiplies by `scale.factor`, and applies `log1p()` (natural logarithm).
#'   Results are stored in `"logcounts"`.
#' - `"TPM"`: Divides counts by spliced transcript lengths in kilobases, then scales
#'   the length-adjusted counts to total one million per nonempty cell.
#'   Cells with zero total counts remain zero. Results are stored in
#'   `"tpmcounts"`.
#' - `"FT"`: Divides counts by the total counts in each cell, multiplies by
#'   `scale.factor`, and applies the Freeman-Tukey transformation
#'   `sqrt(x) + sqrt(x + 1)` to stored sparse entries. Implicit sparse zeros
#'   remain zero. Results are stored in `"ftcounts"`.
#'
#' TPM lengths are the summed widths of non-overlapping exon ranges for each
#' transcript, including retained introns when represented within exon records.
#' Transcript biotype labels do not determine which bases are counted.
#' A supplied GTF takes precedence; otherwise existing exon ranges
#' in a `GRangesList` are used, or reconstructed from the stored GTF when absent.
#' A single genomic span in a `GRanges` is insufficient without exon annotations.
#' Exons must have positive widths and share a chromosome and strand within
#' each transcript. The full GTF and its transcript-ID column are stored together
#' in `metadata(object)$GTF` and `metadata(object)$gtf.transcript.id`.
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
  if (is.null(gtf.transcript.id)) {
    gtf.transcript.id <- metadata(object)$gtf.transcript.id
  }
  assertString(gtf.transcript.id, null.ok = TRUE)
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

    # Resolve exon annotations and preserve transcript metadata.
    if (!is.null(gtf)) {
      if (!quiet) message("\u2139 Using user supplied gtf.")
      object <- .StoreTranscriptGTF(object, gtf, gtf.transcript.id)
    } else if (inherits(rowRanges(object), "GRangesList") &&
               all(lengths(rowRanges(object)) > 0L)) {
      annotations <- rowData(object)
      rowRanges(object) <- .ValidateTranscriptExons(rowRanges(object), rownames(object))
      rowData(object) <- annotations
    } else if (!is.null(metadata(object)$GTF)) {
      if (!quiet) message("\u2139 Using exon annotations from the stored GTF.")
      object <- .StoreTranscriptGTF(object, metadata(object)$GTF, gtf.transcript.id)
    } else {
      stop("Transcript exon ranges are required for TPM. Please provide a GTF with exon annotations.", call. = FALSE)
    }

    kb_widths <- sum(GenomicRanges::width(rowRanges(object))) / 1000
    rpk_counts <- (raw_counts / kb_widths)
    rpk_totals <- colSums(rpk_counts)
    tpm_scale <- numeric(length(rpk_totals))
    nonempty <- rpk_totals > 0
    tpm_scale[nonempty] <- 1e6 / rpk_totals[nonempty]
    norm_counts <- rpk_counts %*% Diagonal(x = tpm_scale)
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
