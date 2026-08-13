#' Subset by transcripts
#'
#' Subsets a `SingleCellExperiment` object by transcript metadata or transcript IDs.
#'
#' @param object A `SingleCellExperiment` object.
#' @param subset Logical expression evaluated in `rowData(object)`. Provide exactly one of `subset` or `transcripts`.
#' @param transcripts Optional vector of transcript IDs. Uses active transcript IDs when `active.transcript.id` is set. Provide exactly one of `subset` or `transcripts`.
#' @param invert Logical; if `TRUE`, remove the selected transcripts instead of keeping them.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns The object after transcript selection, with cell QC columns updated.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix colSums

SubsetTranscripts <- function(
    object,
    subset = NULL,
    transcripts = NULL,
    invert = FALSE,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertCharacter(transcripts, null.ok = TRUE, any.missing = FALSE)
  assertFlag(invert)
  assertFlag(quiet)

  subset_expr <- substitute(subset)
  has_subset <- !missing(subset) && !identical(subset_expr, quote(NULL))
  has_transcripts <- !is.null(transcripts)
  if (has_subset == has_transcripts) {
    stop("Please provide exactly one of the `subset` or `transcripts` parameters.")
  }

  # Transcript IDs
  assertString(metadata(object)$active.transcript.id)
  active.transcript.id <- metadata(object)$active.transcript.id
  use_active_transcript_id <- active.transcript.id != ""

  if (use_active_transcript_id) {
    assertChoice(active.transcript.id, colnames(rowData(object)))
    assertFALSE(any(duplicated(rowData(object)[[active.transcript.id]])))
    assertFALSE(anyMissing(rowData(object)[[active.transcript.id]]))

    active_transcripts <- as.character(rowData(object)[[active.transcript.id]])
    original_rownames <- rownames(object)
    names(original_rownames) <- active_transcripts
    rownames(object) <- active_transcripts
  }

  ntx_before <- nrow(object)

  # Subset expression
  if (has_subset) {
    try_res <- try(
      {
        rdata <- as.data.frame(rowData(object))
        names(rdata) <- names(rowData(object))
        subset_tx <- rdata %>%
          filter({{ subset }}) %>%
          rownames()

        if (!invert) {
          object <- object[rownames(object) %in% subset_tx, , drop = FALSE]
        }
        if (invert) {
          object <- object[!rownames(object) %in% subset_tx, , drop = FALSE]
        }
      },
      silent = TRUE
    )

    if (inherits(try_res, "try-error")) {
      stop("The subset term was not valid. Check logical expression or active.transcript.id?", call. = FALSE)
    }
  }
  # Transcripts were provided
  else if (has_transcripts) {
    if (!quiet) message("Subsetting using provided transcripts.")

    try_res <- try(
      {
        if (!invert) {
          object <- object[rownames(object) %in% transcripts, , drop = FALSE]
        }
        if (invert) {
          object <- object[!rownames(object) %in% transcripts, , drop = FALSE]
        }
      },
      silent = TRUE
    )

    if (inherits(try_res, "try-error")) {
      stop("The provided transcript IDs were not valid.", call. = FALSE)
    }
  }

  # Change back rownames
  if (use_active_transcript_id) {
    rownames(object) <- unname(original_rownames[rownames(object)])
  }

  # Update meta cols
  counts <- assay(object, "counts")
  object[["nCount"]] <- colSums(counts)
  object[["nTranscript"]] <- as.integer(colSums(counts > 0))
  gene_ids <- rowData(object)[[metadata(object)$active.gene.id]]
  object[["nGene"]] <- .NGenesPerCell(counts, gene_ids)

  if (!quiet) message("Transcripts before subset: ", ntx_before)
  if (!quiet) message("Transcripts after subset: ", nrow(object))

  return(object)
}
