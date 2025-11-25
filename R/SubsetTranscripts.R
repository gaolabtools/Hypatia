#' Subset by transcripts
#'
#' Subset a `SingleCellExperiment` object by transcripts.
#'
#' @param object A `SingleCellExperiment` object.
#' @param subset A logical expression evaluated within the `rowData` of the object.
#' Rows (transcripts) for which the expression evaluates to `TRUE` will be retained. Supports column names from `rowData` and standard logical operators.
#' @param transcripts A vector of transcript IDs for filtering.
#' @param invert Logical; if `TRUE`, inverts the selection to remove the transcripts instead.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns Returns the object after transcript selection.
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

  # Transcript IDs
  assertString(metadata(object)$active.transcript.id)
  active.transcript.id <- metadata(object)$active.transcript.id

  if (metadata(object)$active.transcript.id != "") {
    assertChoice(active.transcript.id, colnames(rowData(object)))
    assertFALSE(any(duplicated(rowData(object)[[active.transcript.id]])))
    assertFALSE(anyMissing(rowData(object)[[active.transcript.id]]))
    rowData(object)$rn_temp_holder <- rownames(object)
    rownames(object) <- rowData(object)[[active.transcript.id]]
  }

  ntx_before <- nrow(object)

  # Subset expression
  if (is.null(transcripts)) {
    rdata <- as.data.frame(rowData(object))
    subset_tx <- rdata %>%
      filter({{ subset }}) %>%
      rownames()

    # Subset
    try_res <- try({
      if (!invert) {
        object <- object[rownames(object) %in% subset_tx, , drop = FALSE]
      }
      if (invert) {
        object <- object[!rownames(object) %in% subset_tx, , drop = FALSE]
      }},
      silent = FALSE
    )

    if (inherits(try_res, "try-error")) {
      stop("The subset term was not valid. Check logical expression or active.transcript.id?")
    }

  }
  # Transcripts were provided
  else if (!is.null(transcripts)) {
    if (!quiet) message("Subsetting using provided transcripts.")

    try_res <- try({
      if (!invert) {
        object <- object[rownames(object) %in% transcripts, ]
      }
      if (invert) {
        object <- object[!rownames(object) %in% transcripts, ]
      }},
      silent = FALSE
    )
  } else if (!is.null(transcripts) && !is.null(subset)) {
    stop("Please provide either the `subset` or `transcripts` parameter.")
  }

  # Change back rownames
  if (metadata(object)$active.transcript.id != "") {
    rownames(object) <- rowData(object)$rn_temp_holder
    rowData(object)$rn_temp_holder <- NULL
  }

  # Update meta cols
  counts <- assay(object, "counts")
  object[["nCount"]] <- colSums(counts)
  object[["nTranscript"]] <- diff(counts@p)

  ## nGene
  gene_ids <- rowData(object)[[metadata(object)$active.gene.id]]
  n_genes_per_cell <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- (counts@p[j] + 1) : counts@p[j + 1]
    if (length(idx) == 0) return(0L)
    transcripts <- counts@i[idx] + 1
    length(unique(gene_ids[transcripts]))
  },
  integer(1)
  )
  object[["nGene"]] <- n_genes_per_cell

  if (!quiet) message("Transcripts before subset: ", ntx_before)
  if (!quiet) message("Transcripts after subset: ", nrow(object))

  return(object)
}
