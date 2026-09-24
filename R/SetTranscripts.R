#' Set active transcript IDs
#'
#' Selects the `rowData` column used as transcript IDs in downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of the `rowData` column containing unique transcript IDs. Use `""` to report row names.
#'
#' @returns The object with `metadata(object)$active.transcript.id` updated.
#' @details A nonempty `id` must identify a character column of unique,
#' non-missing transcript IDs. Downstream transcript queries and reported
#' transcript labels use that column. Setting `id = ""` restores the use of
#' object row names. This changes the active setting without renaming the
#' stored assay rows or changing transcript-to-gene assignments.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetTranscripts <- function(
    object,
    id
) {

  assertClass(object, "SingleCellExperiment")
  assertString(id)

  if (id != "") {
    assertChoice(id, colnames(rowData(object)))
    assertCharacter(rowData(object)[[id]], any.missing = FALSE)
    assertFALSE(any(duplicated(rowData(object)[[id]])))
  }

  metadata(object)$active.transcript.id <- id

  return(object)
}
