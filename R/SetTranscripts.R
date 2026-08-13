#' Set active transcript IDs
#'
#' Selects the `rowData` column used as transcript IDs in downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of the `rowData` column containing unique transcript IDs. Use `""` to report row names.
#'
#' @returns The object with `metadata(object)$active.transcript.id` updated.
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
