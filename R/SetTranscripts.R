#' Set the active transcript ID
#'
#' Update the active transcript ID for downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of `rowData` column containing unique transcript IDs.
#' Set to `""` to use the transcript IDs (row names) of original count matrix.
#'
#' @returns Returns object with updated `active.transcript.id` stored in `metadata(object)`.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetTranscripts <- function (
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
