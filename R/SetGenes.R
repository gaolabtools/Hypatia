#' Set the active gene ID.
#'
#' Updates the active gene ID (gene symbols, ENSG IDs, etc.) for downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of `rowData` column containing gene IDs associated with each transcript.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns Returns object with updated `active.gene.id` stored in `metadata(object)`.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetGenes <- function (
    object,
    id,
    quiet = FALSE
) {

  assertClass(object, "SingleCellExperiment")
  assertChoice(id, colnames(rowData(object)))
  assertCharacter(rowData(object)[[id]], any.missing = FALSE)
  assertFlag(quiet)

  prev.active.gene.id <- metadata(object)$active.gene.id
  prev.unique.len <- length(unique(rowData(object)[[prev.active.gene.id]]))
  new.unique.len <- length(unique(rowData(object)[[id]]))
  if (prev.unique.len != new.unique.len) {
    if (!quiet) message("\u2139 Warning: The new active.gene.id '", id, "' has different number of unique genes from the previous '", prev.active.gene.id, "' (", new.unique.len, " vs ", prev.unique.len, ").")
  }

  metadata(object)$active.gene.id <- id
  return(object)
}
