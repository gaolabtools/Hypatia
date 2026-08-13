#' Set active gene IDs
#'
#' Selects the `rowData` column used as gene IDs in downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of the `rowData` column containing gene IDs.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns The object with `metadata(object)$active.gene.id` updated.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetGenes <- function(
    object,
    id,
    quiet = FALSE
) {

  assertClass(object, "SingleCellExperiment")
  assertChoice(id, colnames(rowData(object)))
  assertCharacter(rowData(object)[[id]], any.missing = FALSE)
  assertFlag(quiet)

  prev.active.gene.id <- metadata(object)$active.gene.id
  assertString(prev.active.gene.id)

  prev.unique.len <- length(unique(rowData(object)[[prev.active.gene.id]]))
  new.unique.len <- length(unique(rowData(object)[[id]]))
  if (prev.unique.len != new.unique.len) {
    if (!quiet) {
      message(
        "\u2139 Warning: The new active.gene.id '",
        id,
        "' has different number of unique genes from the previous '",
        prev.active.gene.id,
        "' (",
        new.unique.len,
        " vs ",
        prev.unique.len,
        ")."
      )
    }
  }

  metadata(object)$active.gene.id <- id

  return(object)
}
