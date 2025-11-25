#' Set the active cell group
#'
#' Sets or updates the active cell group for downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of `colData` column containing the desired cell groups (cell types, cell states, sample ID, etc.).
#'
#' @returns Returns object with updated `active.group.id` stored in `metadata(object)`.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetGroups <- function (
    object,
    id
) {

  assertClass(object, "SingleCellExperiment")
  assertString(id)
  assertChoice(id, colnames(colData(object)))
  assertCharacter(colData(object)[[id]], any.missing = FALSE)

  if (length(unique(colData(object)[[id]])) < 2) {
    stop("The grouping variable must contain at least 2 unique variables.")
  }

  metadata(object)$active.group.id <- id
  return(object)
}
