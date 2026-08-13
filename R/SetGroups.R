#' Set active cell groups
#'
#' Selects the `colData` column used as the default cell grouping in downstream functions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param id Name of the `colData` column containing cell groups.
#' Character and factor columns are supported.
#'
#' @returns The object with `metadata(object)$active.group.id` updated.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-

SetGroups <- function(
    object,
    id
) {

  assertClass(object, "SingleCellExperiment")
  assertString(id)
  assertChoice(id, colnames(colData(object)))

  groups <- colData(object)[[id]]
  if (!is.character(groups) && !is.factor(groups)) {
    stop("The grouping variable must be a character or factor column.")
  }
  assertFALSE(anyMissing(groups))

  if (length(unique(groups)) < 2) {
    stop("The grouping variable must contain at least 2 unique variables.")
  }

  metadata(object)$active.group.id <- id

  return(object)
}
