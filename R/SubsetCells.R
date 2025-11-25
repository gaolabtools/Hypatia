#' Subset by cells
#'
#' Subsets a `SingleCellExperiment` object by cells.
#'
#' @param object A `SingleCellExperiment` object.
#' @param subset A logical expression evaluated within the `colData` of the object.
#' Rows (cells) for which the expression evaluates to `TRUE` will be retained. Supports column names from `colData` and standard logical operators.
#' @param cells A vector of cell IDs for filtering.
#' @param invert Logical; if `TRUE`, inverts the selection to remove the cells instead.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns Returns the object after cell selection.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom Matrix rowSums

SubsetCells <- function(
    object,
    subset = NULL,
    cells = NULL,
    invert = FALSE,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertCharacter(cells, null.ok = TRUE, any.missing = FALSE)
  assertFlag(invert)
  assertFlag(quiet)

  ncell_before <- ncol(object)

  # Subset expression
  if (is.null(cells)) {
    cdata <- as.data.frame(colData(object))
    subset_cells <- cdata %>%
      filter({{ subset }}) %>%
      rownames()

    # Subset
    try_res <- try({
      if (!invert) {
        object <- object[, colnames(object) %in% subset_cells, drop = FALSE]
      }
      if (invert) {
        object <- object[, !colnames(object) %in% subset_cells, drop = FALSE ]
      }},
      silent = FALSE
    )

    if (inherits(try_res, "try-error")) {
      stop("The subset term was not valid. Check logical expression?")
    }

  }
  # Cell IDs were provided
  else if (!is.null(cells)) {
    if (!quiet) message("Subsetting using provided cell IDs.")

    try_res <- try({
      if (!invert) {
        object <- object[, colnames(object) %in% cells, drop = FALSE]
      }
      if (invert) {
        object <- object[, !colnames(object) %in% cells, drop = FALSE]
      }},
      silent = FALSE
    )
  } else if (!is.null(cells) && !is.null(subset)) {
    stop("Please provide either the `subset` or `cells` parameter.")
  }

  # Update meta cols
  counts <- assay(object, "counts")
  rowData(object)[["nCell"]] <- rowSums(counts > 0)

  if (!quiet) message("Cells before subset: ", ncell_before)
  if (!quiet) message("Cells after subset: ", ncol(object))

  return(object)
}
