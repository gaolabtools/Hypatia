#' Subset by cells
#'
#' Subsets a `SingleCellExperiment` object by cell metadata or cell IDs.
#'
#' @param object A `SingleCellExperiment` object.
#' @param subset Logical expression evaluated in `colData(object)`. Provide exactly one of `subset` or `cells`.
#' @param cells Optional vector of cell IDs. Provide exactly one of `subset` or `cells`.
#' @param invert Logical; if `TRUE`, remove the selected cells instead of keeping them.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns The object after cell selection, with transcript QC columns updated.
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

  subset_expr <- substitute(subset)
  has_subset <- !missing(subset) && !identical(subset_expr, quote(NULL))
  has_cells <- !is.null(cells)
  if (has_subset == has_cells) {
    stop("Please provide exactly one of the `subset` or `cells` parameters.")
  }

  ncell_before <- ncol(object)

  # Subset expression
  if (has_subset) {
    try_res <- try(
      {
        cdata <- as.data.frame(colData(object))
        names(cdata) <- names(colData(object))
        subset_cells <- cdata %>%
          filter({{ subset }}) %>%
          rownames()

        if (!invert) {
          object <- object[, colnames(object) %in% subset_cells, drop = FALSE]
        }
        if (invert) {
          object <- object[, !colnames(object) %in% subset_cells, drop = FALSE]
        }
      },
      silent = TRUE
    )

    if (inherits(try_res, "try-error")) {
      stop("The subset term was not valid. Check logical expression?", call. = FALSE)
    }
  }
  # Cell IDs were provided
  else if (has_cells) {
    if (!quiet) message("Subsetting using provided cell IDs.")

    try_res <- try(
      {
        if (!invert) {
          object <- object[, colnames(object) %in% cells, drop = FALSE]
        }
        if (invert) {
          object <- object[, !colnames(object) %in% cells, drop = FALSE]
        }
      },
      silent = TRUE
    )

    if (inherits(try_res, "try-error")) {
      stop("The provided cell IDs were not valid.", call. = FALSE)
    }
  }

  # Update meta cols
  counts <- assay(object, "counts")
  rowData(object)[["nCell"]] <- rowSums(counts > 0)

  if (!quiet) message("Cells before subset: ", ncell_before)
  if (!quiet) message("Cells after subset: ", ncol(object))

  return(object)
}
