#' Get isoform expression summaries
#'
#' Summarizes transcript detection and average expression by cell group.
#'
#' @param object A `SingleCellExperiment` object.
#' @param transcripts Vector of active transcript IDs to summarize. Ignored when `genes` is supplied.
#' @param genes Optional vector of active gene IDs. When supplied, all associated transcripts are summarized.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.subset Optional vector of group labels to include.
#' @param assay.use Assay name to use.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell expression across all cells in each group.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The (associated) gene being queried.}
#'   \item{`transcript`}{The transcript being queried.}
#'   \item{`pct`}{Fraction of cells in `group` with expression of the transcript, on the `[0, 1]` scale.}
#'   \item{`avgExpr`}{Average expression of the transcript across all cells in `group`.}
#'   \item{`cell.n`}{Number of cells used for cell-level summaries.}
#'   \item{`cell.expr.median`, `cell.expr.sd`, `cell.expr.iqr`}{Median, sample standard deviation, and interquartile range of cell-level expression.}
#' }
#' Dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#' @details Supply either transcript IDs or gene IDs; when `genes` is supplied,
#' all of its associated transcripts are summarized and `transcripts` is ignored.
#' Queries use the active identifiers. Multiple `group.by` columns are joined
#' with `_`, and `group.subset` selects groups without pooling them.
#'
#' The mean expression and optional cell-level summaries use the selected assay
#' across all cells in each group, including zeros. Detection (`pct`) is the
#' fraction of cells with an assay value above zero. With `assay.use = "logcounts"`,
#' `avgExpr` is the mean of log-normalized values, not the logarithm of a mean.
#' The requested assay must already exist; use [NormalizeCounts()] to create it.
#' This function summarizes expression without statistical testing or a minimum
#' detection filter.
#' @seealso [RunDEI()], [PlotExpression()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom S4Vectors metadata

GetExpression <- function (
    object,
    transcripts = NULL,
    genes = NULL,
    group.by = NULL,
    group.subset = NULL,
    assay.use = "logcounts",
    quiet = FALSE,
    cell.dispersion = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  if (is.null(genes)) {
    assertCharacter(transcripts, any.missing = FALSE, unique = TRUE, null.ok = TRUE)
  }
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(genes, any.missing = FALSE, unique = TRUE, null.ok = TRUE)
  assertCharacter(group.subset, null.ok = TRUE)
  assertTRUE(assay.use %in% assayNames(object))
  assertFlag(quiet)
  assertFlag(cell.dispersion)

  # Active transcript and gene names
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame(optional = TRUE) %>%
    rownames_to_column(var = "transcripts_query") %>%
    rename("gene_query" = all_of(active.gene.id))

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  unique_groups <- unique(colData(object)$group_var)
  ## check group subset
  assertSubset(group.subset, unique_groups, empty.ok = TRUE)
  ## subset cells
  if (!is.null(group.subset)) {
    object <- object[, colData(object)$group_var %in% group.subset, drop = FALSE]
  }

  # Transcript or gene filter
  if (is.null(genes)) {
    transcript_filter <- .FilterTranscripts(object, transcripts, quiet = quiet)
    object <- transcript_filter$object
    transcripts <- transcript_filter$transcripts
  } else {
    gene_filter <- .FilterGenes(object, genes, active.gene.id, quiet = quiet)
    object <- gene_filter$object
    genes <- gene_filter$genes
  }

  # Grouped expression summaries
  expr_mat <- assay(object, assay.use)
  col_group <- colData(object)[["group_var"]]
  n_cells_grp <- table(col_group)
  grp_tx_cts <- t(rowsum(t(expr_mat), group = col_group))
  grp_tx_avg <- sweep(grp_tx_cts, 2, n_cells_grp, FUN = "/")
  grp_tx_pos_cts <- t(rowsum(t(expr_mat > 0) * 1, group = col_group))
  grp_tx_pct <- sweep(grp_tx_pos_cts, 2, n_cells_grp, FUN = "/")

  grp_tx_avg <- grp_tx_avg %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    pivot_longer(-transcripts_query, names_to = "group", values_to = "avgExpr")
  grp_tx_pct <- grp_tx_pct %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    pivot_longer(-transcripts_query, names_to = "group", values_to = "pct")

  if (cell.dispersion) {
    if (!quiet) message("Calculating cell-to-cell dispersion...")
    dispersion_list <- lapply(unique(col_group), function(group) {
      group_dispersion <- .CellExpressionDispersion(
        expr_mat[, col_group == group, drop = FALSE]
      )
      group_dispersion$group <- group
      group_dispersion %>%
        rename(transcripts_query = transcript)
    })
    dispersion_data <- purrr::reduce(dispersion_list, rbind)
  } else {
    dispersion_data <- full_join(
      grp_tx_avg[c("group", "transcripts_query")],
      grp_tx_pct[c("group", "transcripts_query")],
      by = c("group", "transcripts_query")
    ) %>%
      mutate(
        cell.n = NA_integer_,
        cell.expr.median = NA_real_,
        cell.expr.sd = NA_real_,
        cell.expr.iqr = NA_real_
      )
  }

  result <- full_join(grp_tx_avg, grp_tx_pct, by = c("group", "transcripts_query")) %>%
    left_join(dispersion_data, by = c("group", "transcripts_query")) %>%
    left_join(gene.id.df, by = "transcripts_query") %>%
    select(group, gene_query, transcripts_query, pct, avgExpr, cell.n,
           cell.expr.median, cell.expr.sd, cell.expr.iqr) %>%
    arrange(group, transcripts_query) %>%
    rename(
      "gene" = "gene_query",
      "transcript" = "transcripts_query"
    ) %>%
    as.data.frame()

  return(result)
}
