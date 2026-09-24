#' Get isoform usage summaries
#'
#' Summarizes isoform counts and proportions for one or more genes.
#'
#' @param object A `SingleCellExperiment` object.
#' @param genes Vector of active gene IDs to summarize.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.subset Optional vector of group labels to include.
#' @param assay.use Assay name to use.
#' @param min.tx.cts Minimum transcript counts required before proportions are calculated.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell transcript proportions among cells with positive retained-gene counts.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The gene being queried.}
#'   \item{`gene.pct`}{Fraction of cells in `group` with expression of the gene, on the `[0, 1]` scale.}
#'   \item{`transcript`}{The associated transcript.}
#'   \item{`cts`}{Total counts of the transcript across cells in `group`.}
#'   \item{`prop`}{Transcript proportion.}
#'   \item{`cell.n`}{Number of gene-positive cells used for cell-level summaries.}
#'   \item{`cell.prop.mean`, `cell.prop.median`, `cell.prop.sd`, `cell.prop.iqr`}{Mean, median, sample standard deviation, and interquartile range of cell-level transcript proportions.}
#'   \item{`cell.prop.zero.frac`}{Fraction of gene-positive cells in which the transcript is not detected.}
#' }
#' Dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#' @details Counts are pooled across cells separately for each gene and group.
#' Within each group, transcripts below `min.tx.cts` are removed before dividing
#' transcript counts by the retained gene total. Consequently, `prop` describes
#' pooled usage, rather than the average of cell-level proportions. `gene.pct`
#' is measured before transcript count filtering.
#'
#' With `cell.dispersion = TRUE`, proportions are also calculated separately in
#' cells with positive retained-gene counts. Undetected transcripts contribute
#' zero within these eligible cells; cells with no retained-gene counts are
#' excluded. These summaries describe biological cell-to-cell variation.
#'
#' Gene queries use the active gene IDs, and transcript labels use the active
#' transcript IDs or object row names. Multiple `group.by` columns are joined
#' with `_`. `group.subset` selects groups without pooling them.
#' @seealso [RunDIU()], [PlotUsage()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr pivot_longer

GetUsage <- function (
    object,
    genes,
    group.by = NULL,
    group.subset = NULL,
    assay.use = "counts",
    min.tx.cts = 1,
    quiet = FALSE,
    cell.dispersion = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertCharacter(genes, any.missing = FALSE, unique = TRUE)
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.subset, null.ok = TRUE)
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertFlag(quiet)
  assertFlag(cell.dispersion)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    rename("gene_query" = all_of(active.gene.id))

  # Gene filter
  gene_filter <- .FilterGenes(object, genes, active.gene.id, quiet = quiet)
  object <- gene_filter$object
  genes <- gene_filter$genes

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  unique_groups <- unique(colData(object)$group_var)
  ## check group subset
  assertSubset(group.subset, unique_groups, empty.ok = TRUE)
  ## subset cells
  if (!is.null(group.subset)) {
    object <- object[, colData(object)$group_var %in% group.subset, drop = FALSE]
  }

  # Expression mat
  expr_mat <- assay(object[rowData(object)[[active.gene.id]] %in% genes, , drop = FALSE], assay.use)
  col_group <- colData(object)[["group_var"]]
  ## calculate tx counts per group
  grp_tx_cts <- t(rowsum(t(expr_mat), col_group))
  ## calculate gene pct per group
  n_cells_grp <- table(col_group)
  row_group <- rowData(object)[[active.gene.id]]
  gene_cts <- rowsum(expr_mat, group = row_group)
  grp_gene_pos_cts <- t(rowsum(t(gene_cts > 0) * 1, group = col_group))
  grp_gene_pct <- sweep(grp_gene_pos_cts, 2, n_cells_grp, FUN = "/")
  grp_gene_pct <- grp_gene_pct %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_query") %>%
    pivot_longer(-gene_query, values_to = "gene.pct", names_to = "group_var")

  ## output
  grp_tx_cts <- grp_tx_cts %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    ## add gene_ids
    left_join(., gene.id.df, by = "transcripts_query") %>%
    pivot_longer(-c("transcripts_query", "gene_query"), names_to = "group_var", values_to = "grp_cts") %>%
    ## filter isoforms by counts
    filter(grp_cts >= min.tx.cts) %>%
    ## calculate isoform props
    group_by(gene_query, group_var) %>%
    mutate(prop = grp_cts / sum(grp_cts)) %>%
    ungroup()

  ## cell-to-cell transcript proportion dispersion
  if (cell.dispersion && nrow(grp_tx_cts) > 0) {
    if (!quiet) message("Calculating cell-to-cell dispersion...")
    dispersion_list <- lapply(unique(grp_tx_cts$group_var), function(group) {
      group_data <- grp_tx_cts %>%
        filter(group_var == group)
      group_object <- object[
        group_data$transcripts_query,
        colData(object)$group_var == group,
        drop = FALSE
      ]
      group_dispersion <- .CellUsageDispersion(
        assay(group_object, assay.use),
        rowData(group_object)[[active.gene.id]]
      )
      group_dispersion$group_var <- group
      group_dispersion %>%
        rename(transcripts_query = transcript)
    })
    dispersion_data <- purrr::reduce(dispersion_list, rbind)
  } else {
    dispersion_data <- grp_tx_cts %>%
      distinct(group_var, transcripts_query) %>%
      mutate(
        cell.n = NA_integer_,
        cell.prop.mean = NA_real_,
        cell.prop.median = NA_real_,
        cell.prop.sd = NA_real_,
        cell.prop.iqr = NA_real_,
        cell.prop.zero.frac = NA_real_
      )
  }

  result <- grp_tx_cts %>%
    left_join(., grp_gene_pct, by = c("group_var", "gene_query")) %>%
    left_join(., dispersion_data, by = c("group_var", "transcripts_query")) %>%
    select(group_var, gene_query, gene.pct, transcripts_query, grp_cts, prop,
           cell.n, cell.prop.mean, cell.prop.median, cell.prop.sd,
           cell.prop.iqr, cell.prop.zero.frac) %>%
    rename("transcript" = "transcripts_query",
           "gene" = "gene_query",
           "group" = "group_var",
           "cts" = "grp_cts",
           "prop" = "prop") %>%
    arrange(group, gene) %>%
    as.data.frame()

  return(result)
}
