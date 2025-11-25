#' Get summarized isoform expression
#'
#' Retrieves summarized isoform expression data of one or more transcripts.
#'
#' @param object A `SingleCellExperiment` object.
#' @param transcripts A vector containing one or more transcript IDs.
#' @param genes A vector containing one or more gene IDs. Will retrieve data for associated transcripts.
#' @param group.by Name of `colData` variable to group cells. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.subset An optional vector specifying a subset of the elements in `group.by` to include.
#' @param assay.use Which `assay` (counts) to use.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The (associated) gene being queried.}
#'   \item{`transcript`}{The transcript being queried.}
#'   \item{`transcript.pct`}{Percentage of cells in `group` with expression of the transcript.}
#'   \item{`avgExpr`}{Average expression of the transcript across all cells in `group`.}
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr pivot_longer unite
#' @importFrom tibble rownames_to_column
#' @importFrom S4Vectors metadata

GetExpression <- function (
    object,
    transcripts,
    genes = NULL,
    group.by = NULL,
    group.subset = NULL,
    assay.use = "logcounts",
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  if (is.null(genes)) {
    assertCharacter(transcripts, any.missing = FALSE, unique = TRUE, null.ok = TRUE)
  }
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertTRUE(assay.use %in% assayNames(object))
  assertFlag(quiet)

  # Active transcript and gene names
  assertString(metadata(object)$active.transcript.id)
  assertString(metadata(object)$active.gene.id)
  active.transcript.id <- metadata(object)$active.transcript.id
  active.gene.id <- metadata(object)$active.gene.id
  assertChoice(active.gene.id, colnames(rowData(object)))
  assertFALSE(anyMissing(rowData(object)[[active.gene.id]]))

  if (metadata(object)$active.transcript.id != "") {
    assertChoice(active.transcript.id, colnames(rowData(object)))
    assertFALSE(any(duplicated(rowData(object)[[active.transcript.id]])))
    assertFALSE(anyMissing(rowData(object)[[active.transcript.id]]))
    rownames(object) <- rowData(object)[[active.transcript.id]]
  }

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    rename("gene_query" = all_of(active.gene.id))

  # Group structure
  colData(object)$group_var <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE) %>%
    pull(group_var)
  unique_groups <- unique(colData(object)$group_var)
  ## check group subset
  assertSubset(group.subset, unique_groups, empty.ok = TRUE)
  ## subset cells
  if (!is.null(group.subset)) {
    object <- object[, colData(object)$group_var %in% group.subset]
  }

  # Transcripts provided
  if (is.null(genes)) {

    ## transcript filter
    if (any(transcripts %in% rownames(object))) {
      missing_transcripts <- setdiff(transcripts, rownames(object))
      if (length(missing_transcripts) > 0) {
        if (!quiet) message("\u2139 Warning: The following transcripts were not found in the object: '", paste0(missing_transcripts, collapse = "', '"), "'.")
        transcripts <- transcripts[transcripts %in% rownames(object)]
      }
    }
    else {
      stop("None of the transcripts provided were found in the object. (Check active.transcript.id?)")
    }

    ## expression mat
    object <- object[rownames(object) %in% transcripts, , drop = FALSE]
    expr_mat <- assay(object, assay.use)
    ## calculate tx counts, avg expr, and pct per group
    col_group <- colData(object)[["group_var"]]
    n_cells_grp <- table(col_group)
    grp_tx_cts <- t(rowsum(t(expr_mat), group = col_group))
    grp_tx_avg <- sweep(grp_tx_cts, 2, n_cells_grp, FUN = "/")
    grp_tx_pos_cts <- t(rowsum(t(expr_mat > 0) * 1, group = col_group))
    grp_tx_pct <- sweep(grp_tx_pos_cts, 2, n_cells_grp, FUN = "/")

    ## output
    grp_tx_avg <- grp_tx_avg %>%
      as.data.frame() %>%
      rownames_to_column(var = "transcripts_query") %>%
      pivot_longer(-transcripts_query, names_to = "group", values_to = "avgExpr")
    grp_tx_pct <- grp_tx_pct %>%
      as.data.frame() %>%
      rownames_to_column(var = "transcripts_query") %>%
      pivot_longer(-transcripts_query, names_to = "group", values_to = "pct")

    result <- full_join(grp_tx_avg, grp_tx_pct, by = c("group", "transcripts_query"))

    result <- result %>%
      left_join(., gene.id.df, by = "transcripts_query", relationship = "many-to-many") %>%
      select(group, gene_query, transcripts_query, pct, avgExpr) %>%
      arrange(group, transcripts_query) %>%
      rename("gene" = "gene_query",
                    "transcript" = "transcripts_query") %>%
      as.data.frame()

    return(result)

  }

  # Genes provided
  else {

    ## active gene names
    assertString(metadata(object)$active.gene.id)
    active.gene.id <- metadata(object)$active.gene.id
    assertChoice(active.gene.id, colnames(rowData(object)))
    assertFALSE(anyMissing(rowData(object)[[active.gene.id]]))

    ## gene filter
    if (any(genes %in% unique(rowData(object)[[active.gene.id]]))) {
      missing_genes <- setdiff(genes, unique(rowData(object)[[active.gene.id]]))
      if (length(missing_genes) == length(genes)) {
        stop("None of the genes were found in the object. (Check active.gene.id?)")
      }
      if (length(missing_genes) > 0) {
        if (!quiet) message("\u2139 Warning: The following genes were not found in the object: '", paste0(missing_genes, collapse = "', '"), "'.")
        genes <- genes[genes %in% unique(rowData(object)[[active.gene.id]])]
      }

      ## expression mat
      object <- object[rowData(object)[[active.gene.id]] %in% genes, , drop = FALSE]
      expr_mat <- assay(object, assay.use)
      ## calculate tx counts, avg expr, and pct per group
      col_group <- colData(object)[["group_var"]]
      n_cells_grp <- table(col_group)
      grp_tx_cts <- t(rowsum(t(expr_mat), group = col_group))
      grp_tx_avg <- sweep(grp_tx_cts, 2, n_cells_grp, FUN = "/")
      grp_tx_pos_cts <- t(rowsum(t(expr_mat > 0) * 1, group = col_group))
      grp_tx_pct <- sweep(grp_tx_pos_cts, 2, n_cells_grp, FUN = "/")

      ## output
      grp_tx_avg <- grp_tx_avg %>%
        as.data.frame() %>%
        rownames_to_column(var = "transcripts_query") %>%
        pivot_longer(-transcripts_query, names_to = "group", values_to = "avgExpr")
      grp_tx_pct <- grp_tx_pct %>%
        as.data.frame() %>%
        rownames_to_column(var = "transcripts_query") %>%
        pivot_longer(-transcripts_query, names_to = "group", values_to = "pct")

      result <- full_join(grp_tx_avg, grp_tx_pct, by = c("group", "transcripts_query"))

      result <- result %>%
        left_join(., gene.id.df, by = "transcripts_query", relationship = "many-to-many") %>%
        select(group, gene_query, transcripts_query, pct, avgExpr) %>%
        arrange(group, transcripts_query) %>%
        rename(
          "gene" = "gene_query",
          "transcript" = "transcripts_query",
          "transcript.pct" = "pct"
          ) %>%
        as.data.frame()

      return(result)

    } else {
      stop("None of the genes were found in the object. (Check active.gene.id?)")
    }

  }

}
