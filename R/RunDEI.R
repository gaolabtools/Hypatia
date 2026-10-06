#' Run differential isoform expression analysis
#'
#' Tests transcripts for differential isoform expression between cell groups.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.1 Group label(s) for the first side of the comparison. If `NULL`, each group is compared against all others.
#' @param group.2 Optional group label(s) for the second side of the comparison. If `NULL`, `group.1` is compared against all other cells.
#' @param assay.use Assay name to use. The default and recommended assay is `"logcounts"`.
#' @param min.pct Minimum fraction of cells in each group where the transcript must be detected.
#' @param only.pos Logical; if `TRUE`, only transcripts with non-negative log2 fold change will be reported.
#' @param transcripts Optional vector of active transcript IDs to test.
#' @param p.adj P-value adjustment method. Must be one of `stats::p.adjust.methods`.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell expression across all cells in each comparison group.
#'
#' @returns A list containing two data frames:
#'
#' \describe{
#'   \item{`$data`}{A data frame of expression summaries with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene associated with the transcript being tested.}
#'       \item{`transcript`}{The transcript being tested.}
#'       \item{`pct.1`, `pct.2`}{Fraction of cells in each group with expression of the transcript, on the `[0, 1]` scale.}
#'       \item{`avgExpr.1`, `avgExpr.2`}{Mean expression of the transcript across all cells in each group.}
#'       \item{`cell.n.1`, `cell.n.2`}{Number of cells used for cell-level summaries in each group.}
#'       \item{`cell.expr.median.1`, `cell.expr.median.2`, `cell.expr.sd.1`, `cell.expr.sd.2`, `cell.expr.iqr.1`, `cell.expr.iqr.2`}{Median, sample standard deviation, and interquartile range of cell-level expression.}
#'     }
#'     Dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#'   }
#'   \item{`$stats`}{A data frame of statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene associated with the transcript being tested.}
#'       \item{`transcript`}{The transcript being tested.}
#'       \item{`log2FC`}{`log2(avgExpr.1 / avgExpr.2)` on the selected assay, without a pseudocount. With `"logcounts"`, this is the log2 ratio of mean log-normalized values.}
#'       \item{`pval`}{P-value from the Wilcoxon rank-sum test.}
#'       \item{`padj`}{Adjusted p-value, calculated separately for each group comparison.}
#'     }
#'   }
#' }
#' Rows in `$data` are ordered by `group.1`, `group.2`, `gene`, and
#' `transcript`. Rows in `$stats` are ordered by `group.1`, `group.2`,
#' `padj`, `gene`, and `transcript`.
#' @details The selected assay must already exist; the default is `"logcounts"`
#' from [NormalizeCounts()]. Transcripts must be detected in at least `min.pct`
#' of cells in both comparison groups before testing. Detection means a selected
#' assay value above zero. A two-sided Wilcoxon rank-sum test compares the
#' cell-level assay values, including zeros, using `matrixTests` with automatic
#' exact/asymptotic p-value selection.
#'
#' `avgExpr` is the arithmetic mean on the selected assay, and `log2FC` is
#' `log2(avgExpr.1 / avgExpr.2)` without a pseudocount. With `"logcounts"`, this
#' compares mean log-normalized values. Optional cell-level dispersion summaries
#' include every cell in each group. `only.pos = TRUE` retains non-negative
#' `log2FC` values after p-value adjustment.
#'
#' Multiple `group.by` columns are joined with `_`. Omitting both comparison
#' arguments compares each group with the remaining cells; exactly two groups
#' produce one comparison. Multiple labels on either side pool their cells.
#' P-values are adjusted across tested transcripts separately within each
#' comparison using Bonferroni by default. An error is returned if no transcripts
#' pass detection filtering across all comparisons.
#' @seealso [GetExpression()], [PlotExpression()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom purrr reduce
#' @importFrom matrixTests row_wilcoxon_twosample
#' @importFrom stats p.adjust

RunDEI <- function(
  object,
  group.by = NULL,
  group.1 = NULL,
  group.2 = NULL,
  assay.use = "logcounts",
  min.pct = 0.01,
  only.pos = FALSE,
  transcripts = NULL,
  p.adj = "bonferroni",
  quiet = FALSE,
  cell.dispersion = FALSE
  ) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertTRUE(identical(rownames(colData(object)), colnames(object)))
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.1, null.ok = TRUE)
  assertCharacter(group.2, null.ok = TRUE)
  assertTRUE(length(intersect(group.1, group.2)) == 0)
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(min.pct, lower = 0, upper = 1, finite = TRUE)
  assertFlag(only.pos)
  assertCharacter(transcripts, null.ok = TRUE, any.missing = FALSE, unique = TRUE)
  p.adj <- .PAdjustMethod(p.adj)
  assertFlag(quiet)
  assertFlag(cell.dispersion)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame(optional = TRUE) %>%
    rownames_to_column(var = "transcript") %>%
    rename("gene" = all_of(active.gene.id))

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  cell_groups <- colData(object)$group_var
  unique_groups <- unique(cell_groups)

  if (length(unique_groups) < 2) {
    stop("There must be at least 2 groups to compare.")
  }

  ## check group.1 and group.2
  assertSubset(group.1, choices = unique_groups, empty.ok = TRUE)
  assertSubset(group.2, choices = setdiff(unique_groups, group.1), empty.ok = TRUE)


  # Transcript filter
  if (!is.null(transcripts)) {

    transcript_filter <- .FilterTranscripts(
      object,
      transcripts,
      quiet = quiet,
      min.valid = 1,
      none.message = "None of the transcripts provided were found in the object."
    )
    object <- transcript_filter$object
    transcripts <- transcript_filter$transcripts
  }


  # DEI
  object_grp_list <- .BuildGroupComparisons(object, group.1, group.2, unique_groups)
  comparison_mode <- attr(object_grp_list, "mode")
  if (!quiet && comparison_mode == "all") {
    message("Running DEI analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")
  } else if (!quiet && comparison_mode == "one_vs_all") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DEI analysis for ", comparison$grp1.names, " vs all other cells...")
  } else if (!quiet && comparison_mode == "pair") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DEI analysis for ", comparison$grp1.names, " vs ", comparison$grp2.names, "...")
  }

  # Loop through object grp list
  data_list <- list()
  stats_list <- list()

  for (comp in names(object_grp_list)) {

    if (!quiet && comp != "single_test" && length(unique_groups) > 2) message("  ", comp, "... ")

    ## get group objects and names
    object_grp1 <- object_grp_list[[comp]]$grp1.object
    object_grp2 <- object_grp_list[[comp]]$grp2.object
    group.1 <- object_grp_list[[comp]]$grp1.names
    group.2 <- object_grp_list[[comp]]$grp2.names

    ## count mat for each group
    expr_mat_grp1 <- assay(object_grp1, assay.use)
    expr_mat_grp2 <- assay(object_grp2, assay.use)

    ## gene name
    gene_groups_grp1 <- rowData(object_grp1)[[active.gene.id]]

    ## mean expression for each group
    avg_grp1 <- rowMeans(expr_mat_grp1)
    avg_grp2 <- rowMeans(expr_mat_grp2)

    ## expression detection rates (pct) for each group
    pct_grp1 <- rowSums(expr_mat_grp1 > 0) / ncol(expr_mat_grp1)
    pct_grp2 <- rowSums(expr_mat_grp2 > 0) / ncol(expr_mat_grp2)

    ## report expression stats
    stopifnot(all(rownames(gene_groups_grp1) == rownames(avg_grp1)))
    expr_df <- data.frame("gene" = gene_groups_grp1,
                          "pct.grp1" = pct_grp1,
                          "pct.grp2" = pct_grp2,
                          "avg.grp1" = avg_grp1,
                          "avg.grp2" = avg_grp2
                          ) %>%
      # fold change
      mutate(log2FC = log2(avg.grp1 / avg.grp2))

    ## filter transcripts by min.pct before, dense, tests, and correction
    test_transcripts <- expr_df %>%
      filter(pct.grp1 >= min.pct & pct.grp2 >= min.pct) %>%
      rownames()
    if (length(test_transcripts) == 0) {
      next
    }
    expr_df <- expr_df[test_transcripts, , drop = FALSE]
    expr_mat_grp1 <- expr_mat_grp1[test_transcripts, , drop = FALSE]
    expr_mat_grp2 <- expr_mat_grp2[test_transcripts, , drop = FALSE]
    stopifnot(all(rownames(expr_mat_grp1) == rownames(expr_mat_grp2)))
    stopifnot(all(rownames(expr_mat_grp1) == rownames(expr_df)))

    # Wilcox test using matrixTests
    ## numeric matrix is required
    mat_grp1_dense <- suppressWarnings(as.matrix(expr_mat_grp1))
    mat_grp2_dense <- suppressWarnings(as.matrix(expr_mat_grp2))

    ## cell-to-cell expression dispersion for tested transcripts
    if (cell.dispersion) {
      if (!quiet) message("  Calculating cell-to-cell dispersion...")
      dispersion_grp1 <- .CellExpressionDispersion(mat_grp1_dense)
      dispersion_grp2 <- .CellExpressionDispersion(mat_grp2_dense)
    } else {
      dispersion_grp1 <- data.frame(
        transcript = test_transcripts,
        cell.n = rep(NA_integer_, length(test_transcripts)),
        cell.expr.median = rep(NA_real_, length(test_transcripts)),
        cell.expr.sd = rep(NA_real_, length(test_transcripts)),
        cell.expr.iqr = rep(NA_real_, length(test_transcripts))
      )
      dispersion_grp2 <- dispersion_grp1
    }
    stopifnot(identical(dispersion_grp1$transcript, rownames(expr_df)))
    stopifnot(identical(dispersion_grp2$transcript, rownames(expr_df)))
    names(dispersion_grp1)[-1] <- paste0(names(dispersion_grp1)[-1], ".1")
    names(dispersion_grp2)[-1] <- paste0(names(dispersion_grp2)[-1], ".2")
    expr_df <- cbind(
      expr_df,
      dispersion_grp1[, -1, drop = FALSE],
      dispersion_grp2[, -1, drop = FALSE]
    )

    ## test
    test_result <- row_wilcoxon_twosample(
      x = mat_grp1_dense,
      y = mat_grp2_dense,
      alternative = "two.sided",
      exact = NA
      )
    test_result <- test_result[, "pvalue", drop = FALSE]
    stopifnot(all(rownames(test_result) == rownames(expr_df)))

    ## combine test results with expression stats
    results <- cbind(expr_df, test_result)
    ## add group names
    results <- results %>%
      mutate("group.1" = group.1,
             "group.2" = group.2,
             .before = "gene")

    ## calculate adjusted p-values per comparison
    results <- results %>%
      mutate(padj = p.adjust(pvalue, method = p.adj)) %>%
      arrange(padj) %>%
      ungroup()

    ## filter only positive log2FC
    if (only.pos) {
      results <- results %>%
        filter(log2FC >= 0)
    }

    ## output format
    results <- results %>%
      rownames_to_column(var = "transcript") %>%
      rename(
        "pct.1" = "pct.grp1",
        "pct.2" = "pct.grp2",
        "avgExpr.1" = "avg.grp1",
        "avgExpr.2" = "avg.grp2",
        "pval" = "pvalue"
      )

    data_list[[comp]] <- results %>%
      select(group.1, group.2, gene, transcript, pct.1, pct.2,
             avgExpr.1, avgExpr.2, cell.n.1, cell.expr.median.1,
             cell.expr.sd.1, cell.expr.iqr.1, cell.n.2,
             cell.expr.median.2, cell.expr.sd.2, cell.expr.iqr.2)
    stats_list[[comp]] <- results %>%
      select(group.1, group.2, gene, transcript, log2FC, pval, padj)
  }

  # combine results from across comparisons
  if (length(data_list) == 0) {
    stop("0 transcripts passed filtering (check min.pct).", call. = FALSE)
  }
  final_data <- reduce(data_list, rbind) %>%
    arrange(group.1, group.2, gene, transcript)
  final_stats <- reduce(stats_list, rbind) %>%
    arrange(group.1, group.2, padj, gene, transcript)

  if (!quiet) message("Done.")

  return(list(data = final_data, stats = final_stats))
}
