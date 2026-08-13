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
#' @param only.pos Logical; if `TRUE`, only transcripts with positive fold change will be reported.
#' @param transcripts Optional vector of active transcript IDs to test.
#' @param p.adj P-value adjustment method. Must be one of `stats::p.adjust.methods`.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'   \item{`gene`}{The gene associated with the transcript being tested.}
#'   \item{`transcript`}{The transcript being tested.}
#'   \item{`pct.1`}{Percentage of cells in `group.1` with expression of the transcript.}
#'   \item{`pct.2`}{Percentage of cells in `group.2` with expression of the transcript.}
#'   \item{`avgExpr.1`}{Average expression of the transcript in `group.1`.}
#'   \item{`avgExpr.2`}{Average expression of the transcript in `group.2`.}
#'   \item{`log2FC`}{The log2 fold change in transcript expression between the two groups (`group.1` - `group.2`).}
#'   \item{`pval`}{P-value from the Wilcoxon rank-sum test.}
#'   \item{`padj`}{Adjusted p-value, calculated separately for each group comparison.}
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column column_to_rownames
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
  quiet = FALSE
  ) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertTRUE(identical(rownames(colData(object)), colnames(object)))
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
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

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame() %>%
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
      min.valid = 2,
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
      ) %>%
      select(group.1, group.2, gene, transcript, pct.1, pct.2, avgExpr.1, avgExpr.2, log2FC, pval, padj)

    data_list[[comp]] <- results
  }

  # combine results from across comparisons
  if (length(data_list) == 0) {
    stop("0 transcripts passed detection thresholds (check min.pct).", call. = FALSE)
  }
  final_results <- reduce(data_list, rbind)

  if (!quiet) message("Done.")

  return(final_results)
}
