#' Run differential isoform expression analysis
#'
#' Runs differential isoform expression analysis, which compares the mean isoform expression across two groups of cells using Wilcoxon-based testing.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by Name of `colData` variable to group cells for comparisons. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.1 First group in the comparison.
#' @param group.2 Second group in the comparison.
#' @param assay.use Which `assay` (counts) to use. The default and recommended option is `"logcounts"`.
#' @param min.pct Minimum percentage of cells in which a transcript must be expressed in both groups for it to be adjusted for multiple testing and reported in results.
#' @param only.pos Logical; if `TRUE`, only transcripts with positive fold change will be reported.
#' @param transcripts A vector of transcript IDs to test.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group.1` & `group.1`}{The two cell groups being compared.}
#'   \item{`gene`}{The gene associated with the transcript being tested.}
#'   \item{`transcript`}{The transcript being tested.}
#'   \item{`pct.1`}{Average expression of the transcript across all cells in `group.1`.}
#'   \item{`pct.2`}{Average expression of the transcript across all cells in `group.2`.}
#'   \item{`log2FC`}{The log2 fold change in transcript expression between the two groups (`group.1` - `group.2`).}
#'   \item{`pval`}{P-value from the the Wilcoxon rank-sum test.}
#'   \item{`padj`}{Adjusted p-value, calculated separately for each group comparison (default: Bonferroni).}
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom matrixTests row_wilcoxon_twosample

RunDEI <- function(
  object,
  group.by = NULL,
  group.1 = NULL,
  group.2 = NULL,
  assay.use = "logcounts",
  min.pct = 0.01,
  only.pos = FALSE,
  transcripts = NULL,
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
  assertFlag(quiet)

  # Transcript and gene IDs
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
    rownames_to_column(var = "transcript") %>%
    rename("gene" = all_of(active.gene.id))

  # Expression mat
  expr_mat <- assay(object, assay.use)

  # Group structure
  colData(object)$group_var <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE) %>%
    pull(group_var)
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

    # transcripts supplied are detected
    if (any(transcripts %in% rownames(object))) {
      missing_transcripts <- setdiff(transcripts, rownames(object))
      # if missing transcripts are detected
      if (length(missing_transcripts) > 0) {
        if (!quiet) message("\u2139 Warning: The following transcripts were not found in the object: '", paste0(missing_transcripts, collapse = "', '"), "'.")
        transcripts <- transcripts[transcripts %in% rownames(object)]
        if (length(transcripts) < 2) {
          stop("Please provide at least 2 valid transcripts to test.")
        } else {
          expr_mat <- expr_mat[transcripts, , drop = FALSE]
        }
      # all transcripts provided are in detected
      } else {
        if (length(transcripts) < 2) {
          stop("Please provide at least 2 valid transcripts to test.")
        } else {
          expr_mat <- expr_mat[transcripts, , drop = FALSE]
        }
      }
    # no transcripts supplied are detected
    } else {
      stop("None of the transcripts provided were found in the object.")
    }
  }


  # DEI
  object_grp_list <- list()

  ## every group vs all
  if (is.null(group.1) && is.null(group.2)) {

    # loop through each group to make object list
    for (grp in unique_groups) {

      # assign group.2
      group.1 <- grp
      group.2 <- setdiff(unique_groups, group.1)

      # create an object for grp1 and grp2
      object_grp1 <- object[, object$group_var == grp]
      object_grp2 <- object[, object$group_var != grp]

      # update cell groups
      group.1.updated <- paste0(group.1, collapse = ",")
      group.2.updated <- paste0(group.2, collapse = ",")

      # add to list
      object_grp_list[[grp]] <- list("grp1.object" = object_grp1,
                                     "grp2.object" = object_grp2,
                                     "grp1.names" = group.1.updated,
                                     "grp2.names" = group.2.updated)

      # if only two groups, only run one
      if (length(unique_groups) == 2) {
        break
      }

    }

    if (!quiet) message("Running DEI analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")

  }

  ## 1 group vs all
  else if (!is.null(group.1) && is.null(group.2)) {

    # assign group.2
    group.2 <- setdiff(unique_groups, group.1)

    # create an object for grp1 and grp2
    object_grp1 <- object[, object$group_var %in% group.1]
    object_grp2 <- object[, object$group_var %in% group.2]

    # update cell groups
    group.1.updated <- paste0(group.1, collapse = ",")
    group.2.updated <- paste0(group.2, collapse = ",")

    # add to list
    object_grp_list[["single_test"]] <- list("grp1.object" = object_grp1,
                                             "grp2.object" = object_grp2,
                                             "grp1.names" = group.1.updated,
                                             "grp2.names" = group.2.updated)

    if (!quiet) message("Running DEI analysis for ", group.1.updated, " vs all other cells...")
  }

  ## 2 groups comparison
  else if (!is.null(group.1) && !is.null(group.2)) {

    # create an object for grp1 and grp2
    object_grp1 <- object[, object$group_var %in% group.1]
    object_grp2 <- object[, object$group_var %in% group.2]

    # update cell groups
    group.1.updated <- paste0(group.1, collapse = ",")
    group.2.updated <- paste0(group.2, collapse = ",")

    # add to list
    object_grp_list[["single_test"]] <- list("grp1.object" = object_grp1,
                                             "grp2.object" = object_grp2,
                                             "grp1.names" = group.1.updated,
                                             "grp2.names" = group.2.updated)

    if (!quiet) message("Running DEI analysis for ", group.1.updated, " vs ", group.2.updated, "...")

  }

  ## case: group 1 unspecified but group 2 is specified
  else if (is.null(group.1) && !is.null(group.2)) {
    stop("`group.1` must be specified prior to `group.2`")
  }


  # Loop through object grp list
  data_list <- list()

  for (comp in names(object_grp_list)) {

    if (!quiet && comp != "single_test" && length(unique_groups) > 2) message(comp, "... ")

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

    ## calculate bonferroni corrected p-values (per group)
    results <- results %>%
      mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
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
  final_results <- reduce(data_list, rbind)

  if (!quiet) message("Done.")

  return(final_results)
}
