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
#' @importFrom presto wilcoxauc

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

  ## every group vs all
  if (is.null(group.1) && is.null(group.2)) {

    if (!quiet) message("Running DEI analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")
    results <- wilcoxauc(X = expr_mat, y = cell_groups)
    group.2.updated <- NULL

    if (length(unique(cell_groups)) == 2) {
      results <- results %>%
        filter(group == unique(cell_groups)[1])
    }

  }

  ## 1 group vs all
  else if (!is.null(group.1) && is.null(group.2)) {

    # assign group.2
    group.2.updated <- setdiff(unique_groups, group.1)
    # update cell groups
    cell_groups_updated <- ifelse(cell_groups %in% group.1, paste0(group.1, collapse = ","), cell_groups)
    cell_groups_updated <- ifelse(cell_groups_updated %in% group.2.updated, paste0(group.2.updated, collapse = ","), cell_groups_updated)
    group.1.updated <- paste0(group.1, collapse = ",")
    group.2.updated <- paste0(group.2.updated, collapse = ",")
    # wilcox test
    if (!quiet) message("Running DEI analysis for ", group.1.updated, " vs all other cells...")
    results <- wilcoxauc(X = expr_mat, y = cell_groups_updated, groups_use = c(group.1.updated, group.2.updated))
    results <- results %>% filter(group == group.1.updated) # filter result for specified group

  }

  ## 2 groups comparison
  else if (!is.null(group.1) && !is.null(group.2)) {

    # update cell groups
    cell_groups_updated <- ifelse(cell_groups %in% group.1, paste0(group.1, collapse = ","), cell_groups)
    cell_groups_updated <- ifelse(cell_groups_updated %in% group.2, paste0(group.2, collapse = ","), cell_groups_updated)
    group.1.updated <- paste0(group.1, collapse = ",")
    group.2.updated <- paste0(group.2, collapse = ",")
    # wilcox test
    if (!quiet) message("Running DEI analysis for ", group.1.updated, " vs ", group.2.updated, "...")
    results <- wilcoxauc(X = expr_mat, y = cell_groups_updated, groups_use = c(group.1.updated, group.2.updated))
    results <- results %>% filter(group == group.1.updated) # presto wilcoxauc outputs both directions, use group.1 direction

  }

  ## case: group 1 unspecified but group 2 is specified
  else if (is.null(group.1) && !is.null(group.2)) {
    stop("`group.1` must be specified prior to `group.2`")
  }


  # Results
  ## filter results by min.pct before pval correction
  results <- results %>%
    mutate(pct_in = pct_in / 100,
           pct_out = pct_out / 100) %>%
    filter(pct_in >= min.pct & pct_out >= min.pct)
  ## calculate valid log2FC (presto's logFC column is actually diff between means)
  results <- results %>%
    rename("avgExpr.1" = "avgExpr") %>%
    mutate(avgExpr.2 = avgExpr.1 - logFC,
           log2FC = log2(avgExpr.1 / avgExpr.2))
  ## calculate bonferroni corrected p-values (per group)
  results <- results %>%
    select(-c(padj)) %>%
    group_by(group) %>%
    mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
    arrange(desc(padj)) %>%
    ungroup()
  ## filter only positive log2FC
  if (only.pos) {
    results <- results %>%
      filter(log2FC >= 0)
  }

  ## output
  results <- results %>%
    rename(
      "transcript" = "feature",
      "pct.1" = "pct_in",
      "pct.2" = "pct_out",
      "group.1" = "group"
      ) %>%
    mutate("group.2" = ifelse(is.null(group.2.updated), NA, group.2.updated)) %>%
    left_join(., gene.id.df, by = "transcript") %>%
    select(group.1, group.2, gene, transcript, pct.1, pct.2, avgExpr.1, avgExpr.2, log2FC, pval, padj) %>%
    arrange(group.1, padj) %>%
    as.data.frame()

  if (is.null(group.2.updated)) {
    results <- results %>%
      rowwise() %>%
      mutate(group.2 = paste0(setdiff(unique_groups, group.1), collapse = ",")) %>%
      ungroup() %>%
      as.data.frame()

  }

  if (!quiet) message("Done.")

  return(results)
}
