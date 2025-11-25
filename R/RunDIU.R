#' Run differential isoform usage analysis
#'
#' Runs differential isoform usage analysis.
#' For each comparison, genes are first filtered according to coverage, counts are aggregated across groups, then the Chi-square or Fisher's exact test is applied.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by Name of `colData` variable to group cells for comparisons. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.1 First group in the comparison.
#' @param group.2 Second group in the comparison.
#' @param assay.use Which `assay` (counts) to use.
#' @param method.use Statistical test, either "Chisq" for Chi-squared test or "Fisher" for Fisher's exact test.
#' @param min.gene.pct Minimum percentage of cells in which a gene must be expressed in both groups for it to be tested.
#' @param min.gene.cts Minimum total transcript counts in which a gene must have in both groups for it be tested.
#' @param min.tx.cts Minimum transcript counts required for a transcript to be included in contingency tables.
#' @param genes Vector of genes to test. Note: genes will still be subject to filtering.
#' @param only.valid Logical; if `TRUE`, only tests with valid Chi-square approximations will be reported and corrected for multiple testing.
#' @param simulate.p Logical; if `TRUE`, p-values are computed using a Monte Carlo simulation. Note: Fisher's exact test always uses simulation.
#' @param p.adj Method for p-value adjustment. Options are `"BH"` (Benjamini-Hochberg) or `"Bonferroni"`.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A list containing two data frames:
#'
#' \describe{
#'   \item{`$data`}{
#'     A data frame of summarized data with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`gene.pct.1`}{Percentage of cells in `group.1` with expression of the gene.}
#'       \item{`gene.pct.2`}{Percentage of cells in `group.2` with expression of the gene.}
#'       \item{`transcript`}{The associated transcript.}
#'       \item{`cts.1`}{Total counts of the transcript across all cells in `group.1`.}
#'       \item{`cts.2`}{Total counts of the transcript across all cells in `group.2`.}
#'       \item{`prop.1`}{Transcript proportion for `group.1`.}
#'       \item{`prop.2`}{Transcript proportion for `group.2`.}
#'       \item{`delta`}{The difference in transcript proportions between groups (`group.1` - `group.2`).}
#'     }
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`max.delta`}{The largest absolute difference in transcript proportions between `group.1` and `group.2`.}
#'       \item{`transcript`}{The transcript associated with `max.delta`.}
#'       \item{`pval`}{P-value from the statistical test specified in `method.use` (default: Chi-square test).}
#'       \item{`effect.size`}{Effect size of the test, measured as Cramer's V.}
#'       \item{`approx`}{Indicates whether the Chi-square approximation is valid ("valid") or potentially unreliable ("warning"), based on whether at least 80% of transcript counts of the contingency table exceed 5.}
#'     }
#'   }
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr unite
#' @importFrom purrr reduce
#' @importFrom S4Vectors metadata metadata<-

RunDIU <- function(
    object,
    group.by = NULL,
    group.1 = NULL,
    group.2 = NULL,
    assay.use = "counts",
    method.use = "Chisq",
    min.gene.pct = 0.05,
    min.gene.cts = 15,
    min.tx.cts = 1,
    genes = NULL,
    only.valid = FALSE, # if TRUE and method.use is Chisq, removes genes that do not meet sample size for adequate approximation
    simulate.p = FALSE, # Fisher's tests will always be simulated with Monte Carlo
    p.adj = "BH",
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
  assertChoice(assay.use, assayNames(object))
  assertChoice(method.use, c("Chisq", "Fisher"))
  assertNumber(min.gene.pct, lower = 0, upper = 1, finite = TRUE)
  assertNumber(min.gene.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertCharacter(genes, unique = TRUE, null.ok = TRUE, any.missing = FALSE)
  assertFlag(only.valid)
  assertFlag(simulate.p)
  assertChoice(p.adj, c("BH", "Bonferroni"))
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

  # Group structure
  colData(object)$group_var <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE) %>%
    pull(group_var)
  unique_groups <- unique(colData(object)$group_var)

  if (length(unique_groups) < 2) {
    stop("There must be at least 2 groups to compare.")
  }

  ## check group subset
  assertSubset(group.1, choices = unique_groups, empty.ok = TRUE)
  assertSubset(group.2, choices = setdiff(unique_groups, group.1), empty.ok = TRUE)

  # Gene filter
  if (!is.null(genes)) {
    if (any(genes %in% unique(rowData(object)[[active.gene.id]]))) {
      missing_genes <- setdiff(genes, unique(rowData(object)[[active.gene.id]]))
      if (length(missing_genes) == length(genes)) {
        stop("None of the genes were found in the object. (Check active.gene.id?)")
      }
      if (length(missing_genes) > 0) {
        if (!quiet) message("\u2139 Warning: The following genes were not found in the object: '", paste0(missing_genes, collapse = "', '"), "'.")
        genes <- genes[genes %in% unique(rowData(object)[[active.gene.id]])]
      }
      object <- object[rowData(object)[[active.gene.id]] %in% genes, , drop = FALSE]
    } else {
      stop("None of the genes were found in the object. (Check active.gene.id?)")
    }
  }

  # DIU
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

    if (!quiet) message("Running DIU analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")

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

    if (!quiet) message("Running DIU analysis for ", group.1.updated, " vs all other cells...")
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

    if (!quiet) message("Running DIU analysis for ", group.1.updated, " vs ", group.2.updated, "...")

  }

  ## case: group 1 unspecified but group 2 is specified
  else if (is.null(group.1) && !is.null(group.2)) {
    stop("`group.1` must be specified prior to `group.2`")
  }


  # Loop through object grp list
  data_list <- list()
  stats_list <- list()

  for (comp in names(object_grp_list)) {

    if (!quiet && comp != "single_test" && length(unique_groups) > 2) message(comp, ": ")

    ## get group objects and names
    object_grp1 <- object_grp_list[[comp]]$grp1.object
    object_grp2 <- object_grp_list[[comp]]$grp2.object
    group.1 <- object_grp_list[[comp]]$grp1.names
    group.2 <- object_grp_list[[comp]]$grp2.names

    ## count mat for each group
    expr_mat_grp1 <- assay(object_grp1, assay.use)
    expr_mat_grp2 <- assay(object_grp2, assay.use)

    ## gene sums for each group
    gene_groups_grp1 <- rowData(object_grp1)[[active.gene.id]]
    gene_groups_grp2 <- rowData(object_grp2)[[active.gene.id]]
    expr_mat_gene_grp1 <- rowsum(expr_mat_grp1, group = gene_groups_grp1)
    expr_mat_gene_grp2 <- rowsum(expr_mat_grp2, group = gene_groups_grp2)
    gene_cts_grp1 <- rowSums(expr_mat_gene_grp1)
    gene_cts_grp2 <- rowSums(expr_mat_gene_grp2)

    ## gene detection rates for each group
    gene_pct_grp1 <- rowSums(expr_mat_gene_grp1 > 0) / ncol(expr_mat_gene_grp1)
    gene_pct_grp2 <- rowSums(expr_mat_gene_grp2 > 0) / ncol(expr_mat_gene_grp2)

    ## report
    gene_dr_df <- data.frame("gene.pct.grp1" = gene_pct_grp1,
                             "gene.pct.grp2" = gene_pct_grp2,
                             "gene.sum.grp1" = gene_cts_grp1,
                             "gene.sum.grp2" = gene_cts_grp2)

    ## filtering by gene detection and gene counts
    filt_genes <- gene_dr_df %>%
      filter(gene.pct.grp1 >= min.gene.pct &
                      gene.pct.grp2 >= min.gene.pct &
                      gene.sum.grp1 >= min.gene.cts &
                      gene.sum.grp2 >= min.gene.cts)
    filt_genes$gene.id <- rownames(filt_genes)

    ## filter genes from grp objects
    filt_object_grp1 <- object_grp1[rowData(object_grp1)[[active.gene.id]] %in% rownames(filt_genes), drop = FALSE]
    filt_object_grp2 <- object_grp2[rowData(object_grp2)[[active.gene.id]] %in% rownames(filt_genes), drop = FALSE]

    ## aggregate transcript counts
    agg_cts_df <- data.frame("gene.id.1" = rowData(filt_object_grp1)[[active.gene.id]],
                             "gene.id.2" = rowData(filt_object_grp2)[[active.gene.id]],
                             "gene.id" = rowData(filt_object_grp1)[[active.gene.id]],
                             "cts.1" = rowSums(assay(filt_object_grp1, assay.use)),
                             "cts.2" = rowSums(assay(filt_object_grp2, assay.use)))
    agg_cts_df <- agg_cts_df %>%
      select(-gene.id.1, -gene.id.2) %>%
      rownames_to_column(var = "transcript") %>%
      left_join(., filt_genes[, c("gene.pct.grp1", "gene.pct.grp2", "gene.id")], by = "gene.id")

    ## filter transcripts
    agg_cts_df <- agg_cts_df %>%
      filter(cts.1 >= min.tx.cts | cts.2 >= min.tx.cts)

    ## remove genes that have <2 isoforms
    agg_cts_df <- agg_cts_df %>%
      group_by(gene.id) %>%
      filter(n() > 1) %>%
      ungroup()

    ## number of tests to conduct
    n_tests <- length(unique(agg_cts_df$gene.id))
    if (!quiet) message(n_tests, " genes passed detection thresholds.")

    ## test sample size for Chisq approximation
    if (method.use == "Chisq") {
      agg_cts_df <- agg_cts_df %>%
        group_by(gene.id) %>%
        mutate(approx = ifelse(mean(c(cts.1, cts.2) > 5) > 0.80 & all(c(cts.1, cts.2) > 0), "valid", "warning")) %>%
        ungroup()

      warning_genes <- agg_cts_df %>%
        filter(approx == "warning") %>%
        pull(gene.id) %>%
        unique()

      if (!quiet && length(warning_genes) > 0) message("\u2139 Warning: ", length(warning_genes), " genes with inadequate sample size for Chisq approximation.")

      if (only.valid) {
        if (!quiet) message("\u2139 `only.valid` is set to TRUE. Only genes with valid approximations will be considered.")
        agg_cts_df <- agg_cts_df %>%
          filter(approx == "valid")
      }
    } else if (method.use == "Fisher") {
      agg_cts_df <- agg_cts_df %>%
        mutate(approx = NA)
    }

    ## number of tests to conduct after sample size assessment
    n_tests <- length(unique(agg_cts_df$gene.id))
    if (n_tests == 0) {
      next
    }
    if (!quiet) message("Testing DIU for ", n_tests, " genes...")

    ## proportion and delta
    diu_data <- agg_cts_df %>%
      group_by(gene.id) %>%
      mutate(prop.1 = cts.1 / sum(cts.1),
             prop.2 = cts.2 / sum(cts.2),
             dprop = prop.1 - prop.2) %>%
      ungroup()

    ## update list
    data_list[[comp]] <- diu_data %>%
      mutate(grp.1 = group.1,
             grp.2 = group.2) %>%
      select(grp.1, grp.2, gene.id, gene.pct.grp1, gene.pct.grp2, transcript,
             cts.1, cts.2, prop.1, prop.2, dprop) %>%
      rename("group.1" = grp.1,
                    "group.2" = grp.2,
                    "gene" = "gene.id",
                    "gene.pct.1" = "gene.pct.grp1",
                    "gene.pct.2" = "gene.pct.grp2",
                    "cts.1" = cts.1,
                    "cts.2" = cts.2,
                    "prop.1" = prop.1,
                    "prop.2" = prop.2,
                    "delta" = dprop) %>%
      arrange(group.1, gene)

    ## test statistics
    if (method.use == "Chisq") {
      if (!quiet && simulate.p) message("\u2139 p-values from Chi-square tests will be approximated by Monte Carlo simulation.")
      diu_stats <- diu_data %>%
        group_by(gene.id) %>%
        mutate(test_stats = list(suppressWarnings(chisq.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE), simulate.p.value = simulate.p))),
               pval = test_stats[[1]]$p.value,
               effect.size = sqrt(test_stats[[1]]$statistic / (sum(cts.1, cts.2) * 1))) %>%
        ungroup()
    } else if (method.use == "Fisher") {
      if (simulate.p == FALSE) {
        simulate.p <- TRUE
        if (!quiet) message("\u2139 p-values from Fisher's exact tests will be approximated by Monte Carlo simulation.")
      }
      diu_stats <- diu_data %>%
        group_by(gene.id) %>%
        mutate(test_stats = list(suppressWarnings(fisher.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE), simulate.p.value = TRUE))),
               pval = test_stats[[1]]$p.value,
               chisq_stat = suppressWarnings(chisq.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE), correct = FALSE)$statistic),
               effect.size = sqrt(chisq_stat / (sum(cts.1, cts.2) * 1))) %>%
        ungroup()
    }

    ## adjusted pval
    diu_stats <- diu_stats %>%
      group_by(gene.id) %>%
      slice_max(order_by = abs(dprop), with_ties = FALSE) %>%
      ungroup() %>%
      select(gene.id, dprop, transcript, pval, effect.size, approx) %>%
      mutate(padj = p.adjust(pval, method = p.adj))

    ## update list
    stats_list[[comp]] <- diu_stats %>%
      as.data.frame() %>%
      mutate(grp.1 = group.1,
             grp.2 = group.2) %>%
      rename("group.1" = grp.1,
                    "group.2" = grp.2,
                    "gene" = "gene.id",
                    "max.prop.diff" = "dprop") %>%
      select(group.1, group.2, gene, max.prop.diff, transcript, pval, padj, effect.size, approx) %>%
      arrange(padj)
  }

  # Output
  return_list <- list()
  if (length(data_list) > 0) {
    return_list$data <- as.data.frame(reduce(data_list, rbind))
  } else {
    return_list$data <- data.frame("group.1" = character(),
                                   "group.2" = character(),
                                   "gene" = character(),
                                   "gene.pct.1" = numeric(),
                                   "gene.pct.2" = numeric(),
                                   "transcript" = character(),
                                   "cts.1" = numeric(),
                                   "cts.2" = numeric(),
                                   "prop.1" = numeric(),
                                   "prop.2" = numeric(),
                                   "delta" = numeric())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind))
  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "max.delta" = numeric(),
                                    "pval" = numeric(),
                                    "padj" = character(),
                                    "effect.size" = character(),
                                    "approx" = numeric())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("There were 0 genes that passed detection thresholds for all comparisons.")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
