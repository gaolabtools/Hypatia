#' Run differential isoform usage analysis
#'
#' Tests genes for differential isoform usage between cell groups.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.1 Group label(s) for the first side of the comparison. If `NULL`, each group is compared against all others.
#' @param group.2 Optional group label(s) for the second side of the comparison. If `NULL`, `group.1` is compared against all other cells.
#' @param assay.use Assay name to use.
#' @param method.use Statistical test: `"Chisq"` or `"Fisher"`.
#' @param min.gene.pct Minimum fraction of cells in each group where the gene must be detected.
#' @param min.gene.cts Minimum total gene counts required in each group.
#' @param min.tx.cts Minimum transcript counts required for inclusion in contingency tables.
#' @param genes Optional vector of active gene IDs to test. Genes are still subject to filtering.
#' @param only.valid Logical; if `TRUE`, report only genes with valid Chi-square approximations.
#' @param simulate.p Logical; if `TRUE`, use Monte Carlo p-values. Fisher's exact test always uses simulation.
#' @param p.adj P-value adjustment method. Must be one of `stats::p.adjust.methods`.
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
#'       \item{`prop.diff`}{The difference in transcript proportions between groups (`group.1` - `group.2`).}
#'     }
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`max.prop.diff`}{The largest absolute difference in transcript proportions between `group.1` and `group.2`.}
#'       \item{`transcript`}{The transcript associated with `max.prop.diff`.}
#'       \item{`pval`}{P-value from the selected statistical test.}
#'       \item{`padj`}{Adjusted p-value.}
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
#' @importFrom purrr reduce
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats chisq.test fisher.test p.adjust

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
  p.adj <- .PAdjustMethod(p.adj)
  assertFlag(quiet)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  unique_groups <- unique(colData(object)$group_var)

  if (length(unique_groups) < 2) {
    stop("There must be at least 2 groups to compare.")
  }

  ## check group subset
  assertSubset(group.1, choices = unique_groups, empty.ok = TRUE)
  assertSubset(group.2, choices = setdiff(unique_groups, group.1), empty.ok = TRUE)

  # Gene filter
  if (!is.null(genes)) {
    gene_filter <- .FilterGenes(object, genes, active.gene.id, quiet = quiet)
    object <- gene_filter$object
    genes <- gene_filter$genes
  }

  # DIU
  object_grp_list <- .BuildGroupComparisons(object, group.1, group.2, unique_groups)
  comparison_mode <- attr(object_grp_list, "mode")
  if (!quiet && comparison_mode == "all") {
    message("Running DIU analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")
  } else if (!quiet && comparison_mode == "one_vs_all") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DIU analysis for ", comparison$grp1.names, " vs all other cells...")
  } else if (!quiet && comparison_mode == "pair") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DIU analysis for ", comparison$grp1.names, " vs ", comparison$grp2.names, "...")
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
    filt_object_grp1 <- object_grp1[rowData(object_grp1)[[active.gene.id]] %in% rownames(filt_genes), , drop = FALSE]
    filt_object_grp2 <- object_grp2[rowData(object_grp2)[[active.gene.id]] %in% rownames(filt_genes), , drop = FALSE]

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
    if (!quiet) message("  ", n_tests, " genes passed detection thresholds.")

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

      if (!quiet && length(warning_genes) > 0) message("  \u2139 Warning: ", length(warning_genes), " genes with inadequate sample size for Chisq approximation.")

      if (only.valid) {
        if (!quiet) message("  \u2139 `only.valid` is set to TRUE. Only genes with valid approximations will be considered.")
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
    if (!quiet) message("  Performing DIU comparisons...")

    ## proportion difference
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
                    "prop.diff" = dprop) %>%
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
                                   "prop.diff" = numeric())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind))
  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "max.prop.diff" = numeric(),
                                    "transcript" = character(),
                                    "pval" = numeric(),
                                    "padj" = numeric(),
                                    "effect.size" = numeric(),
                                    "approx" = character())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("There were 0 genes that passed detection thresholds for all comparisons.")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
