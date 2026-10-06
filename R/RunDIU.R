#' Run differential isoform usage analysis
#'
#' Tests genes for differential isoform usage between cell groups.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.1 Group label(s) for the first side of the comparison. If `NULL`, each group is compared against all others.
#' @param group.2 Optional group label(s) for the second side of the comparison. If `NULL`, `group.1` is compared against all other cells.
#' @param assay.use Assay name to use.
#' @param method.use Statistical test: uncorrected Pearson Chi-square (`"Chisq"`)
#'   or Fisher's exact test (`"Fisher"`).
#' @param min.gene.pct Minimum fraction of cells in each group where the gene must be detected.
#' @param min.gene.cts Minimum total gene counts required in each group.
#' @param min.tx.cts Minimum total transcript counts in at least one comparison group for inclusion in contingency tables.
#' @param genes Optional vector of active gene IDs to test. Genes are still subject to filtering.
#' @param only.valid Logical; if `TRUE`, report only genes with valid Chi-square approximations.
#' @param simulate.p Logical; if `TRUE`, request Monte Carlo p-values. Fisher's test always requests simulation, which `stats::fisher.test()` uses for tables larger than 2 by 2.
#' @param permutation Logical; if `TRUE`, calculate empirical p-values by permuting cell labels while preserving comparison-group sizes.
#' @param perm.iter Number of cell-label permutations used to construct the empirical null distribution.
#' @param perm.statistics Character vector of gene-level statistics to evaluate
#'   for each cell-label permutation. Options are the uncorrected Pearson
#'   Chi-square statistic (`"chisq"`) and the maximum absolute transcript
#'   proportion difference (`"max_delta"`).
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell transcript proportion dispersion among cells with positive retained-gene counts.
#' @param bootstrap Logical; if `TRUE`, bootstrap cells within each comparison group to estimate uncertainty in transcript proportion differences.
#' @param boot.iter Number of bootstrap iterations. Default: 5000.
#' @param boot.fraction Proportion of cells sampled from each comparison group per bootstrap iteration. Default: 1 (each group's original cell count).
#' @param boot.ncells Optional fixed number of cell draws sampled with replacement from each comparison group per bootstrap iteration. Cannot be supplied together with `boot.fraction`.
#' @param boot.conf Nominal confidence level for percentile bootstrap intervals. Default: 0.95.
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
#'       \item{`gene.pct.1`}{Fraction of cells in `group.1` with expression of the gene, on the `[0, 1]` scale.}
#'       \item{`gene.pct.2`}{Fraction of cells in `group.2` with expression of the gene, on the `[0, 1]` scale.}
#'       \item{`transcript`}{The associated transcript.}
#'       \item{`cts.1`}{Total counts of the transcript across all cells in `group.1`.}
#'       \item{`cts.2`}{Total counts of the transcript across all cells in `group.2`.}
#'       \item{`prop.1`}{Transcript proportion for `group.1`.}
#'       \item{`prop.2`}{Transcript proportion for `group.2`.}
#'       \item{`prop.diff`}{The difference in transcript proportions between groups (`group.1` - `group.2`).}
#'       \item{`cell.n.1`, `cell.n.2`}{Number of gene-positive cells used for cell-level summaries in each group.}
#'       \item{`cell.prop.mean.1`, `cell.prop.mean.2`, `cell.prop.median.1`, `cell.prop.median.2`, `cell.prop.sd.1`, `cell.prop.sd.2`, `cell.prop.iqr.1`, `cell.prop.iqr.2`}{Mean, median, sample standard deviation, and interquartile range of cell-level transcript proportions.}
#'       \item{`cell.prop.zero.frac.1`, `cell.prop.zero.frac.2`}{Fraction of gene-positive cells in which the transcript is not detected.}
#'       \item{`boot.prop.1`, `boot.prop.2`, `boot.prop.diff`}{List-columns containing bootstrap proportions and differences. Each entry is an empty numeric vector when `bootstrap = FALSE`.}
#'     }
#'     Cell-dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`max.prop.diff`}{The signed transcript-proportion difference (`group.1` - `group.2`) for the transcript with the greatest absolute difference.}
#'       \item{`transcript`}{The transcript associated with `max.prop.diff`.}
#'       \item{`pval`}{P-value from the selected statistical test. Chi-square
#'       tests use the uncorrected Pearson statistic.}
#'       \item{`padj`}{Adjusted p-value, calculated separately within each comparison (default: Benjamini-Hochberg).}
#'       \item{`pval.perm`}{Empirical p-value from cell-label permutation using
#'       the uncorrected Pearson Chi-square statistic. Values are `NA` when
#'       `permutation = FALSE` or `"chisq"` is not requested.}
#'       \item{`padj.perm`}{Adjusted empirical Chi-square permutation p-value.}
#'       \item{`pval.perm.delta`}{Empirical gene-level p-value from cell-label
#'       permutation using the maximum absolute transcript proportion difference.
#'       Values are `NA` when `permutation = FALSE` or `"max_delta"` is not
#'       requested.}
#'       \item{`padj.perm.delta`}{Adjusted empirical maximum-delta permutation
#'       p-value.}
#'       \item{`cramers.v`}{Cramer's V calculated from the uncorrected Pearson
#'       Chi-square statistic.}
#'       \item{`approx`}{A count-based heuristic: `"valid"` when more than 80% of observed contingency-table counts exceed 5 and all are positive; `"warning"` otherwise. This does not check expected counts. `NA` for Fisher's test.}
#'       \item{`boot.prop.diff.mean`, `boot.prop.diff.lower`, `boot.prop.diff.upper`}{Bootstrap mean and percentile interval for the signed proportion difference of the transcript associated with `max.prop.diff`. Values are `NA` when `bootstrap = FALSE` or no finite differences are available. Interval bounds are also `NA` when fewer than two finite differences are available.}
#'       \item{`boot.valid.iter`}{Number of finite bootstrap differences used for the summaries.}
#'     }
#'   }
#' }
#' Rows in `$data` are ordered by `group.1`, `group.2`, `gene`, and
#' `transcript`. Rows in `$stats` are ordered by `group.1`, `group.2`,
#' `padj`, `gene`, and `transcript`.
#' @details Genes must meet `min.gene.pct` and `min.gene.cts` in both comparison
#' groups. Transcripts are retained when their total counts meet `min.tx.cts`
#' in either group, and a gene must retain at least two transcripts. Counts are
#' pooled into a transcript-by-group contingency table for each gene.
#'
#' - `"Chisq"`: Pearson Chi-square testing without continuity correction.
#'   `simulate.p = TRUE` requests Monte Carlo p-values. `only.valid = TRUE`
#'   restricts testing to genes passing the observed-count heuristic reported
#'   in `approx`.
#' - `"Fisher"`: Fisher's test, requesting Monte Carlo p-values for tables
#'   larger than 2 by 2. Two-by-two tables use the exact calculation.
#'   `only.valid` is not applied with this method.
#'
#' With `permutation = TRUE`, cell labels are shuffled within each comparison
#' while preserving group sizes and the retained transcript set. `perm.statistics`
#' selects the Pearson Chi-square statistic, the maximum absolute transcript
#' proportion difference, or both. Permutation p-values are returned separately
#' from the selected contingency-table test's p-values.
#'
#' With `bootstrap = TRUE`, cells are sampled with replacement within each
#' group. Each iteration draws `ceiling(boot.fraction * ncol(group))` cells,
#' or `boot.ncells` when supplied. The defaults use 5000 iterations and each
#' group's original sample size (`boot.fraction = 1`). Differences are computed
#' within matching iterations; only finite differences enter the summaries.
#' Interval bounds are unavailable with fewer than two finite differences.
#' With `boot.conf = 0.95`, bounds are the 2.5th and 97.5th percentiles.
#'
#' Percentile intervals summarize the signed proportion difference for the
#' transcript selected from the original data as having the largest absolute
#' difference. This transcript remains fixed across bootstrap iterations. The
#' interval is not for the gene-level maximum absolute difference, is not
#' adjusted for transcript selection, and does not provide simultaneous
#' coverage across transcripts or genes. `max.prop.diff` retains its sign.
#' Confidence levels are nominal; sparse genes, biased estimates, and many
#' unavailable differences can affect interval coverage. Smaller fractions or
#' fixed `boot.ncells` change the resampling sample size, so their unadjusted
#' percentiles are not automatically calibrated as confidence intervals for the
#' original sample. Cells are resampled independently; donor-level or other
#' clustered study designs require resampling that accounts for those units.
#' Bootstrap summaries are separate from the theoretical and permutation
#' p-values. Optional `cell.dispersion` summaries describe proportions
#' among cells with positive retained-gene counts, separately from bootstrapping.
#'
#' Multiple `group.by` columns are joined with `_`. Omitting both comparison
#' arguments compares each group with the remaining cells; exactly two groups
#' produce one comparison. Multiple labels in `group.1` or `group.2` pool cells
#' on that side. P-values are adjusted separately within each comparison and
#' each test family, using Benjamini-Hochberg by default. Call `set.seed()` for
#' reproducible resampling results; no seed is set internally. An error is
#' returned if no genes pass filtering across all comparisons.
#' @seealso [GetUsage()], [PlotUsage()]
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
    simulate.p = FALSE, # Fisher's tests request Monte Carlo for tables larger than 2 by 2
    permutation = FALSE,
    perm.iter = 1000,
    perm.statistics = c("chisq", "max_delta"),
    bootstrap = TRUE,
    boot.iter = 5000,
    boot.fraction = 1,
    boot.ncells = NULL,
    boot.conf = 0.95,
    p.adj = "BH",
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
  assertChoice(assay.use, assayNames(object))
  assertChoice(method.use, c("Chisq", "Fisher"))
  assertNumber(min.gene.pct, lower = 0, upper = 1, finite = TRUE)
  assertNumber(min.gene.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertCharacter(genes, unique = TRUE, null.ok = TRUE, any.missing = FALSE)
  assertFlag(only.valid)
  assertFlag(simulate.p)
  assertFlag(permutation)
  assertCount(perm.iter, positive = TRUE)
  assertCharacter(perm.statistics, min.len = 1, unique = TRUE, any.missing = FALSE)
  assertSubset(perm.statistics, c("chisq", "max_delta"), empty.ok = FALSE)
  assertFlag(cell.dispersion)
  assertFlag(bootstrap)
  assertCount(boot.iter, positive = TRUE)
  boot.fraction.supplied <- !missing(boot.fraction)
  assertNumber(boot.fraction, lower = 0.01, upper = 1, finite = TRUE)
  assertCount(boot.ncells, positive = TRUE, null.ok = TRUE)
  if (!is.null(boot.ncells) && boot.fraction.supplied) {
    stop("Please provide only one of `boot.fraction` or `boot.ncells`.", call. = FALSE)
  }
  assertNumber(boot.conf, lower = 0, upper = 1, finite = TRUE)
  if (boot.conf == 0 || boot.conf == 1) {
    stop("`boot.conf` must be strictly between 0 and 1.", call. = FALSE)
  }
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
    if (!quiet) message("  ", n_tests, " genes passed filtering.")

    ## test sample size for Chisq approximation
    if (method.use == "Chisq") {
      agg_cts_df <- agg_cts_df %>%
        group_by(gene.id) %>%
        mutate(approx = ifelse(mean(c(cts.1, cts.2) > 5) > 0.80 & all(c(cts.1, cts.2) > 0), "valid", "warning")) %>%
        ungroup()

      if (only.valid) {
        if (!quiet) message("  \u2139 `only.valid` is set to TRUE. Only genes with valid approximations will be considered.")
        agg_cts_df <- agg_cts_df %>%
          filter(approx == "valid")
      }
    } else if (method.use == "Fisher") {
      agg_cts_df <- agg_cts_df %>%
        mutate(approx = NA_character_)
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

    ## test statistics
    if (method.use == "Chisq") {
      if (!quiet && simulate.p) message("\u2139 p-values from Chi-square tests will be approximated by Monte Carlo simulation.")
      diu_stats <- diu_data %>%
        group_by(gene.id) %>%
        mutate(test_stats = list(suppressWarnings(chisq.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE), correct = FALSE, simulate.p.value = simulate.p))),
               pval = test_stats[[1]]$p.value,
               cramers.v = sqrt(test_stats[[1]]$statistic / (sum(cts.1, cts.2) * 1))) %>%
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
               cramers.v = sqrt(chisq_stat / (sum(cts.1, cts.2) * 1))) %>%
        ungroup()
    }

    ## Cell-label permutation p-values. The observed filtering and transcript
    ## features are fixed; only cell labels are permuted within each comparison.
    if (permutation) {
      if (!quiet) message("  Performing cell-label permutations (", perm.iter, " iterations)...")

      perm_transcripts <- diu_data$transcript
      perm_object <- object[, c(colnames(object_grp1), colnames(object_grp2)), drop = FALSE]
      perm_expr <- assay(perm_object[perm_transcripts, , drop = FALSE], assay.use)
      perm_gene_ids <- rowData(perm_object[perm_transcripts, , drop = FALSE])[[active.gene.id]]
      perm_grp1_ncells <- ncol(object_grp1)
      perm_gene_ids_unique <- unique(perm_gene_ids)
      perm_gene_factor <- factor(perm_gene_ids, levels = perm_gene_ids_unique)
      perm_total_tx <- rowSums(perm_expr)
      perm_total_gene <- as.numeric(rowsum(perm_total_tx, perm_gene_factor))

      permutation_statistics <- function(grp1_idx) {
        tx_grp1 <- rowSums(perm_expr[, grp1_idx, drop = FALSE])
        gene_grp1 <- as.numeric(rowsum(tx_grp1, perm_gene_factor))
        gene_grp2 <- perm_total_gene - gene_grp1
        tx_grp2 <- perm_total_tx - tx_grp1

        statistics <- list()

        if ("chisq" %in% perm.statistics) {
          expected_grp1 <- perm_total_tx * gene_grp1[perm_gene_factor] /
            perm_total_gene[perm_gene_factor]
          expected_grp2 <- perm_total_tx * gene_grp2[perm_gene_factor] /
            perm_total_gene[perm_gene_factor]
          valid <- expected_grp1 > 0 & expected_grp2 > 0
          contribution <- rep(NA_real_, length(tx_grp1))
          contribution[valid] <- (tx_grp1[valid] - expected_grp1[valid])^2 /
            expected_grp1[valid] +
            (tx_grp2[valid] - expected_grp2[valid])^2 / expected_grp2[valid]
          chisq_statistic <- rowsum(
            contribution, perm_gene_factor, na.rm = FALSE
          )[, 1]
          names(chisq_statistic) <- perm_gene_ids_unique
          statistics$chisq <- chisq_statistic
        }

        if ("max_delta" %in% perm.statistics) {
          gene_total_grp1 <- gene_grp1[perm_gene_factor]
          gene_total_grp2 <- gene_grp2[perm_gene_factor]
          prop_diff <- tx_grp1 / gene_total_grp1 - tx_grp2 / gene_total_grp2
          max_delta_statistic <- vapply(
            split(abs(prop_diff), perm_gene_factor),
            max,
            numeric(1)
          )
          statistics$max_delta <- max_delta_statistic[perm_gene_ids_unique]
        }

        statistics
      }

      observed_perm_stats <- permutation_statistics(seq_len(perm_grp1_ncells))
      null_stats <- lapply(perm.statistics, function(statistic) {
        matrix(
          NA_real_,
          nrow = length(perm_gene_ids_unique),
          ncol = perm.iter,
          dimnames = list(perm_gene_ids_unique, NULL)
        )
      })
      names(null_stats) <- perm.statistics
      for (iter in seq_len(perm.iter)) {
        permuted_idx <- sample.int(ncol(perm_expr))
        perm_grp1_idx <- permuted_idx[seq_len(perm_grp1_ncells)]
        iter_stats <- permutation_statistics(perm_grp1_idx)
        for (statistic in perm.statistics) {
          null_stats[[statistic]][, iter] <- iter_stats[[statistic]]
        }
      }

      perm_pvals <- lapply(perm.statistics, function(statistic) {
        observed_stats <- observed_perm_stats[[statistic]]
        pvals <- vapply(seq_along(observed_stats), function(i) {
          observed <- observed_stats[[i]]
          null <- null_stats[[statistic]][i, ]
          if (!is.finite(observed)) return(NA_real_)
          null <- null[is.finite(null)]
          if (length(null) == 0L) return(NA_real_)
          (1 + sum(null >= observed)) / (1 + length(null))
        }, numeric(1))
        names(pvals) <- perm_gene_ids_unique
        pvals
      })
      names(perm_pvals) <- perm.statistics
    } else {
      perm_pvals <- NULL
    }

    ## adjusted pval
    diu_stats <- diu_stats %>%
      group_by(gene.id) %>%
      slice_max(order_by = abs(dprop), with_ties = FALSE) %>%
      ungroup() %>%
      select(gene.id, dprop, transcript, pval, cramers.v, approx) %>%
      mutate(
        padj = p.adjust(pval, method = p.adj),
        pval.perm = if (permutation && "chisq" %in% perm.statistics) {
          unname(perm_pvals$chisq[gene.id])
        } else {
          NA_real_
        },
        padj.perm = if (permutation && "chisq" %in% perm.statistics) {
          p.adjust(pval.perm, method = p.adj)
        } else {
          NA_real_
        },
        pval.perm.delta = if (permutation && "max_delta" %in% perm.statistics) {
          unname(perm_pvals$max_delta[gene.id])
        } else {
          NA_real_
        },
        padj.perm.delta = if (permutation && "max_delta" %in% perm.statistics) {
          p.adjust(.data$pval.perm.delta, method = p.adj)
        } else {
          NA_real_
        }
      )

    ## bootstrap transcript proportions
    if (bootstrap) {
      if (!quiet) message("  Performing bootstrap sampling...")

      boot_transcripts <- diu_data$transcript
      boot_expr_mat_grp1 <- assay(object_grp1[boot_transcripts, , drop = FALSE], assay.use)
      boot_expr_mat_grp2 <- assay(object_grp2[boot_transcripts, , drop = FALSE], assay.use)
      boot_gene_ids <- rowData(object_grp1[boot_transcripts, , drop = FALSE])[[active.gene.id]]

      if (is.null(boot.ncells)) {
        grp1_boot_ncells <- ceiling(ncol(boot_expr_mat_grp1) * boot.fraction)
        grp2_boot_ncells <- ceiling(ncol(boot_expr_mat_grp2) * boot.fraction)
      } else {
        grp1_boot_ncells <- boot.ncells
        grp2_boot_ncells <- boot.ncells
      }

      ## Pool counts directly to avoid rebuilding grouped data frames per draw.
      boot_gene_factor <- factor(boot_gene_ids, levels = unique(boot_gene_ids))
      boot_gene_index <- as.integer(boot_gene_factor)
      boot_prop_1 <- matrix(NA_real_, nrow = length(boot_transcripts), ncol = boot.iter)
      boot_prop_2 <- matrix(NA_real_, nrow = length(boot_transcripts), ncol = boot.iter)
      for (b in seq_len(boot.iter)) {
        grp1_col_idx <- sample(
          seq_len(ncol(boot_expr_mat_grp1)),
          size = grp1_boot_ncells,
          replace = TRUE
        )
        grp2_col_idx <- sample(
          seq_len(ncol(boot_expr_mat_grp2)),
          size = grp2_boot_ncells,
          replace = TRUE
        )

        cts.1 <- rowSums(boot_expr_mat_grp1[, grp1_col_idx, drop = FALSE])
        cts.2 <- rowSums(boot_expr_mat_grp2[, grp2_col_idx, drop = FALSE])
        gene_cts.1 <- as.numeric(base::rowsum(cts.1, boot_gene_factor, reorder = FALSE))
        gene_cts.2 <- as.numeric(base::rowsum(cts.2, boot_gene_factor, reorder = FALSE))
        boot_prop_1[, b] <- cts.1 / gene_cts.1[boot_gene_index]
        boot_prop_2[, b] <- cts.2 / gene_cts.2[boot_gene_index]
      }
      boot_prop_diff <- boot_prop_1 - boot_prop_2

      boot_data <- data.frame(transcript = boot_transcripts)
      boot_data$boot.prop.1 <- lapply(seq_along(boot_transcripts), function(i) boot_prop_1[i, ])
      boot_data$boot.prop.2 <- lapply(seq_along(boot_transcripts), function(i) boot_prop_2[i, ])
      boot_data$boot.prop.diff <- lapply(seq_along(boot_transcripts), function(i) boot_prop_diff[i, ])
    } else {
      boot_data <- data.frame(transcript = diu_data$transcript)
      boot_data$boot.prop.1 <- rep(list(numeric(0)), nrow(boot_data))
      boot_data$boot.prop.2 <- rep(list(numeric(0)), nrow(boot_data))
      boot_data$boot.prop.diff <- rep(list(numeric(0)), nrow(boot_data))
    }

    ## cell-to-cell transcript proportion dispersion
    dispersion_transcripts <- diu_data$transcript
    if (cell.dispersion) {
      if (!quiet) message("  Calculating cell-to-cell dispersion...")
      dispersion_object_grp1 <- object_grp1[dispersion_transcripts, , drop = FALSE]
      dispersion_object_grp2 <- object_grp2[dispersion_transcripts, , drop = FALSE]
      dispersion_gene_ids <- rowData(dispersion_object_grp1)[[active.gene.id]]
      dispersion_data_grp1 <- .CellUsageDispersion(
        assay(dispersion_object_grp1, assay.use), dispersion_gene_ids
      )
      dispersion_data_grp2 <- .CellUsageDispersion(
        assay(dispersion_object_grp2, assay.use), dispersion_gene_ids
      )
      names(dispersion_data_grp1)[-1] <- paste0(names(dispersion_data_grp1)[-1], ".1")
      names(dispersion_data_grp2)[-1] <- paste0(names(dispersion_data_grp2)[-1], ".2")
      dispersion_data <- left_join(
        dispersion_data_grp1, dispersion_data_grp2, by = "transcript"
      )
    } else {
      dispersion_data <- data.frame(
        transcript = dispersion_transcripts,
        cell.n.1 = rep(NA_integer_, length(dispersion_transcripts)),
        cell.prop.mean.1 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.median.1 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.sd.1 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.iqr.1 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.zero.frac.1 = rep(NA_real_, length(dispersion_transcripts)),
        cell.n.2 = rep(NA_integer_, length(dispersion_transcripts)),
        cell.prop.mean.2 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.median.2 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.sd.2 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.iqr.2 = rep(NA_real_, length(dispersion_transcripts)),
        cell.prop.zero.frac.2 = rep(NA_real_, length(dispersion_transcripts))
      )
    }

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
      left_join(boot_data, by = "transcript") %>%
      left_join(dispersion_data, by = "transcript") %>%
      arrange(group.1, gene)

    ## bootstrap summaries for the observed maximum-difference transcript
    selected_boot_diff <- boot_data$boot.prop.diff[
      match(diu_stats$transcript, boot_data$transcript)
    ]
    selected_boot_diff <- lapply(selected_boot_diff, function(x) x[is.finite(x)])
    boot_valid_iter <- lengths(selected_boot_diff)
    boot_prop_diff_mean <- vapply(selected_boot_diff, function(x) {
      if (length(x) == 0) NA_real_ else mean(x)
    }, numeric(1))
    boot_alpha <- (1 - boot.conf) / 2
    boot_prop_diff_lower <- vapply(selected_boot_diff, function(x) {
      if (length(x) < 2) NA_real_ else unname(stats::quantile(x, probs = boot_alpha))
    }, numeric(1))
    boot_prop_diff_upper <- vapply(selected_boot_diff, function(x) {
      if (length(x) < 2) NA_real_ else unname(stats::quantile(x, probs = 1 - boot_alpha))
    }, numeric(1))

    diu_stats <- diu_stats %>%
      mutate(
        boot.prop.diff.mean = boot_prop_diff_mean,
        boot.prop.diff.lower = boot_prop_diff_lower,
        boot.prop.diff.upper = boot_prop_diff_upper,
        boot.valid.iter = as.integer(boot_valid_iter)
      )

    ## update list
    stats_list[[comp]] <- diu_stats %>%
      as.data.frame() %>%
      mutate(grp.1 = group.1,
             grp.2 = group.2) %>%
      rename("group.1" = grp.1,
                    "group.2" = grp.2,
                    "gene" = "gene.id",
                    "max.prop.diff" = "dprop") %>%
      select(group.1, group.2, gene, max.prop.diff, transcript, pval, padj,
             pval.perm, padj.perm, "pval.perm.delta", "padj.perm.delta",
             cramers.v, approx, boot.prop.diff.mean,
             boot.prop.diff.lower, boot.prop.diff.upper, boot.valid.iter) %>%
      arrange(padj)
  }

  # Output
  return_list <- list()
  if (length(data_list) > 0) {
    return_list$data <- as.data.frame(reduce(data_list, rbind)) %>%
      arrange(group.1, group.2, gene, transcript)
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
                                   "prop.diff" = numeric(),
                                   "boot.prop.1" = I(list()),
                                   "boot.prop.2" = I(list()),
                                   "boot.prop.diff" = I(list()),
                                   "cell.n.1" = integer(),
                                   "cell.prop.mean.1" = numeric(),
                                   "cell.prop.median.1" = numeric(),
                                   "cell.prop.sd.1" = numeric(),
                                   "cell.prop.iqr.1" = numeric(),
                                   "cell.prop.zero.frac.1" = numeric(),
                                   "cell.n.2" = integer(),
                                   "cell.prop.mean.2" = numeric(),
                                   "cell.prop.median.2" = numeric(),
                                   "cell.prop.sd.2" = numeric(),
                                   "cell.prop.iqr.2" = numeric(),
                                   "cell.prop.zero.frac.2" = numeric())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind)) %>%
      arrange(group.1, group.2, padj, gene, transcript)
  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "max.prop.diff" = numeric(),
                                    "transcript" = character(),
                                    "pval" = numeric(),
                                    "padj" = numeric(),
                                    "pval.perm" = numeric(),
                                    "padj.perm" = numeric(),
                                    "pval.perm.delta" = numeric(),
                                    "padj.perm.delta" = numeric(),
                                    "cramers.v" = numeric(),
                                    "approx" = character(),
                                    "boot.prop.diff.mean" = numeric(),
                                    "boot.prop.diff.lower" = numeric(),
                                    "boot.prop.diff.upper" = numeric(),
                                    "boot.valid.iter" = integer())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("There were 0 genes that passed detection thresholds for all comparisons.")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
