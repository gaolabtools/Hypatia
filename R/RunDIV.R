#' Run isoform diversity analysis
#'
#' Assesses differential isoform diversity between cell groups using directional
#' bootstrap support for an effect threshold.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.1 Group label(s) for the first side of the comparison. If `NULL`, each group is compared against all others.
#' @param group.2 Optional group label(s) for the second side of the comparison. If `NULL`, `group.1` is compared against all other cells.
#' @param entropy.use Diversity index: `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `"InverseSimpson"`.
#' @param assay.use Assay name to use.
#' @param entropy.thresh Diversity index threshold used to classify genes as monoform or polyform. If `NULL`, a default is chosen from the entropy index. Default thresholds for Tsallis and Renyi are defined only at orders 3 and 2, respectively; other orders return `NA` classifications unless a threshold is supplied.
#' @param prop.thresh Minimum within-gene transcript proportion used to define an effective isoform. Transcripts with proportions greater than or equal to this value are effective.
#' @param min.gene.pct Minimum fraction of cells in each group where the gene must be detected.
#' @param min.gene.cts Minimum total gene counts required in each group.
#' @param min.tx.cts Minimum total transcript counts in at least one comparison group for inclusion in diversity calculations. Effective isoform counts apply this threshold separately within each group.
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell diversity among cells with positive retained-gene counts.
#' @param boot.iter Number of bootstrap iterations to perform. Default: 5000.
#' @param boot.fraction Proportion of cells sampled from each comparison group per bootstrap iteration. Default: 1 (each group's original cell count).
#' @param boot.ncells Optional fixed number of cell draws sampled with replacement from each comparison group per bootstrap iteration. Cannot be supplied together with `boot.fraction`.
#' @param include.single Logical; if `FALSE`, genes with only one associated transcript after filtering will be excluded from the analysis.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi. At order 1, Tsallis and Renyi use their Shannon entropy limit, and NormalizedRenyi uses normalized Shannon entropy.
#' @param genes Optional vector of active gene IDs to test. Genes are still subject to filtering.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param top.n Optional number of the most abundant isoforms to include in diversity calculations. If `NULL`, all isoforms are included. Must be at least 2 when supplied.
#' @param renormalize Logical; if `TRUE`, rescale the selected isoform proportions to sum to one before calculating diversity.
#' @param div.diff.thresh Nonnegative minimum diversity difference for bootstrap support. Differences must strictly exceed this threshold in magnitude. Default: 0.10; choose a value appropriate for the selected entropy scale. Independent of `entropy.thresh`, which controls monoform/polyform classification.
#' @param support.thresh Directional bootstrap support cutoff, greater than 0.5 and at most 1. Default: 0.975.
#'
#' @returns A list containing two data frames:
#' \describe{
#'   \item{`$data`}{
#'     A data frame of summarized data with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`gene.pct.1`}{Fraction of cells in `group.1` with expression of the gene, on the `[0, 1]` scale.}
#'       \item{`gene.pct.2`}{Fraction of cells in `group.2` with expression of the gene, on the `[0, 1]` scale.}
#'       \item{`n.transcripts`}{Number of transcripts retained for the gene after count filtering, before optional `top.n` selection.}
#'       \item{`div.1`}{A list-column containing bootstrapped isoform diversity values of each gene for `group.1`.}
#'       \item{`div.2`}{A list-column containing bootstrapped isoform diversity values of each gene for `group.2`.}
#'       \item{`cell.n.1`, `cell.n.2`}{Number of gene-positive cells used for cell-level summaries in each group.}
#'       \item{`cell.div.mean.1`, `cell.div.mean.2`, `cell.div.median.1`, `cell.div.median.2`, `cell.div.sd.1`, `cell.div.sd.2`, `cell.div.iqr.1`, `cell.div.iqr.2`}{Mean, median, sample standard deviation, and interquartile range of cell-level diversity.}
#'     }
#'     Cell-dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#'     Includes NA entries derived from bootstrapped sampling that resulted in 0 gene counts.
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`avgDiv.1`}{Average of bootstrapped isoform diversity of the gene for cells in `group.1`.}
#'       \item{`avgDiv.2`}{Average of bootstrapped isoform diversity of the gene for cells in `group.2`.}
#'       \item{`div.diff`}{Mean finite within-iteration difference (`group.1` - `group.2`). Equals the difference of group-wise bootstrap means when both sides are finite in every iteration.}
#'       \item{`div.diff.lower`, `div.diff.upper`}{The 2.5th and 97.5th percentiles of finite within-iteration diversity differences (95% percentile interval). `NA_real_` if fewer than two differences are finite.}
#'       \item{`support.positive`, `support.negative`}{Fractions of finite differences strictly above `div.diff.thresh` and strictly below `-div.diff.thresh`, respectively. `NA_real_` if no differences are finite.}
#'       \item{`boot.valid.iter`}{Number of iterations with finite diversity values on both sides.}
#'       \item{`supported`}{Logical flag indicating that either directional support reaches `support.thresh`. `NA` if fewer than two differences are finite. This flag does not control multiple testing.}
#'       \item{`n.effective.1`, `n.effective.2`}{Number of transcripts with within-gene proportion greater than or equal to `prop.thresh` in each group.}
#'       \item{`div.class.1`}{`"monoform"` when `avgDiv.1` is at or below `entropy.thresh` and `"polyform"` otherwise.}
#'       \item{`div.class.2`}{`"monoform"` when `avgDiv.2` is at or below `entropy.thresh` and `"polyform"` otherwise.}
#'     }
#'   }
#' }
#' Rows in `$data` are ordered by `group.1`, `group.2`, and `gene`. Rows in
#' `$stats` are ordered by `group.1`, `group.2`, and `gene`.
#' @details Genes must meet `min.gene.pct` and `min.gene.cts` in both comparison
#' groups. Transcripts are retained when their total counts meet `min.tx.cts`
#' in either group. This retained set is used for bootstrapping; `include.single`
#' controls whether genes with only one retained transcript are tested.
#'
#' Cells are sampled with replacement within each group for `boot.iter`
#' iterations. Each iteration draws `ceiling(boot.fraction * ncol(group))` cells,
#' or `boot.ncells` when supplied, and calculates diversity from their pooled
#' transcript counts. The entropy indices and optional `top.n` selection and
#' renormalization follow [GetDiversity()]. `avgDiv` averages available bootstrap
#' values, and classification applies `entropy.thresh` to that average.
#'
#' Bootstrap support subtracts group 2 diversity from group 1
#' diversity within each iteration. Only iterations with finite values on both
#' sides enter the support fractions and percentile interval. Genes with no
#' valid differences remain in the output with unavailable support. `supported`
#' is `TRUE` when at least two differences are finite and either directional
#' support is at least `support.thresh`. Results are returned for all retained
#' genes; filter `$stats` by `supported` to select bootstrap-supported effects.
#' Support fractions describe resampling stability, not posterior probabilities
#' or adjusted p-values.
#' The default 0.975 directional cutoff approximately corresponds to a 95%
#' percentile interval wholly outside the effect band when bootstrap interval
#' coverage is valid. No false discovery rate control is implied.
#'
#' A conventional cell bootstrap uses each group's original sample size
#' (`boot.fraction = 1`). Smaller fractions or fixed `boot.ncells` change the
#' resampling sample size; their percentile intervals are not automatically
#' calibrated as confidence intervals for the original sample. This function
#' resamples cells independently within groups; donor-level or other clustered
#' study designs require a resampling procedure accounting for those units.
#'
#' Effective isoform counts are calculated from the original pooled counts,
#' applying `min.tx.cts` separately in each group, then counting proportions
#' at or above `prop.thresh`. They are independent of bootstrap sampling and
#' `top.n`. Optional `cell.dispersion` summaries describe diversity among cells
#' with positive retained-gene counts, rather than bootstrap uncertainty.
#'
#' Multiple `group.by` columns are joined with `_`. Omitting both comparison
#' arguments compares each group with the remaining cells; exactly two groups
#' produce one comparison. Multiple labels on either side pool their cells.
#' Call `set.seed()` for reproducible results; no seed is set internally.
#' An error is returned if no genes pass filtering across all comparisons.
#' @seealso [GetDiversity()], [PlotDiversity()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom purrr reduce

RunDIV <- function (
    object,
    group.by = NULL,
    group.1 = NULL,
    group.2 = NULL,
    entropy.use = "Tsallis",
    assay.use = "counts",
    entropy.thresh = NULL,
    prop.thresh = 0.2,
    min.gene.pct = 0.05,
    min.gene.cts = 15,
    min.tx.cts = 1,
    boot.iter = 5000,
    boot.fraction = 1,
    boot.ncells = NULL,
    include.single = TRUE,
    order = NULL,
    genes = NULL,
    quiet = FALSE,
    cell.dispersion = FALSE,
    top.n = NULL,
    renormalize = FALSE,
    div.diff.thresh = 0.10,
    support.thresh = 0.975
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.1, null.ok = TRUE)
  assertCharacter(group.2, null.ok = TRUE)
  assertChoice(entropy.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(entropy.thresh, lower = 0, finite = TRUE, null.ok = TRUE)
  assertNumber(prop.thresh, lower = 0, upper = 1, finite = TRUE)
  if (prop.thresh == 0) {
    stop("`prop.thresh` must be greater than 0.", call. = FALSE)
  }
  assertNumber(min.gene.pct, lower = 0, upper = 1, finite = TRUE)
  assertNumber(min.gene.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertFlag(cell.dispersion)
  assertNumber(div.diff.thresh, lower = 0, finite = TRUE)
  assertNumber(support.thresh, lower = 0.5, upper = 1, finite = TRUE)
  if (support.thresh <= 0.5) {
    stop("`support.thresh` must be greater than 0.5.", call. = FALSE)
  }
  boot.fraction.supplied <- !missing(boot.fraction)
  assertCount(boot.iter, positive = TRUE)
  assertNumber(boot.fraction, lower = 0.01, upper = 1, finite = TRUE)
  assertCount(boot.ncells, positive = TRUE, null.ok = TRUE)
  if (!is.null(boot.ncells) && boot.fraction.supplied) {
    stop("Please provide only one of `boot.fraction` or `boot.ncells`.", call. = FALSE)
  }
  assertFlag(include.single)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertCharacter(genes, null.ok = TRUE, any.missing = FALSE, unique = TRUE)
  assertFlag(quiet)
  assertCount(top.n, positive = TRUE, null.ok = TRUE)
  if (!is.null(top.n) && top.n < 2) {
    stop("`top.n` must be at least 2.", call. = FALSE)
  }
  assertFlag(renormalize)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  # Diversity functions
  div.func <- .DiversityFunction(entropy.use, order, top.n, renormalize)
  entropy.thresh <- .DiversityThreshold(entropy.use, entropy.thresh, order)
  .DiversityThresholdMessage(entropy.use, order, entropy.thresh, quiet)

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

  # Diversity
  object_grp_list <- .BuildGroupComparisons(object, group.1, group.2, unique_groups)
  comparison_mode <- attr(object_grp_list, "mode")
  if (!quiet && comparison_mode == "all") {
    message("Running DIV analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")
  } else if (!quiet && comparison_mode == "one_vs_all") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DIV analysis for ", comparison$grp1.names, " vs all other cells...")
  } else if (!quiet && comparison_mode == "pair") {
    comparison <- object_grp_list[["single_test"]]
    message("Running DIV analysis for ", comparison$grp1.names, " vs ", comparison$grp2.names, "...")
  }

  # Loop through object grp list
  data_list <- list()
  stats_list <- list()

  for (comp in names(object_grp_list)) {

    if (!quiet && comp != "single_test" && length(unique_groups) > 2) message(comp, ":")

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
    gene_dr_df <- gene_dr_df %>%
      filter(gene.pct.grp1 >= min.gene.pct &
               gene.pct.grp2 >= min.gene.pct &
               gene.sum.grp1 >= min.gene.cts &
               gene.sum.grp2 >= min.gene.cts)
    gene_dr_df$gene.id <- rownames(gene_dr_df)

    ## filter genes from grp objects
    filt_object_grp1 <- object_grp1[rowData(object_grp1)[[active.gene.id]] %in% rownames(gene_dr_df), , drop = FALSE]
    filt_object_grp2 <- object_grp2[rowData(object_grp2)[[active.gene.id]] %in% rownames(gene_dr_df), , drop = FALSE]

    ## aggregate transcript counts
    agg_cts_df <- data.frame("gene.id.1" = rowData(filt_object_grp1)[[active.gene.id]],
                             "gene.id.2" = rowData(filt_object_grp2)[[active.gene.id]],
                             "gene.id" = rowData(filt_object_grp1)[[active.gene.id]],
                             "cts.1" = rowSums(assay(filt_object_grp1, assay.use)),
                             "cts.2" = rowSums(assay(filt_object_grp2, assay.use)))
    agg_cts_df <- agg_cts_df %>%
      select(-gene.id.1, -gene.id.2) %>%
      rownames_to_column(var = "transcript") %>%
      left_join(., gene_dr_df[, c("gene.pct.grp1", "gene.pct.grp2", "gene.id")], by = "gene.id")
    unfiltered_agg_cts_df <- agg_cts_df

    ## filter transcripts
    agg_cts_df <- agg_cts_df %>%
      filter(cts.1 >= min.tx.cts | cts.2 >= min.tx.cts)

    ## remove genes that have <2 isoforms
    if (include.single == FALSE) {
      agg_cts_df <- agg_cts_df %>%
        group_by(gene.id) %>%
        filter(n() > 1) %>%
        ungroup()
    }
    keep_genes <- unique(agg_cts_df$gene.id)
    keep_transcripts <- agg_cts_df$transcript
    n_transcripts_df <- agg_cts_df %>%
      add_count(gene.id, name = "n.transcripts") %>%
      distinct(gene.id, n.transcripts)

    ## effective isoforms are determined independently within each group
    effective_grp1 <- unfiltered_agg_cts_df %>%
      filter(gene.id %in% keep_genes, cts.1 >= min.tx.cts) %>%
      group_by(gene.id) %>%
      summarise(
        n.effective.1 = .EffectiveIsoformCount(cts.1 / sum(cts.1), prop.thresh),
        .groups = "drop"
      )
    effective_grp2 <- unfiltered_agg_cts_df %>%
      filter(gene.id %in% keep_genes, cts.2 >= min.tx.cts) %>%
      group_by(gene.id) %>%
      summarise(
        n.effective.2 = .EffectiveIsoformCount(cts.2 / sum(cts.2), prop.thresh),
        .groups = "drop"
      )
    effective_data <- data.frame(gene.id = keep_genes) %>%
      left_join(effective_grp1, by = "gene.id") %>%
      left_join(effective_grp2, by = "gene.id") %>%
      mutate(
        n.effective.1 = coalesce(n.effective.1, 0L),
        n.effective.2 = coalesce(n.effective.2, 0L)
      )

    ## final filtering
    filt_object_grp1 <- object_grp1[keep_transcripts, , drop = FALSE]
    filt_object_grp2 <- object_grp2[keep_transcripts, , drop = FALSE]

    ## number of tests to conduct
    n_tests <- length(keep_genes)
    if (!quiet) message("  ", n_tests, " genes passed filtering.")
    if (n_tests == 0) {
      next
    }

    ## cell-to-cell diversity dispersion
    if (cell.dispersion) {
      if (!quiet) message("  Calculating cell-to-cell dispersion...")
      dispersion_gene_ids <- rowData(filt_object_grp1)[[active.gene.id]]
      dispersion_data_grp1 <- .CellDiversityDispersion(
        assay(filt_object_grp1, assay.use), dispersion_gene_ids, div.func
      )
      dispersion_data_grp2 <- .CellDiversityDispersion(
        assay(filt_object_grp2, assay.use), dispersion_gene_ids, div.func
      )
      names(dispersion_data_grp1)[-1] <- paste0(names(dispersion_data_grp1)[-1], ".1")
      names(dispersion_data_grp2)[-1] <- paste0(names(dispersion_data_grp2)[-1], ".2")
      dispersion_data <- left_join(
        dispersion_data_grp1, dispersion_data_grp2, by = "gene.id"
      )
    } else {
      dispersion_data <- data.frame(
        gene.id = keep_genes,
        cell.n.1 = rep(NA_integer_, length(keep_genes)),
        cell.div.mean.1 = rep(NA_real_, length(keep_genes)),
        cell.div.median.1 = rep(NA_real_, length(keep_genes)),
        cell.div.sd.1 = rep(NA_real_, length(keep_genes)),
        cell.div.iqr.1 = rep(NA_real_, length(keep_genes)),
        cell.n.2 = rep(NA_integer_, length(keep_genes)),
        cell.div.mean.2 = rep(NA_real_, length(keep_genes)),
        cell.div.median.2 = rep(NA_real_, length(keep_genes)),
        cell.div.sd.2 = rep(NA_real_, length(keep_genes)),
        cell.div.iqr.2 = rep(NA_real_, length(keep_genes))
      )
    }

    # Bootstrap comparisons
    if (!quiet) message("  Performing DIV comparisons...")

    if (is.null(boot.ncells)) {
      grp1_boot_ncells <- ceiling(ncol(filt_object_grp1) * boot.fraction)
      grp2_boot_ncells <- ceiling(ncol(filt_object_grp2) * boot.fraction)
    } else {
      grp1_boot_ncells <- boot.ncells
      grp2_boot_ncells <- boot.ncells
    }

    ## Work on count matrices directly; avoid copying SCE metadata per draw.
    bootstrap_matrices <- .BootstrapDiversityMatrices(
      assay(filt_object_grp1, assay.use), assay(filt_object_grp2, assay.use),
      rowData(filt_object_grp1)[[active.gene.id]],
      grp1_boot_ncells, grp2_boot_ncells, boot.iter, div.func
    )
    mat.div.1 <- bootstrap_matrices$div.1
    mat.div.2 <- bootstrap_matrices$div.2
    bootstrap_meta <- data.frame(
      gene.id = bootstrap_matrices$gene, grp.1 = group.1, grp.2 = group.2
    )

    ## comparison results
    comp_result <- data.frame(
      "gene.id" = bootstrap_meta$gene.id,
      "grp.1" = bootstrap_meta$grp.1,
      "grp.2" = bootstrap_meta$grp.2,
      "avgDiv.1" = rowMeans(mat.div.1, na.rm = TRUE),
      "avgDiv.2" = rowMeans(mat.div.2, na.rm = TRUE))

    bootstrap_stats <- bind_rows(lapply(seq_len(nrow(mat.div.1)), function(i) {
      .BootstrapDiversitySupport(
        mat.div.1[i, ], mat.div.2[i, ], div.diff.thresh, support.thresh
      )
    }))
    comp_result <- bind_cols(comp_result, bootstrap_stats)
    comp_result$avgDiv.1[!is.finite(comp_result$avgDiv.1)] <- NA_real_
    comp_result$avgDiv.2[!is.finite(comp_result$avgDiv.2)] <- NA_real_
    if (!quiet && any(comp_result$boot.valid.iter < boot.iter)) {
      message("  ", sum(comp_result$boot.valid.iter < boot.iter),
              " genes have unavailable bootstrap differences; see boot.valid.iter.")
    }

    comp_result <- comp_result %>%
      left_join(effective_data, by = "gene.id") %>%
      mutate(
        "class.1" = .DiversityClass(avgDiv.1, entropy.thresh),
        "class.2" = .DiversityClass(avgDiv.2, entropy.thresh)
      )

    ## stats and data output
    div_stats <- comp_result %>%
      rename("group.1" = grp.1,
             "group.2" = grp.2,
             "gene" = gene.id,
             "div.class.1" = class.1,
             "div.class.2" = class.2) %>%
      select(group.1, group.2, gene, avgDiv.1, avgDiv.2, div.diff,
             n.effective.1, n.effective.2, div.class.1, div.class.2,
             div.diff.lower, div.diff.upper, support.positive, support.negative,
             boot.valid.iter, supported)

    genes <- bootstrap_meta$gene.id
    grp1 <- bootstrap_meta$grp.1
    grp2 <- bootstrap_meta$grp.2
    div.grp1 <- lapply(seq_along(genes), function(i) mat.div.1[i, ])
    div.grp2 <- lapply(seq_along(genes), function(i) mat.div.2[i, ])

    div_data <- data.frame(
      "gene.id" = genes,
      "grp.1" = grp1,
      "grp.2" = grp2,
      "div.1" = I(div.grp1),
      "div.2" = I(div.grp2))

    div_data <- div_data %>%
      left_join(., gene_dr_df[, c("gene.pct.grp1", "gene.pct.grp2", "gene.id")], by = "gene.id") %>%
      left_join(., n_transcripts_df, by = "gene.id") %>%
      left_join(., dispersion_data, by = "gene.id") %>%
      rename("gene" = gene.id,
             "group.1" = grp.1,
             "group.2" = grp.2,
             "gene.pct.1" = gene.pct.grp1,
             "gene.pct.2" = gene.pct.grp2) %>%
      select(group.1, group.2, gene, gene.pct.1, gene.pct.2, n.transcripts,
             div.1, div.2, cell.n.1, cell.div.mean.1, cell.div.median.1,
             cell.div.sd.1, cell.div.iqr.1,
             cell.n.2, cell.div.mean.2, cell.div.median.2,
             cell.div.sd.2, cell.div.iqr.2)

    ## update list
    data_list[[comp]] <- div_data
    stats_list[[comp]] <- div_stats
  }

  # Output
  return_list <- list()
  if (length(data_list) > 0) {
    return_list$data <- as.data.frame(reduce(data_list, rbind)) %>%
      arrange(group.1, group.2, gene)
  } else {
    return_list$data <- data.frame("group.1" = character(),
                                   "group.2" = character(),
                                   "gene" = character(),
                                   "gene.pct.1" = numeric(),
                                   "gene.pct.2" = numeric(),
                                   "n.transcripts" = integer(),
                                   "div.1" = I(list()),
                                   "div.2" = I(list()),
                                   "cell.n.1" = integer(),
                                   "cell.div.mean.1" = numeric(),
                                   "cell.div.median.1" = numeric(),
                                   "cell.div.sd.1" = numeric(),
                                   "cell.div.iqr.1" = numeric(),
                                   "cell.n.2" = integer(),
                                   "cell.div.mean.2" = numeric(),
                                   "cell.div.median.2" = numeric(),
                                   "cell.div.sd.2" = numeric(),
                                   "cell.div.iqr.2" = numeric())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind)) %>%
      arrange(group.1, group.2, gene)

  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "avgDiv.1" = numeric(),
                                    "avgDiv.2" = numeric(),
                                    "div.diff" = numeric(),
                                    "n.effective.1" = integer(),
                                    "n.effective.2" = integer(),
                                    "div.class.1" = character(),
                                    "div.class.2" = character(),
                                    "div.diff.lower" = numeric(),
                                    "div.diff.upper" = numeric(),
                                    "support.positive" = numeric(),
                                    "support.negative" = numeric(),
                                    "boot.valid.iter" = integer(),
                                    "supported" = logical())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("0 genes passed detection thresholds (check min. parameters).")
  }

  if (!quiet) message("Done.")
  return(return_list)

}

# Summarize differences from matching bootstrap iterations. Missing values
# must be removed jointly so that draws from different iterations never pair.
.BootstrapDiversitySupport <- function(div.1, div.2, div.diff.thresh,
                                       support.thresh) {
  delta <- div.1 - div.2
  delta <- delta[is.finite(div.1) & is.finite(div.2) & is.finite(delta)]
  n_valid <- length(delta)
  result <- data.frame(
    div.diff = NA_real_,
    div.diff.lower = NA_real_,
    div.diff.upper = NA_real_,
    support.positive = NA_real_,
    support.negative = NA_real_,
    boot.valid.iter = as.integer(n_valid),
    supported = NA
  )
  if (n_valid == 0L) return(result)

  result$div.diff <- mean(delta)
  result$support.positive <- mean(delta > div.diff.thresh)
  result$support.negative <- mean(delta < -div.diff.thresh)
  if (n_valid >= 2L) {
    interval <- stats::quantile(delta, probs = c(0.025, 0.975), names = FALSE)
    result$div.diff.lower <- interval[1]
    result$div.diff.upper <- interval[2]
    result$supported <- result$support.positive >= support.thresh ||
      result$support.negative >= support.thresh
  }
  result
}

# Retain whole-cell transcript vectors and share each cell resample across genes.
.BootstrapDiversityMatrices <- function(counts.1, counts.2, gene.ids,
                                        ncells.1, ncells.2, boot.iter, div.func) {
  gene_rows <- split(seq_along(gene.ids), factor(gene.ids, levels = unique(gene.ids)))
  div.1 <- matrix(NA_real_, nrow = length(gene_rows), ncol = boot.iter)
  div.2 <- matrix(NA_real_, nrow = length(gene_rows), ncol = boot.iter)
  for (b in seq_len(boot.iter)) {
    cells.1 <- sample(seq_len(ncol(counts.1)), size = ncells.1, replace = TRUE)
    cells.2 <- sample(seq_len(ncol(counts.2)), size = ncells.2, replace = TRUE)
    cts.1 <- rowSums(counts.1[, cells.1, drop = FALSE])
    cts.2 <- rowSums(counts.2[, cells.2, drop = FALSE])
    div.1[, b] <- vapply(gene_rows, function(rows) {
      div.func(cts.1[rows] / sum(cts.1[rows]))
    }, numeric(1))
    div.2[, b] <- vapply(gene_rows, function(rows) {
      div.func(cts.2[rows] / sum(cts.2[rows]))
    }, numeric(1))
  }
  list(gene = names(gene_rows), div.1 = div.1, div.2 = div.2)
}
