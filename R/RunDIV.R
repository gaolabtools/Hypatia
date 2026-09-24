#' Run isoform diversity analysis
#'
#' Tests genes for differential isoform diversity between cell groups.
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
#' @param boot.iter Number of bootstrap iterations to perform.
#' @param boot.fraction Proportion of cells sampled from each comparison group per bootstrap iteration.
#' @param boot.ncells Optional fixed number of cell draws sampled with replacement from each comparison group per bootstrap iteration. Cannot be supplied together with `boot.fraction`.
#' @param include.single Logical; if `FALSE`, genes with only one associated transcript after filtering will be excluded from the analysis.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi. At order 1, Tsallis and Renyi use their Shannon entropy limit, and NormalizedRenyi uses normalized Shannon entropy.
#' @param genes Optional vector of active gene IDs to test. Genes are still subject to filtering.
#' @param p.adj P-value adjustment method. Must be one of `stats::p.adjust.methods`.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param top.n Optional number of the most abundant isoforms to include in diversity calculations. If `NULL`, all isoforms are included. Must be at least 2 when supplied.
#' @param renormalize Logical; if `TRUE`, rescale the selected isoform proportions to sum to one before calculating diversity.
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
#'       \item{`div.diff`}{The difference in averaged isoform diversities between groups (`group.1` - `group.2`).}
#'       \item{`pval`}{P-value from the Wilcoxon rank-sum test.}
#'       \item{`padj`}{Adjusted p-value, calculated separately within each comparison.}
#'       \item{`n.effective.1`, `n.effective.2`}{Number of transcripts with within-gene proportion greater than or equal to `prop.thresh` in each group.}
#'       \item{`div.class.1`}{`"monoform"` when `avgDiv.1` is at or below `entropy.thresh` and `"polyform"` otherwise.}
#'       \item{`div.class.2`}{`"monoform"` when `avgDiv.2` is at or below `entropy.thresh` and `"polyform"` otherwise.}
#'     }
#'   }
#' }
#' Rows in `$data` are ordered by `group.1`, `group.2`, and `gene`. Rows in
#' `$stats` are ordered by `group.1`, `group.2`, `padj`, and `gene`.
#' @details Genes must meet `min.gene.pct` and `min.gene.cts` in both comparison
#' groups. Transcripts are retained when their total counts meet `min.tx.cts`
#' in either group. This retained set is used for bootstrapping; `include.single`
#' controls whether genes with only one retained transcript are tested.
#'
#' Cells are sampled with replacement within each group for `boot.iter`
#' iterations. Each iteration draws `ceiling(boot.fraction * ncol(group))` cells,
#' or `boot.ncells` when supplied, and calculates diversity from their pooled
#' transcript counts. The entropy indices and optional `top.n` selection and
#' renormalization follow [GetDiversity()]. A two-sided, unpaired Wilcoxon
#' rank-sum test with an asymptotic p-value compares the bootstrap diversity
#' values. `avgDiv` averages available bootstrap values, and classification
#' applies `entropy.thresh` to that average.
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
#' P-values are adjusted within each comparison using Bonferroni by default.
#' Call `set.seed()` for reproducible results; no seed is set internally.
#' An error is returned if no genes pass filtering across all comparisons.
#' @seealso [GetDiversity()], [PlotDiversity()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom purrr reduce
#' @importFrom stats p.adjust wilcox.test

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
    boot.iter = 100,
    boot.fraction = 0.6,
    boot.ncells = NULL,
    include.single = TRUE,
    order = NULL,
    genes = NULL,
    p.adj = "bonferroni",
    quiet = FALSE,
    cell.dispersion = FALSE,
    top.n = NULL,
    renormalize = FALSE
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
  assertNumber(boot.iter, lower = 1, finite = TRUE)
  boot.fraction.supplied <- !missing(boot.fraction)
  assertNumber(boot.fraction, lower = 0.01, upper = 1, finite = TRUE)
  assertCount(boot.ncells, positive = TRUE, null.ok = TRUE)
  if (!is.null(boot.ncells) && boot.fraction.supplied) {
    stop("Please provide only one of `boot.fraction` or `boot.ncells`.", call. = FALSE)
  }
  assertFlag(include.single)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertCharacter(genes, null.ok = TRUE, any.missing = FALSE, unique = TRUE)
  p.adj <- .PAdjustMethod(p.adj)
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

    ## perform bootstrap subsampling on filtered object
    boot_result <- replicate(boot.iter, {
      ## col index
      grp1_col_idx <- sample(seq_len(ncol(filt_object_grp1)), size = grp1_boot_ncells, replace = TRUE)
      grp2_col_idx <- sample(seq_len(ncol(filt_object_grp2)), size = grp2_boot_ncells, replace = TRUE)

      ## subsample objects
      boot_object_grp1 <- filt_object_grp1[, grp1_col_idx, drop = FALSE]
      boot_object_grp2 <- filt_object_grp2[, grp2_col_idx, drop = FALSE]

      ## aggregate transcript counts
      boot_agg_cts_df <- data.frame("gene.id.1" = rowData(boot_object_grp1)[[active.gene.id]],
                               "gene.id.2" = rowData(boot_object_grp2)[[active.gene.id]],
                               "gene.id" = rowData(boot_object_grp1)[[active.gene.id]],
                               "cts.1" = rowSums(assay(boot_object_grp1, assay.use)),
                               "cts.2" = rowSums(assay(boot_object_grp2, assay.use)))

      if (!identical(boot_agg_cts_df$gene.id.1, boot_agg_cts_df$gene.id.2)) {
        stop("Gene ids are inconsistent.")
      }

      boot_agg_cts_df <- boot_agg_cts_df %>%
        select(-gene.id.1, -gene.id.2) %>%
        rownames_to_column(var = "transcript")

      ## calculate diversity
      div_data <- boot_agg_cts_df %>%
        group_by(gene.id) %>%
        mutate(grp.1 = group.1,
               grp.2 = group.2,
               prop.1 = cts.1 / sum(cts.1),
               prop.2 = cts.2 / sum(cts.2),
               div.1 = div.func(x = prop.1),
               div.2 = div.func(x = prop.2)) %>%
        ungroup() %>%
        distinct(gene.id, grp.1, grp.2, div.1, div.2)

    },
    simplify = FALSE
    )

    ## comparison results
    mat.div.1 <- do.call(cbind, lapply(boot_result, `[[`, "div.1"))
    mat.div.2 <- do.call(cbind, lapply(boot_result, `[[`, "div.2"))

    comp_result <- data.frame(
      "gene.id" = boot_result[[1]]$gene.id,
      "grp.1" = boot_result[[1]]$grp.1,
      "grp.2" = boot_result[[1]]$grp.2,
      "avgDiv.1" = rowMeans(mat.div.1, na.rm = TRUE),
      "avgDiv.2" = rowMeans(mat.div.2, na.rm = TRUE))

    comp_result <- comp_result %>%
      mutate(
        "div.diff" = avgDiv.1 - avgDiv.2,
      )

    if (sum(is.na(comp_result$div.diff)) > 0) {
      if (!quiet) message("  \u2139 Warning: ", sum(is.na(comp_result$div.diff)), " genes have 0 counts after bootstrap sampling.\n    Try increasing filtering parameters or iterations.")
    }

    ## Wilcox test
    pvals <- sapply(1:nrow(mat.div.1), function(i) {
      tryCatch(
        {wilcox.test(mat.div.1[i, ], mat.div.2[i, ], paired = FALSE, exact = FALSE)$p.value},
        error = function(e) {NA}
      )
    })
    comp_result <- comp_result %>%
      mutate("pval" = pvals)

    ## remove failed bootstrap genes
    comp_result <- comp_result %>%
      filter(!is.na(div.diff))

    ## p value adjustment
    comp_result <- comp_result %>%
      left_join(effective_data, by = "gene.id") %>%
      mutate(
        "padj" = p.adjust(pval, method = p.adj),
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
             pval, padj, n.effective.1, n.effective.2,
             div.class.1, div.class.2)

    genes <- boot_result[[1]]$gene.id
    grp1 <- boot_result[[1]]$grp.1
    grp2 <- boot_result[[1]]$grp.2
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
      arrange(group.1, group.2, padj, gene)

  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "avgDiv.1" = numeric(),
                                    "avgDiv.2" = numeric(),
                                    "div.diff" = numeric(),
                                    "pval" = numeric(),
                                    "padj" = numeric(),
                                    "n.effective.1" = integer(),
                                    "n.effective.2" = integer(),
                                    "div.class.1" = character(),
                                    "div.class.2" = character())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("0 genes passed detection thresholds (check min. parameters).")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
