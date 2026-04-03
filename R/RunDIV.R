#' Run isoform diversity analysis
#'
#' Performs isoform diversity analysis.
#' For each comparison, genes are first filtered according to coverage.
#' Gene-level diversity is then calculated by bootstrapping cells in each group. For each bootstrap iteration, transcript counts are aggregated to compute diversity metrics.
#' These diversities are then compared across cell groups using the Wilcox test.
#' @param object A `SingleCellExperiment` object.
#' @param group.by Name of `colData` variable to group cells for comparisons. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.1 First group in the comparison.
#' @param group.2 Second group in the comparison.
#' @param entropy.use The diversity index to calculate.
#' Options include `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `InverseSimpson`. The default is `"Tsallis"`.
#' @param assay.use Which `assay` (counts) to use.
#' @param entropy.thresh Diversity index threshold to use for monoform and polyform classifications.
#' Default cutoffs are 0.243 for Tsallis, 0.500 for Shannon, 0 for normalized Shannon, 0.435 for Renyi, 0 for normalized Renyi, 0.348 for Gini-Simpson, and 1.533 for inverse Simpson.
#' @param min.gene.pct Minimum percentage of cells in which a gene must be expressed in both groups for it to be tested.
#' @param min.gene.cts Minimum total transcript counts in which a gene must be have in both groups for it be tested.
#' @param min.tx.cts Minimum transcript counts required for a transcript to be included in the contingency table.
#' @param boot.iter Number of bootstrap iterations to perform.
#' @param boot.size Proportion of cells per group from 0.01 to 1 to use for bootstrap subsampling.
#' @param include.single Logical; if `FALSE`, genes with only one associated transcript after filtering will be excluded from the analysis.
#' @param order Value specifying the order of entropy. Corresponds to `q` for Tsallis and `alpha` for Renyi.
#' @param genes A vector containing one or more gene IDs to test.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A list containing two data frames:
#' \describe{
#'   \item{`$data`}{
#'     A data frame of summarized data with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`div.1`}{Bootstrapped isoform diversity values of each gene for `group.1`.}
#'       \item{`div.2`}{Bootstrapped isoform diversity values of each gene for `group.2`.}
#'     }
#'     Includes NA entries derived from bootstrapped sampling that resulted in 0 gene counts.
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`gene`}{The gene being tested.}
#'       \item{`gene.pct.1`}{Percentage of cells in `group.1` with expression of the gene.}
#'       \item{`gene.pct.2`}{Percentage of cells in `group.2` with expression of the gene.}
#'       \item{`n.transcripts`}{Number of transcripts associated with the gene.}
#'       \item{`avgDiv.1`}{Average of bootstrapped isoform diversity of the gene for cells in `group.1`.}
#'       \item{`avgDiv.2`}{Average of bootstrapped isoform diversity of the gene for cells in `group.2`.}
#'       \item{`delta.div`}{The difference in averaged isoform diversities between groups (`group.1` - `group.2`).}
#'       \item{`log2FC`}{The log2 fold-change of averaged isoform diversities between groups (`group.1` - `group.2`).}
#'       \item{`pval`}{P-value from the Wilcox test comparing bootstrapped isoform diversities across groups.}
#'       \item{`padj`}{Adjusted p-value (Bonferroni correction).}
#'       \item{`div.class.1`}{Isoform diversity classification of the gene for `group.1`.}
#'       \item{`div.class.2`}{Isoform diversity classification of the gene for `group.2`.}
#'     }
#'   }
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom purrr reduce map transpose
#' @importFrom stats wilcox.test

RunDIV <- function (
    object,
    group.by = NULL,
    group.1 = NULL,
    group.2 = NULL,
    entropy.use = "Tsallis",
    assay.use = "counts",
    entropy.thresh = NULL,
    min.gene.pct = 0.05,
    min.gene.cts = 15,
    min.tx.cts = 1,
    boot.iter = 100,
    boot.size = 0.6,
    include.single = TRUE,
    order = NULL,
    genes = NULL,
    quiet = FALSE
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
  assertChoice(entropy.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(min.gene.pct, lower = 0, upper = 1, finite = TRUE)
  assertNumber(min.gene.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(boot.iter, lower = 1, finite = TRUE)
  assertNumeric(boot.size, lower = 0.01, upper = 1, finite = TRUE)
  assertFlag(include.single)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertTRUE(order != 1 || is.null(order))
  assertCharacter(genes, null.ok = TRUE, any.missing = FALSE, unique = TRUE)
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

  # Diversity functions
  div.func <- function(x) {

    if (entropy.use == "Shannon") {
      x <- head(sort(x, decreasing = TRUE), 2)
      -sum(x[x > 0] * log(x[x > 0]))
    }
    else if (entropy.use == "NormalizedShannon") {
      n_x <- sum(x > 0)
      (-sum(x[x > 0] * log(x[x > 0]))) / (log(n_x))
    }
    else if (entropy.use == "Renyi") {
      if (is.null(order)) {order <- 2}
      x <- head(sort(x, decreasing = TRUE), 2)
      (1 / (1 - order)) * log( sum( (x[x > 0])^order ) )
    }
    else if (entropy.use == "NormalizedRenyi") {
      if (is.null(order)) {order <- 2}
      n_x <- sum(x > 0)
      (1 / (1 - order)) * log( sum( (x[x > 0])^order ) ) / (log(n_x))
    }
    else if (entropy.use == "GiniSimpson") {
      1 - sum( (x[x > 0])^2 )
    }
    else if (entropy.use == "Tsallis") {
      if (is.null(order)) {order <- 3}
      (1 - sum(x[x > 0]^order)) / (order - 1)
    }
    else if (entropy.use == "InverseSimpson") {
      1 / sum( (x[x > 0])^2 )
    }
  }

  # Diversity thresholds
  if (is.null(entropy.thresh)) {
    if (entropy.use == "Shannon") {entropy.thresh <- 0.500}
    else if (entropy.use == "NormalizedShannon") {entropy.thresh <- 0}
    else if (entropy.use == "Renyi") {entropy.thresh <- 0.435}
    else if (entropy.use == "NormalizedRenyi") {entropy.thresh <- 0}
    else if (entropy.use == "GiniSimpson") {entropy.thresh <- 0.348}
    else if (entropy.use == "Tsallis") {entropy.thresh <- 0.243}
    else if (entropy.use == "InverseSimpson") {entropy.thresh <- 1.533}
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

  # Diversity
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

    if (!quiet) message("Running DIV analysis for all groups in '", paste0(group.by, collapse = "_"), "'...")

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

    if (!quiet) message("Running DIV analysis for ", group.1.updated, " vs all other cells...")
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

    if (!quiet) message("Running DIV analysis for ", group.1.updated, " vs ", group.2.updated, "...")

  }

  ## case: group 1 unspecified but group 2 is specified
  else if (is.null(group.1) && !is.null(group.2)) {
    stop("`group.1` must be specified prior to `group.2`")
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
    filt_object_grp1 <- object_grp1[rowData(object_grp1)[[active.gene.id]] %in% rownames(gene_dr_df), drop = FALSE]
    filt_object_grp2 <- object_grp2[rowData(object_grp2)[[active.gene.id]] %in% rownames(gene_dr_df), drop = FALSE]

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
    n_transcripts_df <- agg_cts_df %>%
      add_count(gene.id, name = "n.transcripts") %>%
      distinct(gene.id, n.transcripts)

    ## final filtering
    filt_object_grp1 <- object_grp1[rowData(object_grp1)[[active.gene.id]] %in% keep_genes, drop = FALSE]
    filt_object_grp2 <- object_grp2[rowData(object_grp2)[[active.gene.id]] %in% keep_genes, drop = FALSE]

    ## number of tests to conduct
    n_tests <- length(keep_genes)
    if (!quiet) message("  ", n_tests, " genes passed detection thresholds.")
    if (n_tests == 0) {
      next
    }

    # Bootstrap comparisons
    if (!quiet) message("  Performing DIV comparisons...")

    grp1_boot_ncells <- ceiling(ncol(filt_object_grp1) * boot.size)
    grp2_boot_ncells <- ceiling(ncol(filt_object_grp2) * boot.size)

    ## perform bootstrap subsampling on filtered object
    boot_result <- replicate(boot.iter, {
      ## col index
      grp1_col_idx <- sample(seq_len(ncol(filt_object_grp1)), size = grp1_boot_ncells, replace = TRUE)
      grp2_col_idx <- sample(seq_len(ncol(filt_object_grp2)), size = grp2_boot_ncells, replace = TRUE)

      ## subsample objects
      boot_object_grp1 <- filt_object_grp1[, grp1_col_idx]
      boot_object_grp2 <- filt_object_grp2[, grp2_col_idx]

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

    })

    ## comparison results
    mat.div.1 <- do.call(cbind, boot_result["div.1", ])
    mat.div.2 <- do.call(cbind, boot_result["div.2", ])

    comp_result <- data.frame(
      "gene.id" = boot_result["gene.id", ][[1]],
      "grp.1" = boot_result["grp.1", ][[1]],
      "grp.2" = boot_result["grp.2", ][[1]],
      "avgDiv.1" = rowMeans(mat.div.1, na.rm = TRUE),
      "avgDiv.2" = rowMeans(mat.div.2, na.rm = TRUE))

    comp_result <- comp_result %>%
      mutate(
      "delta.div" = avgDiv.1 - avgDiv.2,
      "log2FC" = log2(avgDiv.1 / avgDiv.2))

    if (sum(is.na(comp_result$delta.div)) > 0) {
      if (!quiet) message("  \u2139 Warning: ", sum(is.na(comp_result$delta.div)), " genes have 0 counts after bootstrap sampling.\n    Try increasing filtering parameters or iterations.")
    }

    ## Wilcox test
    pvals <- sapply(1:nrow(mat.div.1), function(i) {
      tryCatch(
        {wilcox.test(mat.div.1[i, ], mat.div.2[i, ], paired = FALSE, exact = FALSE)$p.value},
        error = function(e) {NA})

      })
    comp_result <- comp_result %>%
      mutate("pval" = pvals)

    ## remove failed bootstrap genes
    comp_result <- comp_result %>%
      filter(!is.na(delta.div))

    ## p value adjustment
    comp_result <- comp_result %>%
      mutate("padj" = p.adjust(pval, method = "bonferroni"),
             "class.1" = ifelse(avgDiv.1 <= entropy.thresh, "monoform", "polyform"),
             "class.2" = ifelse(avgDiv.2 <= entropy.thresh, "monoform", "polyform")
             )

    ## stats and data output
    div_stats <- comp_result %>%
      left_join(., gene_dr_df[, c("gene.pct.grp1", "gene.pct.grp2", "gene.id")], by = "gene.id") %>%
      left_join(., n_transcripts_df, by = "gene.id") %>%
      rename("group.1" = grp.1,
             "group.2" = grp.2,
             "gene" = "gene.id",
             "gene.pct.1" = "gene.pct.grp1",
             "gene.pct.2" = "gene.pct.grp2",
             "div.class.1" = class.1,
             "div.class.2" = class.2) %>%
      select(group.1, group.2, gene, gene.pct.1, gene.pct.2, n.transcripts, avgDiv.1, avgDiv.2, delta.div, log2FC, everything())
    tmp <- boot_result
    boot_result <- boot_result %>%
      t() %>%
      as.data.frame() %>%
      mutate(grp.1 = unique(unlist(grp.1)),
             grp.2 = unique(unlist(grp.2)))

    genes <- boot_result$gene.id[[1]]
    grp1 <- rep(boot_result$grp.1[1], length(genes))
    grp2 <- rep(boot_result$grp.2[1], length(genes))
    div.grp1 <- map(transpose(boot_result$div.1), unlist)
    div.grp2 <- map(transpose(boot_result$div.2), unlist)

    div_data <- data.frame(
      "gene.id" = genes,
      "grp.1" = grp1,
      "grp.2" = grp2,
      "div.1" = I(div.grp1),
      "div.2" = I(div.grp2)) %>%
      select(grp.1, grp.2, gene.id, div.1, div.2) %>%
      rename("gene" = "gene.id",
             "group.1" = grp.1,
             "group.2" = grp.2)

    ## update list
    data_list[[comp]] <- div_data
    stats_list[[comp]] <- div_stats
  }

  # Output
  return_list <- list()
  if (length(data_list) > 0) {
    return_list$data <- as.data.frame(reduce(data_list, rbind)) %>%
      arrange(group.1)
  } else {
    return_list$data <- data.frame("group.1" = character(),
                                   "group.2" = character(),
                                   "gene" = character(),
                                   "div.1" = numeric(),
                                   "div.2" = numeric())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind)) %>%
      arrange(padj)

  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "gene" = character(),
                                    "gene.pct.1" = numeric(),
                                    "gene.pct.2" = numeric(),
                                    "n.transcripts" = numeric(),
                                    "avgDiv.1" = numeric(),
                                    "avgDiv.2" = numeric(),
                                    "log2FC" = numeric(),
                                    "delta.div" = character(),
                                    "pval" = numeric(),
                                    "padj" = numeric(),
                                    "div.class.1" = character(),
                                    "div.class.2" = character())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("0 genes passed detection thresholds (check min. parameters).")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
