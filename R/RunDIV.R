#' Run isoform diversity analysis
#'
#' Runs isoform diversity analysis.
#' For each comparison, genes are first filtered according to coverage, counts are aggregated across groups, and isoform diversity is calculated.
#' The isoform diversity of genes are compared using pairwise Wilcoxon test and heterogeneity classification are compared using the pairwise McNemar test.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by Name of `colData` variable to group cells for comparisons. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.1 First group in the comparison.
#' @param group.2 Second group in the comparison.
#' @param method.use The diversity index to calculate.
#' Options include `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `InverseSimpson`. The default is `"Tsallis"`.
#' @param assay.use Which `assay` (counts) to use.
#' @param diversity.cutoff The cutoff of the diversity index used for monoform and polyform classification.
#' Default cutoffs are 0.243 for Tsallis, 0.500 for Shannon, 0 for normalized Shannon, 0.435 for Renyi, 0 for normalized Renyi, 0.348 for Gini-Simpson, and 1.533 for inverse Simpson.
#' @param min.gene.pct Minimum percentage of cells in which a gene must be expressed in both groups for it to be tested.
#' @param min.gene.cts Minimum total transcript counts in which a gene must be have in both groups for it be tested.
#' @param min.tx.cts Minimum transcript counts required for a transcript to be included in the contingency table.
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
#'       \item{`gene.pct.1`}{Percentage of cells in `group.1` with expression of the gene.}
#'       \item{`gene.pct.2`}{Percentage of cells in `group.2` with expression of the gene.}
#'       \item{`n.transcripts`}{Number of transcripts associated with the gene.}
#'       \item{`div.1`}{Isoform diversity of the gene for cells in `group.1`.}
#'       \item{`div.2`}{Isoform diversity of the gene for cells in `group.2`.}
#'       \item{`div.diff`}{The difference in isoform diversity between groups (`group.1` - `group.2`).}
#'       \item{`div.class.1`}{Isoform diversity classification of the gene for `group.1`.}
#'       \item{`div.class.2`}{Isoform diversity classification of the gene for `group.2`.}
#'       \item{`div.class.diff`}{A logical indicating difference between isoform diversity classifications between `group.1` and `group.2.`}
#'     }
#'   }
#'
#'   \item{`$stats`}{
#'     A data frame containing statistical results with columns:
#'     \describe{
#'       \item{`group.1` & `group.2`}{The two cell groups being compared.}
#'       \item{`avgDiv.1`}{Average isoform diversity of genes from cells in `group.1`.}
#'       \item{`avgDiv.2`}{Average isoform diversity of genes from cells in `group.2`.}
#'       \item{`wilcox.pval`}{P-value from the paired Wilcoxon rank sum test comparing isoform diversity across the two groups.}
#'       \item{`mono.poly.1`}{Number of genes classified monoform to polyform from cells in `group.1`.}
#'       \item{`mono.poly.2`}{Number of genes classified monoform to polyform from cells in `group.1`.}
#'       \item{`mcnemar.pval`}{P-value from the pairwise McNemar test comparing isoform diversity classifications across the two groups.}
#'     }
#'   }
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom purrr reduce

RunDIV <- function (
    object,
    group.by = NULL,
    group.1 = NULL,
    group.2 = NULL,
    method.use = "Tsallis",
    assay.use = "counts",
    diversity.cutoff = NULL,
    min.gene.pct = 0.05,
    min.gene.cts = 15,
    min.tx.cts = 1,
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
  assertChoice(method.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(min.gene.pct, lower = 0, upper = 1, finite = TRUE)
  assertNumber(min.gene.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
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
    if (include.single == FALSE) {
      agg_cts_df <- agg_cts_df %>%
        group_by(gene.id) %>%
        filter(n() > 1) %>%
        ungroup()
    }

    ## number of tests to conduct
    n_tests <- length(unique(agg_cts_df$gene.id))
    if (!quiet) message(n_tests, " genes passed detection thresholds.")
    if (n_tests == 0) {
      next
    }

    ## calculate diversity
    if (!quiet) message("Calculating isoform diversity...")

    div.func <- function(x) {

      if (method.use == "Shannon") {
        x <- head(sort(x, decreasing = TRUE), 2)
        -sum(x[x > 0] * log(x[x > 0]))
      }
      else if (method.use == "NormalizedShannon") {
        n_x <- sum(x > 0)
        (-sum(x[x > 0] * log(x[x > 0]))) / (log(n_x))
      }
      else if (method.use == "Renyi") {
        if (is.null(order)) {order <- 2}
        x <- head(sort(x, decreasing = TRUE), 2)
        (1 / (1 - order)) * log( sum( (x[x > 0])^order ) )
      }
      else if (method.use == "NormalizedRenyi") {
        if (is.null(order)) {order <- 2}
        n_x <- sum(x > 0)
        (1 / (1 - order)) * log( sum( (x[x > 0])^order ) ) / (log(n_x))
      }
      else if (method.use == "GiniSimpson") {
        1 - sum( (x[x > 0])^2 )
      }
      else if (method.use == "Tsallis") {
        if (is.null(order)) {order <- 3}
        (1 - sum(x[x > 0]^order)) / (order - 1)
      }
      else if (method.use == "InverseSimpson") {
        1 / sum( (x[x > 0])^2 )
      }
    }

    ## diversity threshold
    if (is.null(diversity.cutoff)) {
      if (method.use == "Shannon") {diversity.cutoff <- 0.500}
      else if (method.use == "NormalizedShannon") {diversity.cutoff <- 0}
      else if (method.use == "Renyi") {diversity.cutoff <- 0.435}
      else if (method.use == "NormalizedRenyi") {diversity.cutoff <- 0}
      else if (method.use == "GiniSimpson") {diversity.cutoff <- 0.348}
      else if (method.use == "Tsallis") {diversity.cutoff <- 0.243}
      else if (method.use == "InverseSimpson") {diversity.cutoff <- 1.533}
    }

    div_data <- agg_cts_df %>%
      group_by(gene.id) %>%
      mutate(grp.1 = group.1,
             grp.2 = group.2,
             n.transcripts = n(),
             prop.1 = cts.1 / sum(cts.1),
             prop.2 = cts.2 / sum(cts.2),
             div.1 = div.func(x = prop.1),
             div.2 = div.func(x = prop.2)) %>%
      ungroup() %>%
      mutate(ddiv = div.1 - div.2,
             class.1 = ifelse(div.1 <= diversity.cutoff, "monoform", "polyform"),
             class.2 = ifelse(div.2 <= diversity.cutoff, "monoform", "polyform"),
             dclass = class.1 != class.2) %>%
      distinct(grp.1, grp.2, gene.id, gene.pct.grp1, gene.pct.grp2, n.transcripts,
               div.1, div.2, ddiv, class.1, class.2, dclass)

    ## pairwise wilcox test
    if(!all(dim(table(div_data$class.1, div_data$class.2)) == c(2, 2))) {
      if (!quiet) message("\u2139 Warning: Not enough monoform/polyform classifications for tests.")
      next
    }
    mean.div.1 <- mean(div_data$div.1)
    mean.div.2 <- mean(div_data$div.2)
    n.mono.1 <- sum(div_data$class.1 == "monoform")
    n.poly.1 <- sum(div_data$class.1 == "polyform")
    n.mono.2 <- sum(div_data$class.2 == "monoform")
    n.poly.2 <- sum(div_data$class.2 == "polyform")
    wilcox.p <- wilcox.test(div_data$div.1, div_data$div.2, paired = TRUE, exact = FALSE)$p.value
    mcnemar.p <- mcnemar.test(table(div_data$class.1, div_data$class.2))$p.value

    div_data <- div_data %>%
      rename("group.1" = grp.1,
             "group.2" = grp.2,
             "gene" = "gene.id",
             "gene.pct.1" = "gene.pct.grp1",
             "gene.pct.2" = "gene.pct.grp2",
             "div.1" = div.1,
             "div.2" = div.2,
             "div.class.1" = class.1,
             "div.class.2" = class.2,
             "div.diff" = ddiv,
             "div.class.diff" = dclass)

    div_stats <- data.frame("group.1" = group.1,
                            "group.2" = group.2,
                            "avgDiv.1" = mean.div.1,
                            "avgDiv.2" = mean.div.2,
                            "wilcox.pval" = wilcox.p,
                            "mono.poly.1" = paste0(n.mono.1, "/", n.poly.1),
                            "mono.poly.2" = paste0(n.mono.2, "/", n.poly.2),
                            "mcnemar.pval" = mcnemar.p)

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
                                   "n.transcripts" = numeric(),
                                   "div.1" = numeric(),
                                   "div.2" = numeric(),
                                   "div.diff" = numeric(),
                                   "div.class.1" = character(),
                                   "div.class.2" = character(),
                                   "div.class.diff" = logical())
  }
  if (length(stats_list) > 0) {
    return_list$stats <- as.data.frame(reduce(stats_list, rbind)) %>%
      arrange(mcnemar.pval)

  } else {
    return_list$stats <- data.frame("group.1" = character(),
                                    "group.2" = character(),
                                    "avgDiv.1" = numeric(),
                                    "avgDiv.2" = numeric(),
                                    "wilcox.pval" = numeric(),
                                    "mono.poly.1" = character(),
                                    "mono.poly.2" = character(),
                                    "mcnemar.pval" = numeric())
  }

  if (length(stats_list) == 0 && length(data_list) == 0) {
    stop("Either 0 genes passed detection thresholds (check min. parameters) \nor there were insufficient polyform/monoform classifications for comparisons (check diversity.cutoff).")
  }

  if (!quiet) message("Done.")
  return(return_list)

}
