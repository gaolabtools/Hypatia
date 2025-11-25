#' Get summarized isoform diversity data
#'
#' Retrieves summarized isoform diversity data of one or more genes.
#'
#' @param object A `SingleCellExperiment` object.
#' @param genes A vector of one or more gene IDs.
#' @param group.by Name of `colData` variable to group cells. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.subset An optional vector specifying a subset of the elements in `group.by` to include.
#' @param method.use The diversity index to calculate.
#' Options include `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `InverseSimpson`.
#' @param assay.use Which `assay` (counts) to use.
#' @param diversity.cutoff The cutoff of the diversity index used for monoform and polyform classification.
#' Default cutoffs are 0.243 for Tsallis, 0.500 for Shannon, 0 for normalized Shannon, 0.435 for Renyi, 0 for normalized Renyi, 0.348 for Gini-Simpson, and 1.533 for inverse Simpson.
#' @param min.tx.cts Minimum transcript counts required for a transcript to be included in the contingency table.
#' @param order Value specifying the order of entropy. Corresponds to `q` for Tsallis (default: 3) and `alpha` for Renyi (default: 2).
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The gene being queried}
#'   \item{`gene.pct`}{Percentage of cells in `group` with expression of the gene.}
#'   \item{`n.transcripts`}{Number of associated transcripts for the gene.}
#'   \item{`transcript`}{The associated transcript.}
#'   \item{`cts`}{Total counts of the transcript in `group`.}
#'   \item{`prop`}{The transcript proportion in `group`.}
#'   \item{`div`}{Isoform diversity of the gene in `group`.}
#' }
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom purrr reduce
#' @importFrom tidyr unite
#' @importFrom Matrix rowSums

GetDiversity <- function (
    object,
    genes,
    group.by = NULL,
    group.subset = NULL,
    method.use = "Tsallis",
    assay.use = "counts",
    diversity.cutoff = NULL,
    min.tx.cts = 1,
    order = NULL,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertCharacter(genes, any.missing = FALSE, unique = TRUE)
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.subset, null.ok = TRUE)
  assertChoice(method.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(diversity.cutoff, null.ok = TRUE, finite = TRUE, lower = 0)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertTRUE(order != 1 || is.null(order))
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

  # Gene filter
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

  # Group structure
  colData(object)$group_var <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE) %>%
    pull(group_var)
  unique_groups <- unique(colData(object)$group_var)
  ## check groups
  if (!is.null(group.subset)) {
    assertSubset(group.subset, unique_groups)
    ## subset object for groups
    object <- object[, object$group_var %in% group.subset]
    unique_groups <- unique(colData(object)$group_var)
  }

  # Diversity
  ## loop through each group
  res_list <- list()
  for (group in unique_groups) {

    ## subset group
    object_grp <- object[, object$group_var == group]

    ## gene pct
    gene_groups <- rowData(object_grp)[[active.gene.id]]
    expr_mat_gene <- assay(object_grp, assay.use)
    expr_mat_gene <- rowsum(expr_mat_gene, group = gene_groups)
    gene_pct <- rowSums(expr_mat_gene > 0) / ncol(expr_mat_gene)
    gene_pct_df <- data.frame("gene.pct" = gene_pct) %>%
      rownames_to_column(var = "gene_query")

    ## aggregate transcript counts
    agg_cts_df <- data.frame("gene_query" = rowData(object_grp)[[active.gene.id]],
                             "cts" = rowSums(assay(object_grp, assay.use))) %>%
      rownames_to_column(var = "transcripts_query")
    agg_cts_df <- left_join(agg_cts_df, gene_pct_df, by = "gene_query")

    ## filter transcripts
    agg_cts_df <- agg_cts_df %>%
      filter(cts >= min.tx.cts) %>%
      mutate("group_var" = group)

    ## calculate diversity
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
        # x <- head(sort(x, decreasing = TRUE), 2)
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

    div_res <- agg_cts_df %>%
      group_by(gene_query) %>%
      mutate(n.transcripts = n(),
             prop = cts / sum(cts),
             diversity = div.func(x = prop),
             class = ifelse(diversity <= diversity.cutoff, "monoform", "polyform")
      ) %>%
      ungroup() %>%
      mutate(prop = ifelse(is.nan(prop), NA, prop),
             diversity = ifelse(is.na(prop), NA, diversity)) %>%
      select(group_var, gene_query, gene.pct, n.transcripts, transcripts_query,
             cts, prop, diversity, class) %>%
      rename("group" = "group_var",
             "transcript" = "transcripts_query",
             "gene" = "gene_query",
             "cts" = cts,
             "prop" = prop,
             "div" = diversity,
             "div.class" = class) %>%
      filter(gene.pct > 0)


    res_list[[group]] <- div_res
  }

  return(
    as.data.frame(reduce(res_list, rbind)) %>%
      arrange(group, gene)
  )

}



