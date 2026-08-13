#' Get isoform diversity summaries
#'
#' Summarizes transcript proportions and isoform diversity for one or more genes.
#'
#' @param object A `SingleCellExperiment` object.
#' @param genes Vector of active gene IDs to summarize.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.subset Optional vector of group labels to include.
#' @param entropy.use Diversity index: `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `"InverseSimpson"`.
#' @param assay.use Assay name to use.
#' @param entropy.thresh Threshold used to classify genes as `"monoform"` or `"polyform"`. If `NULL`, a method-specific default is used.
#' @param min.tx.cts Minimum transcript counts required before diversity is calculated.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The gene being queried.}
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
#' @importFrom Matrix rowSums

GetDiversity <- function (
    object,
    genes,
    group.by = NULL,
    group.subset = NULL,
    entropy.use = "Tsallis",
    assay.use = "counts",
    entropy.thresh = NULL,
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
  assertChoice(entropy.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(entropy.thresh, null.ok = TRUE, finite = TRUE, lower = 0)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertTRUE(order != 1 || is.null(order))
  assertFlag(quiet)

  div.func <- .DiversityFunction(entropy.use, order)
  entropy.thresh <- .DiversityThreshold(entropy.use, entropy.thresh)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.gene.id <- active_ids$active.gene.id

  # Gene filter
  gene_filter <- .FilterGenes(object, genes, active.gene.id, quiet = quiet)
  object <- gene_filter$object
  genes <- gene_filter$genes

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  unique_groups <- unique(colData(object)$group_var)
  ## check groups
  if (!is.null(group.subset)) {
    assertSubset(group.subset, unique_groups)
    ## subset object for groups
    object <- object[, object$group_var %in% group.subset, drop = FALSE]
    unique_groups <- unique(colData(object)$group_var)
  }

  # Diversity
  ## loop through each group
  res_list <- list()
  for (group in unique_groups) {

    ## subset group
    object_grp <- object[, object$group_var == group, drop = FALSE]

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

    div_res <- agg_cts_df %>%
      group_by(gene_query) %>%
      mutate(n.transcripts = n(),
             prop = cts / sum(cts),
             diversity = div.func(x = prop),
             class = ifelse(diversity <= entropy.thresh, "monoform", "polyform")
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
