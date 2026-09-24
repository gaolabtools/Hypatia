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
#' @param entropy.thresh Diversity index threshold used to classify genes as monoform or polyform. If `NULL`, a default is chosen from the entropy index. Default thresholds for Tsallis and Renyi are defined only at orders 3 and 2, respectively; other orders return `NA` classifications unless a threshold is supplied.
#' @param prop.thresh Minimum within-gene transcript proportion used to define an effective isoform. Transcripts with proportions greater than or equal to this value are effective.
#' @param min.tx.cts Minimum transcript counts required before diversity is calculated.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi. At order 1, Tsallis and Renyi use their Shannon entropy limit, and NormalizedRenyi uses normalized Shannon entropy.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param cell.dispersion Logical; if `TRUE`, summarize cell-to-cell diversity among cells with positive retained-gene counts.
#' @param top.n Optional number of the most abundant isoforms to include in diversity calculations. If `NULL`, all isoforms are included. Must be at least 2 when supplied.
#' @param renormalize Logical; if `TRUE`, rescale the selected isoform proportions to sum to one before calculating diversity.
#'
#' @returns A data frame with the following columns:
#' \describe{
#'   \item{`group`}{The cell group being queried.}
#'   \item{`gene`}{The gene being queried.}
#'   \item{`gene.pct`}{Fraction of cells in `group` with expression of the gene, on the `[0, 1]` scale.}
#'   \item{`n.transcripts`}{Number of transcripts retained for the gene after count filtering in this group, before optional `top.n` selection.}
#'   \item{`n.effective`}{Number of transcripts with within-gene proportion greater than or equal to `prop.thresh`.}
#'   \item{`transcript`}{The associated transcript.}
#'   \item{`cts`}{Total counts of the transcript in `group`.}
#'   \item{`prop`}{The transcript proportion in `group`.}
#'   \item{`div`}{Isoform diversity of the gene in `group`.}
#'   \item{`div.class`}{`"monoform"` when `div` is at or below `entropy.thresh` and `"polyform"` otherwise.}
#'   \item{`cell.n`}{Number of gene-positive cells used for cell-level summaries.}
#'   \item{`cell.div.mean`, `cell.div.median`, `cell.div.sd`, `cell.div.iqr`}{Mean, median, sample standard deviation, and interquartile range of cell-level diversity.}
#' }
#' Dispersion columns are present as typed `NA` values when `cell.dispersion = FALSE`.
#' @details Counts are pooled within each cell group. Transcripts below
#' `min.tx.cts` in that group are removed before calculating within-gene
#' proportions and diversity. `gene.pct` is measured before this filtering.
#' Multiple `group.by` columns are joined with `_`; `group.subset` selects groups
#' without pooling them. Genes with no detection in a group are omitted.
#'
#' The diversity indices use transcript proportions `p` and natural logarithms:
#' - `"Shannon"`: `-sum(p * log(p))` over positive proportions.
#' - `"NormalizedShannon"`: Shannon entropy divided by `log(k)`, where `k` is
#'   the number of positive proportions.
#' - `"Tsallis"`: `(1 - sum(p^q)) / (q - 1)`, with default order `q = 3`.
#' - `"Renyi"`: `log(sum(p^alpha)) / (1 - alpha)`, with default order `alpha = 2`.
#' - `"NormalizedRenyi"`: Renyi entropy divided by `log(k)`.
#' - `"GiniSimpson"`: `1 - sum(p^2)`.
#' - `"InverseSimpson"`: `1 / sum(p^2)`.
#'
#' Power sums use positive proportions. At order 1, Tsallis and Renyi use
#' Shannon entropy, and NormalizedRenyi uses normalized Shannon entropy.
#' Normalized indices return `NA` when at most one proportion is positive.
#' If `top.n` is supplied, diversity uses the largest proportions; their sum
#' is rescaled to one only when `renormalize = TRUE`. Neither option changes
#' the reported transcript rows, `n.transcripts`, or `n.effective`.
#'
#' Monoform/polyform classes use `entropy.thresh`; effective isoform counts
#' instead count proportions at or above `prop.thresh`. Default classification
#' thresholds are 0.500 for Shannon, 0.243 for Tsallis at order 3, 0.435 for
#' Renyi at order 2, 0.348 for GiniSimpson, 1.533 for InverseSimpson, and zero
#' for normalized indices. Other Tsallis and unnormalized Renyi orders require
#' an explicit threshold to obtain classifications.
#'
#' With `cell.dispersion = TRUE`, diversity is additionally calculated in each
#' cell with positive retained-gene counts. For these cell-level summaries,
#' undefined normalized diversity from a single detected isoform is set to zero.
#' Pooled summaries are calculated directly, without bootstrap sampling.
#' @seealso [RunDIV()], [PlotDiversity()]
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
    prop.thresh = 0.2,
    min.tx.cts = 1,
    order = NULL,
    quiet = FALSE,
    cell.dispersion = FALSE,
    top.n = NULL,
    renormalize = FALSE
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
  assertNumber(entropy.thresh, lower = 0, finite = TRUE, null.ok = TRUE)
  assertNumber(prop.thresh, lower = 0, upper = 1, finite = TRUE)
  if (prop.thresh == 0) {
    stop("`prop.thresh` must be greater than 0.", call. = FALSE)
  }
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertFlag(quiet)
  assertFlag(cell.dispersion)
  assertCount(top.n, positive = TRUE, null.ok = TRUE)
  if (!is.null(top.n) && top.n < 2) {
    stop("`top.n` must be at least 2.", call. = FALSE)
  }
  assertFlag(renormalize)

  div.func <- .DiversityFunction(entropy.use, order, top.n, renormalize)
  entropy.thresh <- .DiversityThreshold(entropy.use, entropy.thresh, order)
  .DiversityThresholdMessage(entropy.use, order, entropy.thresh, quiet)

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
  if (cell.dispersion && !quiet) {
    message("Calculating cell-to-cell dispersion...")
  }
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

    ## cell-to-cell diversity dispersion
    if (cell.dispersion) {
      dispersion_object <- object_grp[agg_cts_df$transcripts_query, , drop = FALSE]
      dispersion_data <- .CellDiversityDispersion(
        assay(dispersion_object, assay.use),
        rowData(dispersion_object)[[active.gene.id]],
        div.func
      ) %>%
        rename(gene_query = gene.id)
    } else {
      dispersion_data <- data.frame(
        gene_query = unique(agg_cts_df$gene_query),
        cell.n = rep(NA_integer_, length(unique(agg_cts_df$gene_query))),
        cell.div.mean = rep(NA_real_, length(unique(agg_cts_df$gene_query))),
        cell.div.median = rep(NA_real_, length(unique(agg_cts_df$gene_query))),
        cell.div.sd = rep(NA_real_, length(unique(agg_cts_df$gene_query))),
        cell.div.iqr = rep(NA_real_, length(unique(agg_cts_df$gene_query)))
      )
    }

    div_res <- agg_cts_df %>%
      group_by(gene_query) %>%
      mutate(n.transcripts = n(),
            prop = cts / sum(cts),
            diversity = as.numeric(div.func(x = prop)),
            n.effective = .EffectiveIsoformCount(prop, prop.thresh),
            class = as.character(.DiversityClass(diversity, entropy.thresh))
      ) %>%
      ungroup() %>%
      left_join(dispersion_data, by = "gene_query") %>%
      mutate(prop = replace(prop, is.nan(prop), NA_real_),
             diversity = replace(diversity, is.na(prop), NA_real_)) %>%
      select(group_var, gene_query, gene.pct, n.transcripts, n.effective, transcripts_query,
             cts, prop, diversity, class, cell.n, cell.div.mean,
             cell.div.median, cell.div.sd, cell.div.iqr) %>%
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
