#' Visualize isoform diversity
#'
#' Plots isoform diversity for one or more genes across cell groups.
#'
#' @param object A `SingleCellExperiment` object.
#' @param genes Vector of active gene IDs to plot.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.subset Optional vector of group labels to include.
#' @param group.order Optional vector of group labels specifying plotting order.
#' @param plot.type Plot type: `"lollipop"`, `"density"`, or `"pcoord"`.
#' @param entropy.use Diversity index: `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `"InverseSimpson"`.
#' @param assay.use Assay name to use.
#' @param entropy.thresh Diversity index threshold used to classify genes as monoform or polyform. If `NULL`, a default is chosen from the entropy index. Default thresholds for Tsallis and Renyi are defined only at orders 3 and 2, respectively; other orders return `NA` classifications unless a threshold is supplied.
#' @param prop.thresh Minimum within-gene transcript proportion used to define an effective isoform. Transcripts with proportions greater than or equal to this value are effective.
#' @param min.tx.cts Minimum transcript counts required before diversity is calculated.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi. At order 1, Tsallis and Renyi use their Shannon entropy limit, and NormalizedRenyi uses normalized Shannon entropy.
#' @param colors A vector of colors, optionally named by group label (lollipop
#'   and density) or gene ID (parallel coordinates). If `NULL`, cell groups use
#'   the shared categorical group palette: ten fixed colors, or the qualitative
#'   HCL `"Dark 3"` palette for more than ten groups. Parallel-coordinate plots
#'   use a separate gene palette.
#' @param text.size Text size.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#' @param top.n Optional number of the most abundant isoforms to include in diversity calculations. If `NULL`, all isoforms are included. Must be at least 2 when supplied.
#' @param renormalize Logical; if `TRUE`, rescale the selected isoform proportions to sum to one before calculating diversity.
#'
#' @returns A ggplot object. The `class` aesthetic reflects whether `diversity` is at or below `entropy.thresh`.
#' @details Diversity is calculated from pooled transcript counts within each
#' group after applying `min.tx.cts` separately in that group. The entropy
#' indices, `top.n` selection, optional renormalization, and classification
#' thresholds follow [GetDiversity()]. These are direct pooled estimates,
#' rather than the bootstrap means returned by [RunDIV()]. Multiple `group.by`
#' columns are joined with `_`.
#'
#' - `"lollipop"`: Diversity by gene, with colored points for each group.
#' - `"pcoord"`: Diversity across groups, with a line for each gene.
#' - `"density"`: The distribution of gene-level diversity values in each group.
#'
#' Point shapes in lollipop and parallel-coordinate plots indicate the
#' entropy-based monoform/polyform class, not the number of effective isoforms.
#' Supply an appropriate `entropy.thresh` when using an entropy order without
#' a default classification threshold. Density plots summarize variation across
#' genes, not cell-to-cell variation or bootstrap uncertainty.
#' @seealso [GetDiversity()], [RunDIV()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @import ggplot2
#' @importFrom Matrix rowSums

PlotDiversity <- function (
    object,
    genes,
    group.by = NULL,
    group.subset = NULL,
    group.order = NULL,
    plot.type = "lollipop",
    entropy.use = "Tsallis",
    assay.use = "counts",
    entropy.thresh = NULL,
    prop.thresh = 0.2,
    min.tx.cts = 1,
    order = NULL,
    colors = NULL,
    text.size = 12,
    quiet = FALSE,
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
  assertCharacter(group.order, null.ok = TRUE)
  assertChoice(plot.type, c("lollipop", "density", "pcoord"))
  assertTRUE(assay.use %in% assayNames(object))
  assertChoice(entropy.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertNumber(entropy.thresh, lower = 0, finite = TRUE, null.ok = TRUE)
  assertNumber(prop.thresh, lower = 0, upper = 1, finite = TRUE)
  if (prop.thresh == 0) {
    stop("`prop.thresh` must be greater than 0.", call. = FALSE)
  }
  assertCharacter(colors, null.ok = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertFlag(quiet)
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
  group_label <- paste0(group.by, collapse = "_")
  ## check groups
  if (!is.null(group.subset)) {
    assertSubset(group.subset, unique_groups)
    ## subset object for groups
    object <- object[, object$group_var %in% group.subset, drop = FALSE]
    unique_groups <- unique(colData(object)$group_var)
  }
  ## order groups
  if (!is.null(group.order)) {
    assertSetEqual(group.order, unique(colData(object)$group_var))
    group_var_order <- group.order
  } else {
    if (!is.null(group.subset)) {
      group_var_order <- group.subset
    } else {
      group_var_order <- unique(colData(object)$group_var)
    }
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
      dplyr::filter(cts >= min.tx.cts) %>%
      mutate("group_var" = group)

     div_res <- agg_cts_df %>%
      group_by(gene_query) %>%
      mutate(prop = cts / sum(cts),
             diversity = div.func(x = prop),
             n.effective = .EffectiveIsoformCount(prop, prop.thresh),
             class = .DiversityClass(diversity, entropy.thresh)) %>%
      ungroup() %>%
      mutate(prop = ifelse(is.nan(prop), NA, prop),
             diversity = ifelse(is.na(prop), NA, diversity)) %>%
      distinct(group_var, gene_query, diversity, n.effective, class)

    res_list[[group]] <- div_res
  }

  plotdata <- purrr::reduce(res_list, rbind)

  # Group order
  plotdata <- plotdata %>%
    mutate(group_var = factor(group_var, levels = group_var_order),
           gene_query = factor(gene_query, levels = genes))

  # Colors
  if (is.null(colors)) {
    n_group_colors <- length(unique(plotdata$group_var))
    group_colors <- .DefaultGroupColors(n_group_colors)

    gene_colors <- c("#FBB463", "#80B1D3", "#F47F72", "#BDBAD8", "#FBF8B4", "#8DD1C6")
    n_gene_colors <- length(unique(plotdata$gene_query))
    gene_colors <- .DefaultDiscreteColors(n_gene_colors, gene_colors)
  }

  # Lollipop plot
  if (plot.type == "lollipop") {

    p1 <- plotdata %>%
      ggplot() +
      geom_linerange(aes(x = gene_query, ymin = 0, ymax = diversity, color = group_var),
                     position = position_dodge(width = 0.8)) +
      geom_point(aes(x = gene_query, y = diversity, color = group_var, shape = class),
                 size = 3.5, fill = "white", position = position_dodge(width = 0.8)) +
      scale_shape_manual(values = c("monoform" = 19, "polyform" = 21),
                         breaks = c("monoform", "polyform"), na.translate = FALSE) +
      labs(color = group_label, shape = "class") +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2)) +
      xlab(active.gene.id) +
      ylab("Diversity") +
      theme_linedraw(base_size = text.size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank())

    if (is.null(colors)) {
      p1 <- p1 +
        scale_color_manual(values = group_colors)
    } else {
      p1 <- p1 +
        scale_color_manual(values = colors)
    }

  }

  # Parallel coord plot
  if (plot.type == "pcoord") {
    p1 <- plotdata %>%
      ggplot() +
      geom_line(aes(x = group_var, y = diversity, color = gene_query, group = gene_query)) +
      geom_point(aes(x = group_var, y = diversity, color = gene_query, shape = class),
                 size = 3.5, fill = "white") +
      scale_shape_manual(values = c("monoform" = 19, "polyform" = 21),
                         breaks = c("monoform", "polyform"), na.translate = FALSE) +
      labs(color = active.gene.id, shape = "class") +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2)) +
      xlab(group_label) +
      ylab("Diversity") +
      theme_linedraw(base_size = text.size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank())

    if (is.null(colors)) {
      p1 <- p1 +
        scale_color_manual(values = gene_colors)
    } else {
      p1 <- p1 +
        scale_color_manual(values = colors)
    }
  }

  # Density plot
  if (plot.type == "density") {
    p1 <- plotdata %>%
      ggplot() +
      geom_density(aes(x = diversity, color = group_var), linewidth = 1) +
      labs(color = group_label) +
      xlab("Diversity") +
      ylab("Density") +
      theme_linedraw(base_size = text.size) +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(color = "black"))

    if (is.null(colors)) {
      p1 <- p1 +
        scale_color_manual(values = group_colors)
    } else {
      p1 <- p1 +
        scale_color_manual(values = colors)
    }

  }

  p1

}
