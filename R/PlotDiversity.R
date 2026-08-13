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
#' @param entropy.thresh Threshold used to classify genes as `"monoform"` or `"polyform"`. If `NULL`, a method-specific default is used.
#' @param min.tx.cts Minimum transcript counts required before diversity is calculated.
#' @param order Entropy order. Corresponds to `q` for Tsallis and `alpha` for Renyi.
#' @param colors A vector of colors to use for the plot.
#' @param text.size Text size.
#' @param nrow Number of facet rows.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A ggplot object.
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
    min.tx.cts = 1,
    order = NULL,
    colors = NULL,
    text.size = 12,
    nrow = 1,
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
  assertCharacter(group.order, null.ok = TRUE)
  assertChoice(plot.type, c("lollipop", "density", "pcoord"))
  assertNumber(nrow, lower = 1, finite = TRUE)
  assertTRUE(assay.use %in% assayNames(object))
  assertChoice(entropy.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertCharacter(colors, null.ok = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertTRUE(order != 1 || is.null(order))
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertNumber(nrow, lower = 1, finite = TRUE)
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
             diversity = div.func(x = prop)) %>%
      ungroup() %>%
      mutate(prop = ifelse(is.nan(prop), NA, prop),
             diversity = ifelse(is.na(prop), NA, diversity)) %>%
      distinct(group_var, gene_query, diversity) %>%
      mutate(class = ifelse(diversity <= entropy.thresh, "monoform", "polyform"))

    res_list[[group]] <- div_res
  }

  plotdata <- purrr::reduce(res_list, rbind)

  # Group order
  plotdata <- plotdata %>%
    mutate(group_var = factor(group_var, levels = group_var_order),
           gene_query = factor(gene_query, levels = genes))

  # Colors
  if (is.null(colors)) {
    group_colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", "#ADD3F4",
                      "#F7C9CF", "#FEE4E8", "#7CA3B8", "#BFB8D6", "#FCCB8E")
    n_group_colors <- length(unique(plotdata$group_var))
    group_colors <- .DefaultDiscreteColors(n_group_colors, group_colors)

    gene_colors <- c("#FBB463", "#80B1D3", "#F47F72", "#BDBAD8", "#FBF8B4", "#8DD1C6")
    n_gene_colors <- length(unique(plotdata$gene_query))
    gene_colors <- .DefaultDiscreteColors(n_gene_colors, gene_colors)
  }

  # Lollipop plot
  if (plot.type == "lollipop") {

    p1 <- plotdata %>%
      ggplot() +
      geom_linerange(aes(x = gene_query, ymin = 0, ymax = diversity, color = group_var, linetype = class),
                     position = position_dodge(width = 0.8)) +
      geom_point(aes(x = gene_query, y = diversity, group = group_var),
                 size = 3.5, color = "black", position = position_dodge(width = 0.8)) +
      geom_point(aes(x = gene_query, y = diversity, color = group_var),
                 size = 3, position = position_dodge(width = 0.8)) +
      scale_linetype_manual(values = c("monoform" = "solid", "polyform" = "dashed"),
                            breaks = c("monoform", "polyform"), na.translate = FALSE) +
      labs(color = group_label, linetype = "class") +
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
      labs(color = active.gene.id, linetype = "class") +
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
      geom_density(aes(x = diversity, fill = group_var)) +
      facet_wrap(~ group_var, nrow = nrow) +
      labs(fill = group_label) +
      xlab("Diversity") +
      ylab("Density") +
      theme_linedraw(base_size = text.size) +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(color = "black"))

    if (is.null(colors)) {
      p1 <- p1 +
        scale_fill_manual(values = group_colors)
    } else {
      p1 <- p1 +
        scale_fill_manual(values = colors)
    }

  }

  p1

}
