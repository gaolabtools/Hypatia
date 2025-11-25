#' Visualize isoform diversity
#'
#' Generates plots showing the isoform diversity of one or more genes.
#'
#' @param object A `SingleCellExperiment` object.
#' @param genes A vector of one or more gene IDs.
#' @param group.by Name of `colData` variable to group cells. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.subset An optional vector specifying a subset of the elements in `group.by` to include.
#' @param group.order An optional vector specifying the order of group elements.
#' @param plot.type Which type of plot to generate. Options include `"lollipop"`, `"density"`, and `"pcoord"`.
#' @param method.use The diversity index to calculate.
#' Options include `"Tsallis"`, `"Shannon"`, `"NormalizedShannon"`, `"Renyi"`, `"NormalizedRenyi"`, `"GiniSimpson"`, or `InverseSimpson`.
#' @param assay.use Which `assay` (counts) to use.
#' @param diversity.cutoff The cutoff of the diversity index used for monoform and polyform classification.
#' Default cutoffs are 0.243 for Tsallis, 0.500 for Shannon, 0 for normalized Shannon, 0.435 for Renyi, 0 for normalized Renyi, 0.348 for Gini-Simpson, and 1.533 for inverse Simpson.
#' @param min.tx.cts Minimum transcript counts required for a transcript to be included in the contingency table.
#' @param order Value specifying the order of entropy. Corresponds to `q` for Tsallis (default: 3) and `alpha` for Renyi (default: 2).
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
#' @importFrom tidyr unite

PlotDiversity <- function (
    object,
    genes,
    group.by = NULL,
    group.subset = NULL,
    group.order = NULL,
    plot.type = "lollipop",
    method.use = "Tsallis",
    assay.use = "counts",
    diversity.cutoff = NULL,
    min.tx.cts = 1,
    order = NULL,
    colors = NULL,
    text.size = 12,
    nrow = 1,
    verbose = TRUE
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
  assertChoice(method.use, c("Tsallis", "Shannon", "NormalizedShannon", "Renyi", "NormalizedRenyi", "GiniSimpson", "InverseSimpson"))
  assertCharacter(colors, null.ok = TRUE)
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(order, lower = 0, finite = TRUE, null.ok = TRUE)
  assertTRUE(order != 1 || is.null(order))
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertNumber(nrow, lower = 1, finite = TRUE)
  assertFlag(verbose)

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
      if (verbose) message("\u2139 Warning: The following genes were not found in the object: '", paste0(missing_genes, collapse = "', '"), "'.")
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
      dplyr::filter(cts >= min.tx.cts) %>%
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
      mutate(prop = cts / sum(cts),
             diversity = div.func(x = prop)) %>%
      ungroup() %>%
      mutate(prop = ifelse(is.nan(prop), NA, prop),
             diversity = ifelse(is.na(prop), NA, diversity)) %>%
      distinct(group_var, gene_query, diversity) %>%
      mutate(class = ifelse(diversity <= diversity.cutoff, "monoform", "polyform"))

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
    if (n_group_colors > length(group_colors)) {
      group_cont_palette <- colorRampPalette(group_colors)
      group_colors <- group_cont_palette(n_group_colors)
    }

    gene_colors <- c("#FBB463", "#80B1D3", "#F47F72", "#BDBAD8", "#FBF8B4", "#8DD1C6")
    n_gene_colors <- length(unique(plotdata$gene_query))
    gene_cont_palette <- colorRampPalette(gene_colors)
    gene_colors <- gene_cont_palette(n_gene_colors)
    if (n_gene_colors > length(gene_colors)) {
      group_cont_palette <- colorRampPalette(group_colors)
      group_colors <- group_cont_palette(n_group_colors)
    }
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
      labs(color = group.by, linetype = "class") +
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
      xlab(group.by) +
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
      labs(fill = group.by) +
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
