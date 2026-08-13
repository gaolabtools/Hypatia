#' Visualize isoform expression
#'
#' Plots transcript expression across cell groups or embeddings.
#'
#' @param object A `SingleCellExperiment` object.
#' @param transcripts Vector of active transcript IDs to plot.
#' @param group.by One or more `colData` column names used to define cell groups. If `NULL`, `metadata(object)$active.group.id` is used.
#' @param group.subset Optional vector of group labels to include.
#' @param group.order Optional vector of group labels specifying plotting order.
#' @param plot.type Plot type: `"violin"`, `"reducedDim"`, or `"heatmap"`.
#' @param colors A vector of colors to use for the plot.
#' @param colors.heatmap A vector of colors for the heatmap color bar.
#' @param dim.use Name of the `reducedDim` slot to use when `plot.type = "reducedDim"`.
#' @param assay.use Assay name to use.
#' @param scale.heatmap Logical; if `TRUE`, standardizes and clips expression values for the heatmap.
#' @param pt.size Point size for the `violin` and `reducedDim` plots.
#' @param pt.alpha Point alpha for the `violin` and `reducedDim` plots.
#' @param text.size Text size.
#' @param label Logical; if `TRUE`, adds cell group labels to the `reducedDim` and `heatmap` plots.
#' @param show.gene Logical; if `TRUE`, shows the active gene ID alongside transcript names.
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
#' @importFrom tidyr pivot_longer
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stats median

PlotExpression <- function(
    object,
    transcripts,
    group.by = NULL,
    group.subset = NULL,
    group.order = NULL,
    plot.type = "violin",
    colors = NULL,
    colors.heatmap = NULL,
    dim.use = NULL,
    assay.use = "logcounts",
    scale.heatmap = TRUE,
    pt.size = 0.1,
    pt.alpha = 1,
    text.size = 12,
    label = TRUE,
    show.gene = TRUE,
    nrow = 1,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertCharacter(transcripts, any.missing = FALSE, unique = TRUE)
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.subset, null.ok = TRUE)
  assertCharacter(group.order, null.ok = TRUE)
  assertChoice(plot.type, c("violin", "reducedDim", "heatmap"))
  assertCharacter(colors, unique = TRUE, null.ok = TRUE)
  assertCharacter(colors.heatmap, unique = TRUE, null.ok = TRUE)
  assertChoice(dim.use, reducedDimNames(object), null.ok = TRUE)
  assertTRUE(assay.use %in% assayNames(object))
  assertFlag(scale.heatmap)
  assertNumber(pt.size, lower = 0, finite = TRUE)
  assertNumber(pt.alpha, lower = 0, upper = 1, finite = TRUE)
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertFlag(label)
  assertFlag(show.gene)
  assertNumber(nrow, lower = 1, finite = TRUE)
  assertFlag(quiet)

  # Transcript and gene IDs
  active_ids <- .ActiveIds(object)
  object <- active_ids$object
  active.transcript.id <- active_ids$active.transcript.id
  active.gene.id <- active_ids$active.gene.id

  # Transcript selection
  transcript_filter <- .FilterTranscripts(object, transcripts, quiet = quiet)
  transcripts <- transcript_filter$transcripts

  # Group structure
  colData(object)$group_var <- .GroupVar(object, group.by)
  unique_groups <- unique(colData(object)$group_var)
  group_label <- paste0(group.by, collapse = "_")
  ## check group subset
  assertSubset(group.subset, unique_groups, empty.ok = TRUE)

  ## subset cells
  if (!is.null(group.subset)) {
    object <- object[, colData(object)$group_var %in% group.subset, drop = FALSE]
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

  # Expression mat
  ## gene name option
  if (show.gene) {
    filt_object <- object[transcripts, , drop = FALSE]
    expr_mat <- assay(filt_object, assay.use)
    if (plot.type == "heatmap") {
      rownames(expr_mat) <- paste0(rowData(filt_object)[[active.gene.id]], ":", rownames(expr_mat))
    } else {
      rownames(expr_mat) <- paste0(rowData(filt_object)[[active.gene.id]], ":\n", rownames(expr_mat))
    }
  } else {
    expr_mat <- assay(object, assay.use)[transcripts, , drop = FALSE]
  }
  ## preserve transcript order
  transcripts <- rownames(expr_mat)

  # Violin plot
  if (plot.type == "violin") {
    if (!quiet) message("Generating violin plot...")

    plotdata <- expr_mat %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      mutate(grouping_var = colData(object)$group_var) %>%
      mutate(grouping_var = factor(grouping_var, levels = group_var_order)) %>%
      pivot_longer(cols = all_of(transcripts), names_to = "transcript", values_to = "expression") %>%
      mutate(transcript = factor(transcript, levels = transcripts))

    p1 <- ggplot(plotdata, aes(x = grouping_var, y = expression)) +
      geom_violin(aes(fill = grouping_var)) +
      geom_point(size = pt.size, alpha = pt.alpha, position = position_jitter(width = 0.2), show.legend = FALSE) +
      facet_wrap(~ transcript, scales = "fixed", axes = "all", axis.labels = "all", nrow = nrow) +
      labs(fill = group_label) +
      ylab(assay.use) +
      xlab(group_label) +
      theme_linedraw(base_size = text.size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black"),
            strip.background = element_blank())

    if (!is.null(colors)) {
      p1 <- p1 +
        scale_fill_manual(values = colors)
    } else {
      colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", "#ADD3F4",
                  "#F7C9CF", "#FEE4E8", "#7CA3B8", "#BFB8D6", "#FCCB8E")
      colors <- .DefaultDiscreteColors(length(group_var_order), colors)
      p1 <- p1 +
        scale_fill_manual(values = colors)
    }

  }

  # Reduced dim plot
  if (plot.type == "reducedDim") {
    if (!quiet) message("Generating reducedDim plot...")

    assertChoice(dim.use, reducedDimNames(object), null.ok = FALSE)
    coords <- reducedDim(object, dim.use)
    coords <- coords[colnames(expr_mat), , drop = FALSE]
    if (ncol(coords) < 2) {
      stop("At least two dimensions are required.")
    }
    dim_names <- colnames(coords)[1:2]
    plotdata <- expr_mat %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      mutate(grouping_var = colData(object)$group_var) %>%
      cbind(coords) %>%
      pivot_longer(cols = all_of(transcripts), names_to = "transcript", values_to = "expression") %>%
      mutate(grouping_var = factor(grouping_var, levels = group_var_order)) %>%
      group_by(grouping_var) %>%
      mutate(
        label_x = mean(.data[[dim_names[1]]]),
        label_y = mean(.data[[dim_names[2]]])
      ) %>%
      ungroup() %>%
      mutate(transcript = factor(transcript, levels = transcripts))

    p1 <- plotdata %>%
      ggplot(aes(x = .data[[dim_names[1]]], y = .data[[dim_names[2]]])) +
      geom_point(aes(color = expression), size = pt.size, alpha = pt.alpha) +
      facet_wrap(~ transcript, scales = "fixed", axes = "all", axis.labels = "all", nrow = nrow) +
      coord_cartesian(clip = "off") +
      labs(color = assay.use) +
      theme_linedraw(base_size = text.size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black"),
            strip.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

    if (is.null(colors)) {
      p1 <- p1 +
        scale_color_gradient(low = "#E0E0E0", high = "#D32F2F")
    } else {
      p1 <- p1 +
        scale_color_gradientn(colours = colors)
    }
    if (label) {
      p1 <- p1 +
        geom_text(data = distinct(plotdata, label_x, label_y, grouping_var),
                  mapping = aes(x = label_x, y = label_y, label = grouping_var), size = text.size / .pt * 0.8)
    }
  }

  # Heatmap
  if (plot.type == "heatmap") {
    if (!quiet) message("Generating heatmap...")

    plotdata <- expr_mat %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      mutate(grouping_var = colData(object)$group_var) %>%
      rownames_to_column(var = "heatmap_cellIDs") %>%
      pivot_longer(cols = all_of(transcripts), names_to = "transcript", values_to = "expression") %>%
      mutate(grouping_var = factor(grouping_var, levels = group_var_order),
             transcript = factor(transcript, levels = rev(transcripts))) %>%
      arrange(grouping_var, heatmap_cellIDs) %>%
      mutate(heatmap_cellIDs = factor(heatmap_cellIDs, levels = unique(heatmap_cellIDs)))

    ## scaling
    if (scale.heatmap) {
      plotdata <- plotdata %>%
        ungroup() %>%
        group_by(transcript) %>%
        mutate(expression = scale(expression),
               expression = pmax(pmin(expression, 2), -2)) %>%
        ungroup()
      gradient_label <- paste0(assay.use, "\nscaled")
    } else {
      gradient_label <- assay.use
    }
    ## group label positions
    group_label_pos <- plotdata %>%
      group_by(grouping_var) %>%
      summarize(x.pos = heatmap_cellIDs[floor(median(seq_along(heatmap_cellIDs)))])
    ## group annotations
    anno_bar_height <- 0.025 * length(transcripts) # 2.5% of the heatmap height
    group_anno <- plotdata %>%
      mutate(x_num = as.numeric(heatmap_cellIDs)) %>%
      distinct(grouping_var, x_num, heatmap_cellIDs) %>%
      mutate(y_pos = max(as.numeric(plotdata$transcript)) + 0.5 + anno_bar_height)
    ## colors
    if (is.null(colors)) {
      colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", "#ADD3F4",
                  "#F7C9CF", "#FEE4E8", "#7CA3B8", "#BFB8D6", "#FCCB8E")
      colors <- .DefaultDiscreteColors(length(group_var_order), colors)
    }
    if (is.null(colors.heatmap)) {
      colors.heatmap <- c("#0D47A1", "black", "#FFEB3B")
      if (scale.heatmap == FALSE) {
        colors.heatmap <- c("#F5F5F5", "#0D47A1")
      }
    }

    p1 <- ggplot(plotdata) +
      geom_tile(aes(x = heatmap_cellIDs, y = transcript, fill = expression)) +
      scale_fill_gradientn(colours = colors.heatmap, name = gradient_label) +
      new_scale_fill() +
      geom_tile(data = group_anno,
                mapping = aes(x = heatmap_cellIDs, y = y_pos, fill = grouping_var),
                height = anno_bar_height) +
      scale_fill_manual(values = colors, name = group_label) +

      ylab(active.transcript.id) +
      xlab(group_label) +
      theme_linedraw(base_size = text.size) +
      theme(axis.text.x.bottom = element_blank(),
            axis.text.x.top = element_text(angle = 45, hjust = 0),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank())

    if (label) {
      p1 <- p1 +
        scale_x_discrete(breaks = group_label_pos$x.pos,
                         labels = group_label_pos$grouping_var, position = "top")
    }
  }

  p1

}
