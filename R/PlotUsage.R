#' Visualize isoform usage across cell groups
#'
#' Generates plots showing the isoform usage of a gene.
#'
#' @param object A `SingleCellExperiment` object.
#' @param gene A gene ID.
#' @param group.by Name of `colData` variable to group cells. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param group.subset An optional vector specifying a subset of the elements in `group.by` to include.
#' @param group.order An optional vector specifying the order of group elements.
#' @param plot.type Which type of plot to generate. Options include `"stackedbar"`, `"bar"`, `"pie"`, and `"heatmap"`.
#' @param colors A vector of colors to use for the plot.
#' @param show.prop Logical; if `TRUE`, proportion labels will be included on the plot.
#' @param nrow Number of facet rows.
#' @param assay.use Which `assay` (counts) to use.
#' @param min.tx.cts Minimum transcript counts in at least 1 group required for a transcript to be included in proportion calculations.
#' @param min.tx.prop Minimum transcript proportion required for a transcript to be plotted. Note: This parameter is intended for use with the `"heatmap"` plot type.
#' @param text.size Text size.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A ggplot object.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @import patchwork

PlotUsage <- function(
    object,
    gene,
    group.by = NULL,
    group.subset = NULL,
    group.order = NULL,
    plot.type = "stackedbar",
    colors = NULL,
    show.prop = FALSE,
    nrow = 1,
    assay.use = "counts",
    min.tx.cts = 1,
    min.tx.prop = 0,
    text.size = 12,
    quiet = FALSE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertString(gene)
  if (is.null(group.by)) {
    group.by <- metadata(object)$active.group.id
    assertChoice(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
    assertFALSE(anyMissing(colData(object)[[group.by]]))
  } else {
    assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  }
  assertCharacter(group.subset, null.ok = TRUE)
  assertCharacter(group.order, null.ok = TRUE)
  assertChoice(plot.type, c("bar", "stackedbar", "pie", "heatmap"))
  assertCharacter(colors, null.ok = TRUE)
  assertFlag(show.prop)
  assertNumber(nrow, lower = 1, finite = TRUE)
  assertTRUE(assay.use %in% assayNames(object))
  assertNumber(min.tx.cts, lower = 0, finite = TRUE)
  assertNumber(min.tx.prop, lower = 0, upper = 1, finite = TRUE)
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertNumber(nrow, lower = 1, finite = TRUE)
  assertFlag(quiet)

  # Transcript and gene IDs
  assertString(metadata(object)$active.transcript.id)
  assertString(metadata(object)$active.gene.id)
  active.transcript.id <- metadata(object)$active.transcript.id
  active.gene.id <- metadata(object)$active.gene.id
  assertChoice(active.gene.id, colnames(rowData(object)))
  assertFALSE(anyMissing(rowData(object)[[active.gene.id]]))
  assertTRUE(gene %in% rowData(object)[[active.gene.id]])

  if (metadata(object)$active.transcript.id != "") {
    assertChoice(active.transcript.id, colnames(rowData(object)))
    assertFALSE(any(duplicated(rowData(object)[[active.transcript.id]])))
    assertFALSE(anyMissing(rowData(object)[[active.transcript.id]]))
    rownames(object) <- rowData(object)[[active.transcript.id]]
  }

  gene.id.df <- rowData(object)[active.gene.id] %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    rename("gene_query" = all_of(active.gene.id))

  # Group structure
  colData(object)$group_var <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE) %>%
    pull(group_var)
  unique_groups <- unique(colData(object)$group_var)
  ## check group subset
  assertSubset(group.subset, unique_groups, empty.ok = TRUE)

  ## subset cells
  if (!is.null(group.subset)) {
    object <- object[, colData(object)$group_var %in% group.subset]
  }

  # Expression mat
  expr_mat <- assay(object[rowData(object)[[active.gene.id]] == gene], assay.use)
  group <- colData(object)[["group_var"]]
  grp_cts <- t(rowsum(t(expr_mat), group))

  ## prepare plot data
  plotdata <- grp_cts %>%
    as.data.frame() %>%
    rownames_to_column(var = "transcripts_query") %>%
    ## add gene_ids
    left_join(., gene.id.df, by = "transcripts_query") %>%
    pivot_longer(-c("transcripts_query", "gene_query"), names_to = "group_var", values_to = "grp_cts") %>%
    ## filter isoforms by counts
    filter(any(grp_cts >= min.tx.cts)) %>%
    ## calculate isoform props
    group_by(group_var) %>%
    mutate(prop = grp_cts / sum(grp_cts)) %>%
    ungroup()

  missing_groups <- plotdata %>%
    filter(is.nan(prop)) %>%
    pull(group_var) %>%
    unique()

  if (!quiet && length(missing_groups) > 0) message("\u2139 Warning: Zero isoform counts detected for '", paste0(missing_groups, collapse = "', "), "'.")

  ## order groups
  if (!is.null(group.order)) {
    avail_groups <- unique(colData(object)$group_var)
    if (length(missing_groups) > 0) {
      avail_groups <- avail_groups[-which(avail_groups %in% missing_groups)]
    }
    assertSetEqual(group.order, avail_groups)
    group_var_order <- group.order
  } else {
    if (!is.null(group.subset)) {
      group_var_order <- group.subset
    } else {
      group_var_order <- unique(colData(object)$group_var)
    }
  }
  if (plot.type != "heatmap") {
      group_var_order <- rev(group_var_order)
  }

  plotdata <- plotdata %>%
    filter(!is.nan(prop)) %>%
    ## assign 'Other'
    group_by(transcripts_query) %>%
    mutate(transcripts_query = ifelse(any(prop >= min.tx.prop), transcripts_query, "Other")) %>%
    ungroup() %>%
    group_by(transcripts_query, group_var) %>%
    mutate(prop = ifelse(transcripts_query == "Other", sum(prop), prop)) %>%
    ungroup() %>%
    distinct(transcripts_query, gene_query, group_var, prop) %>%
    ## group order
    mutate(group_var = factor(group_var, levels = rev(group_var_order)),
           group_var = droplevels(group_var))
  ## transcript order
  tx_order <- plotdata %>%
    filter(transcripts_query != "Other") %>%
    arrange(desc(prop)) %>%
    pull(transcripts_query) %>%
    unique()

  plotdata$transcripts_query <- factor(plotdata$transcripts_query, levels = c(tx_order, "Other"))

  # bar plot
  if (plot.type == "bar") {
    p1 <- ggplot(plotdata) +
      geom_col(aes(y = prop, x = transcripts_query, fill = transcripts_query)) +
      facet_wrap(~ group_var, scales = "fixed", nrow = nrow) +
      ylim(c(0, 1)) +
      xlab(active.transcript.id) +
      ylab("Proportion") +
      labs(title = gene, fill = active.transcript.id) +
      theme_linedraw(base_size = text.size) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black"),
            strip.background = element_blank())

    if (show.prop) {
      p1 <- p1 +
        geom_text(aes(y = prop, x = transcripts_query, label = ifelse(prop == 0, "", round(prop, 3))),
                  vjust = -0.1, size = text.size / .pt * 0.8)
    }
  }

  # stacked bar plot
  if (plot.type == "stackedbar") {
    # y position
    plotdata <- plotdata %>%
      group_by(group_var) %>%
      arrange(transcripts_query, .by_group = TRUE) %>%
      mutate(ypos = cumsum(prop) - prop / 2) %>%
      ungroup()

    p1 <- ggplot(plotdata) +
      geom_col(aes(y = prop, x = group_var, fill = transcripts_query), position = position_stack(reverse = TRUE)) +
      ylim(c(0, 1.01)) +
      xlab(group.by) +
      ylab("Proportion") +
      theme_linedraw(base_size = text.size) +
      labs(title = gene, fill = active.transcript.id) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))

    if (show.prop) {

      p1 <- p1 +
        geom_text(aes(y = ypos, x = group_var, label = ifelse(prop == 0, "", round(prop, 3))),
                  size = text.size / .pt * 0.8)
    }
  }

  # pie
  if (plot.type == "pie") {
    # y position
    plotdata <- plotdata %>%
      group_by(group_var) %>%
      arrange(transcripts_query, .by_group = TRUE) %>%
      mutate(ypos = cumsum(prop) - 0.5 * prop) %>%
      ungroup()

    p1 <- ggplot(plotdata) +
      geom_col(aes(y = prop, x = "", fill = transcripts_query), position = position_stack(reverse = TRUE), width = 1) +
      coord_polar("y", start = 0) +
      facet_wrap(~ group_var, scales = "fixed", nrow = nrow) +
      xlab("") +
      ylab("") +
      labs(title = gene, fill = active.transcript.id) +
      theme_linedraw(base_size = text.size) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_blank(),
            strip.text = element_text(color = "black"),
            strip.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank())

    if (show.prop) {
      p1 <- p1 +
        geom_text(aes(y = ypos, x = "", label = ifelse(prop == 0 | prop == 1, "", round(prop, 3))),
                  size = text.size / .pt * 0.8)
    }
  }

  # heatmap
  if (plot.type == "heatmap") {
    p1 <- ggplot(plotdata) +
      geom_tile(aes(y = group_var, x = transcripts_query, fill = prop), color = "black") +
      # coord_fixed() +
      xlab(active.transcript.id) +
      ylab(group.by) +
      labs(title = gene, fill = NULL) +
      theme_linedraw(base_size = text.size) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank())

    if (show.prop) {
      p1 <- p1 +
        geom_text(aes(y = group_var, x = transcripts_query, label = round(prop, 3)),
                  size = text.size / .pt)
    }
  }

  # colors
  if (plot.type != "heatmap") {
    if (!is.null(colors)) {
      p1 <- p1 +
        scale_fill_manual(values = colors, name = active.transcript.id)
    } else {
      n_colors <- length(unique(plotdata$transcripts_query))
      colors <- c("#FBB463", "#80B1D3", "#F47F72", "#BDBAD8", "#FBF8B4", "#8DD1C6")
      if (n_colors > length(colors)) {
        cont_palette <- colorRampPalette(colors)
        colors <- cont_palette(n_colors)
      }
      if ("Other" %in% plotdata$transcripts_query) {
        colors <- colors[1:n_colors]
        colors[length(colors)] <- "#9E9E9E"
      }
      p1 <- p1 +
        scale_fill_manual(values = colors, name = active.transcript.id)
    }
  }
  if (plot.type == "heatmap") {
    if (!is.null(colors)) {
      suppressMessages(
        p1 <- p1 +
          scale_fill_gradientn(colours = colors, breaks = seq(0, 1, 0.25), limits = c(0, 1), name = "Proportion",
                               guide = guide_colorbar(direction = "vertical", frame.colour = "black", ticks.colour = "black"))
      )
    } else {
      colors <- c("white", "#9C9BE9", "#FFEB3B")
      p1 <- p1 +
        scale_fill_gradientn(colours = colors, breaks = seq(0, 1, 0.25), limits = c(0, 1), name = "Proportion",
                             guide = guide_colorbar(direction = "vertical", frame.colour = "black", ticks.colour = "black"))
    }
  }

  p1

}
