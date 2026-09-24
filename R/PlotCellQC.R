#' Visualize single-cell QC metrics
#'
#' Plots `nTranscript`, `nGene`, and `nCount` QC metrics by cell group.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by One or more `colData` column names used to define cell groups.
#' @param colors A vector of cell-group colors, optionally named by group label.
#'   If `NULL`, uses the shared categorical group palette: ten fixed colors,
#'   or the qualitative HCL `"Dark 3"` palette for more than ten groups.
#' @param pt.size Point size.
#' @param pt.alpha Point alpha.
#' @param text.size Text size.
#' @param show.legend Logical; if `TRUE`, the legend will be shown.
#' @param combine Logical; if `TRUE`, combines plots using `patchwork`.
#'
#' @returns A patchwork object when `combine = TRUE`; otherwise a named list of ggplot objects.
#' @details Three violin plots show the distributions of detected transcripts
#' (`nTranscript`), detected genes (`nGene`), and total isoform counts (`nCount`)
#' with individual cells overlaid. A fourth panel plots `nTranscript` against
#' `nCount` and reports their Pearson correlation across all plotted cells.
#'
#' The plots use the QC values already stored in `colData(object)`. Grouping
#' defaults to `project`; multiple `group.by` columns are joined with `_`.
#' Use `combine = FALSE` to customize the four panels independently.
#' @seealso [CreateSCE()], [SubsetCells()], [SubsetTranscripts()]
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom stats cor

PlotCellQC <- function(
    object,
    group.by = "project",
    colors = NULL,
    pt.size = 0.2,
    pt.alpha = 1,
    text.size = 12,
    show.legend = TRUE,
    combine = TRUE
) {

  # Check inputs
  assertClass(object, "SingleCellExperiment")
  assertSubset(group.by, c(setdiff(names(colData(object)), c("nCount", "nTranscript", "nGene"))))
  assertCharacter(colors, null.ok = TRUE)
  assertNumber(pt.size, lower = 0, finite = TRUE)
  assertNumber(pt.alpha, lower = 0, upper = 1, finite = TRUE)
  assertNumber(text.size, lower = 0, finite = TRUE)
  assertFlag(show.legend)
  assertFlag(combine)

  # Group structure
  plotdata <- as.data.frame(colData(object))
  names(plotdata) <- names(colData(object))
  plotdata[["group_var"]] <- .GroupVar(object, group.by)

  # nTranscript plot
  p_ntranscript <- ggplot(plotdata) +
    geom_violin(aes(y = nTranscript, x = group_var, fill = group_var)) +
    geom_point(aes(y = nTranscript, x = group_var),
               position = "jitter", pch = 16, size = pt.size, alpha = pt.alpha) +
    labs(title = "nTranscript", fill = paste0(group.by, collapse = "_")) +
    xlab(paste0(group.by, collapse = "_")) +
    ylab("nTranscript") +
    theme_linedraw(base_size = text.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # nGene plot
  p_ngene <- ggplot(plotdata) +
    geom_violin(aes(y = nGene, x = group_var, fill = group_var)) +
    geom_point(aes(y = nGene, x = group_var),
               position = "jitter", pch = 16, size = pt.size, alpha = pt.alpha) +
    labs(title = "nGene", fill = paste0(group.by, collapse = "_")) +
    xlab(paste0(group.by, collapse = "_")) +
    ylab("nGene") +
    theme_linedraw(base_size = text.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # nCount plot
  p_ncount <- ggplot(plotdata) +
    geom_violin(aes(y = nCount, x = group_var, fill = group_var)) +
    geom_point(aes(y = nCount, x = group_var),
               position = "jitter", pch = 16, size = pt.size, alpha = pt.alpha) +
    labs(title = "nCount", fill = paste0(group.by, collapse = "_")) +
    xlab(paste0(group.by, collapse = "_")) +
    ylab("nCount") +
    theme_linedraw(base_size = text.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Scatter plot
  pcor <- round(cor(y = plotdata$nTranscript, x = plotdata$nCount, method = "pearson"), 3)

  p_scatt <- ggplot(plotdata) +
    geom_point(aes(y = nTranscript, x = nCount, color = group_var),
               pch = 16, size = pt.size) +
    labs(title = paste0("r = ", pcor), color = paste0(group.by, collapse = "_")) +
    xlab("nCount") +
    ylab("nTranscript") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_linedraw(base_size = text.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Legend
  if (!show.legend) {
    p_scatt <- p_scatt + theme(legend.position = "none")
  }
  if (!is.null(colors)) {
    p_ntranscript <- p_ntranscript + scale_fill_manual(values = colors)
    p_ngene <- p_ngene + scale_fill_manual(values = colors)
    p_ncount <- p_ncount + scale_fill_manual(values = colors)
    p_scatt <- p_scatt + scale_color_manual(values = colors)
  } else {
    n_colors <- length(unique(plotdata$group_var))
    colors <- .DefaultGroupColors(n_colors)
    p_ntranscript <- p_ntranscript + scale_fill_manual(values = colors)
    p_ngene <- p_ngene + scale_fill_manual(values = colors)
    p_ncount <- p_ncount + scale_fill_manual(values = colors)
    p_scatt <- p_scatt + scale_color_manual(values = colors)

  }

  # Combine
  if (combine) {
    wrap_plots(p_ntranscript, p_ngene, p_ncount, p_scatt, nrow = 1) + plot_layout(guides = "collect")
  } else {
    list("nTranscript" = p_ntranscript, "nGene" = p_ngene, "nCount" = p_ncount, "scatter" = p_scatt)
  }
}
