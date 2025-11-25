#' Visualize single-cell QC metrics
#'
#' Generates plots of single-cell QC metrics `nTranscript`, `nGene`, and `nCount` for cell filtering. `nCount` corresponds to transcript counts per cell.
#'
#' @param object A `SingleCellExperiment` object.
#' @param group.by Name of `colData` variable to group cells. If `NULL`, `metadata(object)$active.group.id` will be used.
#' @param colors A vector of colors to use for the plot.
#' @param pt.size Point size.
#' @param pt.alpha Point alpha.
#' @param text.size Text size.
#' @param show.legend Logical; if `TRUE`, the legend will be shown.
#' @param combine Logical; if `TRUE`, combines plots using `patchwork`.
#'
#' @returns A patchwork/ggplot object.
#' @export
#' @import checkmate
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom tidyr unite
#' @import ggplot2
#' @import patchwork

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
  plotdata <- colData(object) %>%
    as.data.frame() %>%
    unite("group_var", all_of(group.by), sep = "_", remove = FALSE)

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
    xlab("nTranscript") +
    ylab("nCount") +
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
    colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", "#ADD3F4",
                "#F7C9CF", "#FEE4E8", "#7CA3B8", "#BFB8D6", "#FCCB8E")
    n_colors <- length(unique(plotdata$group_var))
    if (n_colors > length(colors)) {
      cont_palette <- colorRampPalette(colors)
      colors <- cont_palette(n_colors)
    }
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
