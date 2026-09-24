.ValidateTranscriptExons <- function(exons, transcript.ids) {
  assertClass(exons, "GRangesList")
  assertTRUE(identical(names(exons), transcript.ids))
  if (any(lengths(exons) == 0L)) {
    stop("Missing exon ranges for transcripts: ",
         paste(transcript.ids[lengths(exons) == 0L], collapse = ", "), call. = FALSE)
  }
  flat_exons <- unlist(exons, use.names = FALSE)
  if (any(GenomicRanges::width(flat_exons) <= 0L)) {
    stop("Exon ranges must have positive widths.", call. = FALSE)
  }

  transcript_index <- rep(seq_along(exons), lengths(exons))
  chromosomes <- split(as.character(GenomicRanges::seqnames(flat_exons)), transcript_index)
  strands <- split(as.character(GenomicRanges::strand(flat_exons)), transcript_index)
  consistent <- lengths(lapply(chromosomes, unique)) == 1L &
    lengths(lapply(strands, unique)) == 1L
  if (any(!consistent)) {
    stop("Exons of each transcript must share a chromosome and strand.", call. = FALSE)
  }

  # Store the union of exon bases within each transcript in genomic order.
  GenomicRanges::reduce(exons)
}

.TranscriptExons <- function(gtf, transcript.ids, gtf.transcript.id) {
  assertClass(gtf, "GRanges")
  assertString(gtf.transcript.id)
  assertChoice(gtf.transcript.id, names(mcols(gtf)))
  assertCharacter(transcript.ids, any.missing = FALSE, unique = TRUE)

  if ("type" %in% names(mcols(gtf))) {
    feature_type <- as.character(mcols(gtf)$type)
    exons <- gtf[!is.na(feature_type) & feature_type == "exon"]
  } else {
    # Untyped annotations are interpreted as exon ranges supplied by the caller.
    exons <- gtf
  }
  exon_ids <- as.character(mcols(exons)[[gtf.transcript.id]])
  if (anyNA(exon_ids) || any(!nzchar(exon_ids))) {
    stop("Exon transcript IDs must be non-missing and nonempty.", call. = FALSE)
  }
  missing_ids <- setdiff(transcript.ids, exon_ids)
  if (length(missing_ids) > 0L) {
    stop("Missing exon annotations for transcripts: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  keep <- exon_ids %in% transcript.ids
  exons <- exons[keep]
  names(exons) <- NULL
  mcols(exons) <- NULL
  exons <- GenomicRanges::split(exons, factor(exon_ids[keep], levels = transcript.ids))
  .ValidateTranscriptExons(exons, transcript.ids)
}

.StoreTranscriptGTF <- function(object, gtf, gtf.transcript.id) {
  exons <- .TranscriptExons(gtf, rownames(object), gtf.transcript.id)
  annotations <- rowData(object)
  rowRanges(object) <- exons
  rowData(object) <- annotations
  metadata(object)$GTF <- gtf
  metadata(object)$gtf.transcript.id <- gtf.transcript.id
  object
}

.GroupVar <- function(object, group.by, sep = "_") {
  group_data <- as.data.frame(colData(object))
  names(group_data) <- names(colData(object))
  group_data <- group_data[, group.by, drop = FALSE]
  do.call(paste, c(group_data, sep = sep))
}

.ActiveIds <- function(object, use.transcript.id = TRUE) {
  assertString(metadata(object)$active.transcript.id)
  assertString(metadata(object)$active.gene.id)

  active.transcript.id <- metadata(object)$active.transcript.id
  active.gene.id <- metadata(object)$active.gene.id

  assertChoice(active.gene.id, colnames(rowData(object)))
  assertFALSE(anyMissing(rowData(object)[[active.gene.id]]))

  if (use.transcript.id && active.transcript.id != "") {
    assertChoice(active.transcript.id, colnames(rowData(object)))
    assertFALSE(any(duplicated(rowData(object)[[active.transcript.id]])))
    assertFALSE(anyMissing(rowData(object)[[active.transcript.id]]))
    rownames(object) <- rowData(object)[[active.transcript.id]]
  }

  list(
    object = object,
    active.transcript.id = active.transcript.id,
    active.gene.id = active.gene.id
  )
}

.FilterGenes <- function(object, genes, active.gene.id, quiet = FALSE) {
  available_genes <- unique(rowData(object)[[active.gene.id]])
  if (!any(genes %in% available_genes)) {
    stop("None of the genes were found in the object. (Check active.gene.id?)")
  }

  missing_genes <- setdiff(genes, available_genes)
  if (length(missing_genes) > 0) {
    if (!quiet) {
      message("\u2139 Warning: The following genes were not found in the object: '", paste0(missing_genes, collapse = "', '"), "'.")
    }
    genes <- genes[genes %in% available_genes]
  }

  list(
    object = object[rowData(object)[[active.gene.id]] %in% genes, , drop = FALSE],
    genes = genes
  )
}

.FilterTranscripts <- function(
    object,
    transcripts,
    quiet = FALSE,
    min.valid = 1,
    none.message = "None of the transcripts provided were found in the object. (Check active.transcript.id?)"
) {
  if (!any(transcripts %in% rownames(object))) {
    stop(none.message)
  }

  missing_transcripts <- setdiff(transcripts, rownames(object))
  if (length(missing_transcripts) > 0) {
    if (!quiet) {
      message("\u2139 Warning: The following transcripts were not found in the object: '", paste0(missing_transcripts, collapse = "', '"), "'.")
    }
    transcripts <- transcripts[transcripts %in% rownames(object)]
  }

  if (length(transcripts) < min.valid) {
    stop("Please provide at least ", min.valid, " valid transcripts to test.")
  }

  list(
    object = object[transcripts, , drop = FALSE],
    transcripts = transcripts
  )
}

.BuildGroupComparisons <- function(object, group.1 = NULL, group.2 = NULL, unique_groups = NULL) {
  if (is.null(unique_groups)) {
    unique_groups <- unique(colData(object)$group_var)
  }

  object_grp_list <- list()
  mode <- NULL

  if (is.null(group.1) && is.null(group.2)) {
    mode <- "all"
    for (grp in unique_groups) {
      group.1.current <- grp
      group.2.current <- setdiff(unique_groups, group.1.current)

      object_grp_list[[grp]] <- list(
        "grp1.object" = object[, object$group_var == grp, drop = FALSE],
        "grp2.object" = object[, object$group_var != grp, drop = FALSE],
        "grp1.names" = paste0(group.1.current, collapse = ","),
        "grp2.names" = paste0(group.2.current, collapse = ",")
      )

      if (length(unique_groups) == 2) {
        break
      }
    }
  } else if (!is.null(group.1) && is.null(group.2)) {
    mode <- "one_vs_all"
    group.2.current <- setdiff(unique_groups, group.1)

    object_grp_list[["single_test"]] <- list(
      "grp1.object" = object[, object$group_var %in% group.1, drop = FALSE],
      "grp2.object" = object[, object$group_var %in% group.2.current, drop = FALSE],
      "grp1.names" = paste0(group.1, collapse = ","),
      "grp2.names" = paste0(group.2.current, collapse = ",")
    )
  } else if (!is.null(group.1) && !is.null(group.2)) {
    mode <- "pair"
    object_grp_list[["single_test"]] <- list(
      "grp1.object" = object[, object$group_var %in% group.1, drop = FALSE],
      "grp2.object" = object[, object$group_var %in% group.2, drop = FALSE],
      "grp1.names" = paste0(group.1, collapse = ","),
      "grp2.names" = paste0(group.2, collapse = ",")
    )
  } else {
    stop("`group.1` must be specified prior to `group.2`")
  }

  valid_comparisons <- vapply(object_grp_list, function(comparison) {
    ncol(comparison$grp1.object) > 0 && ncol(comparison$grp2.object) > 0
  },
  logical(1)
  )
  if (any(!valid_comparisons)) {
    stop("Each comparison must contain at least one cell in both groups.", call. = FALSE)
  }

  attr(object_grp_list, "mode") <- mode
  object_grp_list
}

.DiversityFunction <- function(entropy.use, order = NULL, top.n = NULL,
                               renormalize = FALSE) {
  force(entropy.use)
  force(order)
  force(top.n)
  force(renormalize)

  function(x) {
    if (anyNA(x)) return(NA_real_)
    if (!is.null(top.n)) {
      x <- head(sort(x, decreasing = TRUE), top.n)
    }
    if (renormalize) {
      selected_mass <- sum(x)
      if (!is.finite(selected_mass) || selected_mass <= 0) return(NA_real_)
      x <- x / selected_mass
    }

    if (entropy.use == "Shannon") {
      -sum(x[x > 0] * log(x[x > 0]))
    } else if (entropy.use == "NormalizedShannon") {
      n_x <- sum(x > 0)
      if (n_x <= 1) return(NA_real_)
      (-sum(x[x > 0] * log(x[x > 0]))) / log(n_x)
    } else if (entropy.use == "Renyi") {
      order.use <- order
      if (is.null(order.use)) order.use <- 2
      if (order.use == 1) {
        -sum(x[x > 0] * log(x[x > 0]))
      } else {
        (1 / (1 - order.use)) * log(sum((x[x > 0])^order.use))
      }
    } else if (entropy.use == "NormalizedRenyi") {
      order.use <- order
      if (is.null(order.use)) order.use <- 2
      n_x <- sum(x > 0)
      if (n_x <= 1) return(NA_real_)
      if (order.use == 1) {
        (-sum(x[x > 0] * log(x[x > 0]))) / log(n_x)
      } else {
        (1 / (1 - order.use)) * log(sum((x[x > 0])^order.use)) / log(n_x)
      }
    } else if (entropy.use == "GiniSimpson") {
      1 - sum((x[x > 0])^2)
    } else if (entropy.use == "Tsallis") {
      order.use <- order
      if (is.null(order.use)) order.use <- 3
      if (order.use == 1) {
        -sum(x[x > 0] * log(x[x > 0]))
      } else {
        (1 - sum(x[x > 0]^order.use)) / (order.use - 1)
      }
    } else if (entropy.use == "InverseSimpson") {
      1 / sum((x[x > 0])^2)
    }
  }
}

.DiversityThreshold <- function(entropy.use, entropy.thresh = NULL,
                                order = NULL) {
  if (!is.null(entropy.thresh)) return(entropy.thresh)
  if (entropy.use == "Tsallis" && !is.null(order) && order != 3) return(NA_real_)
  if (entropy.use == "Renyi" && !is.null(order) && order != 2) return(NA_real_)

  switch(
    entropy.use,
    Shannon = 0.500,
    NormalizedShannon = 0,
    Renyi = 0.435,
    NormalizedRenyi = 0,
    GiniSimpson = 0.348,
    Tsallis = 0.243,
    InverseSimpson = 1.533,
    0
  )
}

.DiversityThresholdMessage <- function(entropy.use, order, entropy.thresh,
                                       quiet) {
  if (quiet) return(invisible(NULL))

  if (is.na(entropy.thresh)) {
    message(
      "No default entropy.thresh is defined for ", entropy.use,
      " at order = ", order,
      "; monoform/polyform classifications will be NA unless ",
      "entropy.thresh is supplied."
    )
  } else {
    message(
      "Using entropy.thresh = ", entropy.thresh,
      " for monoform/polyform classification."
    )
  }
  invisible(NULL)
}

.EffectiveIsoformCount <- function(prop, prop.thresh) {
  if (any(!is.finite(prop))) {
    return(NA_integer_)
  }
  as.integer(sum(prop >= prop.thresh))
}

.EffectiveIsoformClass <- function(n.effective) {
  ifelse(
    is.na(n.effective),
    NA_character_,
    ifelse(n.effective == 1L, "monoform", "polyform")
  )
}

.DiversityClass <- function(div, entropy.thresh) {
  ifelse(
    is.na(div),
    NA_character_,
    ifelse(div <= entropy.thresh, "monoform", "polyform")
  )
}

#' @importFrom grDevices colorRampPalette
.DefaultDiscreteColors <- function(n, colors) {
  if (n > length(colors)) {
    colorRampPalette(colors)(n)
  } else {
    colors[seq_len(n)]
  }
}

.DefaultGroupColors <- function(n) {
  colors <- c("#4E79A7", "#E15759", "#59A14F", "#B07AA1", "#F28E2B",
              "#76B7B2", "#B6992D", "#9C755F", "#D37295", "#79706E")
  if (n > length(colors)) {
    grDevices::hcl.colors(n, palette = "Dark 3")
  } else {
    colors[seq_len(n)]
  }
}

.NGenesPerCell <- function(countData, gene_ids) {
  detected <- summary(countData > 0)
  n_genes <- integer(ncol(countData))
  if (nrow(detected) == 0) {
    return(n_genes)
  }

  genes_by_cell <- split(gene_ids[detected$i], detected$j)
  n_genes[as.integer(names(genes_by_cell))] <- lengths(lapply(genes_by_cell, unique))
  n_genes
}

.CellUsageDispersion <- function(expr_mat, gene_ids) {
  result <- data.frame(
    transcript = rownames(expr_mat),
    cell.n = integer(nrow(expr_mat)),
    cell.prop.mean = rep(NA_real_, nrow(expr_mat)),
    cell.prop.median = rep(NA_real_, nrow(expr_mat)),
    cell.prop.sd = rep(NA_real_, nrow(expr_mat)),
    cell.prop.iqr = rep(NA_real_, nrow(expr_mat)),
    cell.prop.zero.frac = rep(NA_real_, nrow(expr_mat)),
    check.names = FALSE
  )

  for (idx in split(seq_along(gene_ids), gene_ids)) {
    gene_mat <- expr_mat[idx, , drop = FALSE]
    gene_totals <- colSums(gene_mat)
    eligible <- which(gene_totals > 0)
    result$cell.n[idx] <- as.integer(length(eligible))

    if (length(eligible) == 0) {
      next
    }

    cell_props <- as.matrix(gene_mat[, eligible, drop = FALSE])
    cell_props <- sweep(cell_props, 2, gene_totals[eligible], "/")
    result$cell.prop.mean[idx] <- rowMeans(cell_props)
    result$cell.prop.median[idx] <- apply(cell_props, 1, stats::median)
    result$cell.prop.sd[idx] <- apply(cell_props, 1, stats::sd)
    result$cell.prop.iqr[idx] <- apply(cell_props, 1, stats::IQR)
    result$cell.prop.zero.frac[idx] <- rowMeans(cell_props == 0)
  }

  result
}

.CellDiversityDispersion <- function(expr_mat, gene_ids, div.func) {
  gene_levels <- unique(gene_ids)
  result <- data.frame(
    gene.id = gene_levels,
    cell.n = integer(length(gene_levels)),
    cell.div.mean = rep(NA_real_, length(gene_levels)),
    cell.div.median = rep(NA_real_, length(gene_levels)),
    cell.div.sd = rep(NA_real_, length(gene_levels)),
    cell.div.iqr = rep(NA_real_, length(gene_levels)),
    check.names = FALSE
  )

  gene_rows <- split(seq_along(gene_ids), gene_ids)
  for (gene in gene_levels) {
    idx <- gene_rows[[gene]]
    gene_mat <- expr_mat[idx, , drop = FALSE]
    gene_totals <- colSums(gene_mat)
    eligible <- which(gene_totals > 0)
    out_idx <- match(gene, result$gene.id)
    result$cell.n[out_idx] <- as.integer(length(eligible))

    if (length(eligible) == 0) {
      next
    }

    cell_counts <- as.matrix(gene_mat[, eligible, drop = FALSE])
    cell_props <- sweep(cell_counts, 2, gene_totals[eligible], "/")
    detected <- colSums(cell_counts > 0)
    cell_div <- vapply(seq_len(ncol(cell_props)), function(i) {
      value <- div.func(cell_props[, i])
      if (detected[i] == 1 && !is.finite(value)) 0 else value
    }, numeric(1))

    result$cell.div.mean[out_idx] <- mean(cell_div)
    result$cell.div.median[out_idx] <- stats::median(cell_div)
    result$cell.div.sd[out_idx] <- stats::sd(cell_div)
    result$cell.div.iqr[out_idx] <- stats::IQR(cell_div)
  }

  result
}

.CellExpressionDispersion <- function(expr_mat) {
  result <- data.frame(
    transcript = rownames(expr_mat),
    cell.n = rep(as.integer(ncol(expr_mat)), nrow(expr_mat)),
    cell.expr.median = rep(NA_real_, nrow(expr_mat)),
    cell.expr.sd = rep(NA_real_, nrow(expr_mat)),
    cell.expr.iqr = rep(NA_real_, nrow(expr_mat)),
    check.names = FALSE
  )

  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    return(result)
  }

  cell_expr <- suppressWarnings(as.matrix(expr_mat))
  result$cell.expr.median <- apply(cell_expr, 1, stats::median)
  result$cell.expr.sd <- apply(cell_expr, 1, stats::sd)
  result$cell.expr.iqr <- apply(cell_expr, 1, stats::IQR)
  result
}

.PAdjustMethod <- function(p.adj) {
  assertChoice(p.adj, stats::p.adjust.methods)
  p.adj
}
