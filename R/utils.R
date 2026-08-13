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

.DiversityFunction <- function(entropy.use, order = NULL) {
  force(entropy.use)
  force(order)

  function(x) {
    if (entropy.use == "Shannon") {
      x <- head(sort(x, decreasing = TRUE), 2)
      -sum(x[x > 0] * log(x[x > 0]))
    } else if (entropy.use == "NormalizedShannon") {
      n_x <- sum(x > 0)
      (-sum(x[x > 0] * log(x[x > 0]))) / log(n_x)
    } else if (entropy.use == "Renyi") {
      order.use <- order
      if (is.null(order.use)) order.use <- 2
      x <- head(sort(x, decreasing = TRUE), 2)
      (1 / (1 - order.use)) * log(sum((x[x > 0])^order.use))
    } else if (entropy.use == "NormalizedRenyi") {
      order.use <- order
      if (is.null(order.use)) order.use <- 2
      n_x <- sum(x > 0)
      (1 / (1 - order.use)) * log(sum((x[x > 0])^order.use)) / log(n_x)
    } else if (entropy.use == "GiniSimpson") {
      1 - sum((x[x > 0])^2)
    } else if (entropy.use == "Tsallis") {
      order.use <- order
      if (is.null(order.use)) order.use <- 3
      (1 - sum(x[x > 0]^order.use)) / (order.use - 1)
    } else if (entropy.use == "InverseSimpson") {
      1 / sum((x[x > 0])^2)
    }
  }
}

.DiversityThreshold <- function(entropy.use, entropy.thresh = NULL) {
  if (!is.null(entropy.thresh)) {
    return(entropy.thresh)
  }

  switch(
    entropy.use,
    "Shannon" = 0.500,
    "NormalizedShannon" = 0,
    "Renyi" = 0.435,
    "NormalizedRenyi" = 0,
    "GiniSimpson" = 0.348,
    "Tsallis" = 0.243,
    "InverseSimpson" = 1.533
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

.PAdjustMethod <- function(p.adj) {
  assertChoice(p.adj, stats::p.adjust.methods)
  p.adj
}
