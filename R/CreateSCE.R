#' Create a _SingleCellExperiment_ object
#'
#' Creates a `SingleCellExperiment` object from isoform counts, cell-level metadata, and isoform-level metadata.
#'
#' @param countData A numeric matrix containing raw single-cell isoform counts where row names are transcript IDs and column names are cell IDs.
#' @param colData A data frame of cell-level metadata. Row names must match the column names of `countData`.
#' @param rowData A data frame of isoform-level metadata. Row names must match the column names of `countData`.
#'  Must include a column specifying gene IDs associated with each isoform as defined by the `active.gene.id` parameter.
#' @param active.gene.id Name of the column in `rowData` containing associated gene IDs associated with each isoform.
#' Defaults to `"gene_id"`.
#' @param active.group.id Name of the column in `colData` that defines cell grouping used in downstream functions. If `NULL`, no default grouping is assigned to the object metadata.
#' @param active.transcript.id Name of the column in `rowData` containing alternative transcript IDs (transcript names, ENST IDs, etc.) to use instead of row names of `countData`. Values in this column must be unique.
#' @param gtf A `GRanges` object containing structural annotations for all transcripts in `countData`.
#' @param gtf.transcript.id Name of metadata column in `gtf` that corresponds to the transcript IDs in `countData`.
#' @param project Project name to be stored in the object's `colData` slot.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A `SingleCellExperiment` object.
#'
#' @export
#' @import checkmate
#' @import Matrix
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata metadata<-

CreateSCE <- function(
    countData,
    colData,
    rowData,
    active.gene.id = "gene_id",
    active.group.id = NULL,
    active.transcript.id = NULL,
    gtf = NULL,
    gtf.transcript.id = NULL,
    project = "Project",
    quiet = FALSE
) {

  # Check inputs
  if(!quiet) message("\nChecking inputs... ")

  ## countData
  assertMultiClass(countData, c("matrix", "dgCMatrix", "data.frame"))
  assertCharacter(rownames(countData))
  assertCharacter(colnames(countData))
  if(!quiet) message("\u2714 countData")

  ## colData
  assertClass(colData, "data.frame")
  assertTRUE(identical(rownames(colData), colnames(countData)))
  if (!any(vapply(colData, function(x) length(unique(na.omit(x))) >= 2, logical(1)))) {
    stop("colData must have at least one column with more than 2 unique values")
  }
  if (!is.null(active.group.id)) {
    assertChoice(active.group.id, colnames(colData), null.ok = TRUE)
    assertFALSE(anyMissing(colData[[active.group.id]]))
  }
  if(!quiet) message("\u2714 colData")

  ## rowData
  assertClass(rowData, "data.frame")
  assertTRUE(identical(rownames(rowData), rownames(countData)))
  assertChoice(active.gene.id, colnames(rowData))
  assertCharacter(rowData[[active.gene.id]], any.missing = FALSE)
  if (!is.null(active.transcript.id)) {
    assertChoice(active.transcript.id, colnames(rowData), null.ok = TRUE)
    assertCharacter(rowData[[active.transcript.id]], any.missing = FALSE)
    assertFALSE(any(duplicated(rowData[[active.transcript.id]])))
  }
  if(!quiet) message("\u2714 rowData")

  ## gtf
  assertClass(gtf, "GRanges", null.ok = TRUE)
  if (!is.null(gtf)) {
    assertTRUE(gtf.transcript.id %in% names(mcols(gtf)))
    assertTRUE(all(rownames(countData) %in% mcols(gtf)[[gtf.transcript.id]]))
    if (quiet) message("\u2714 gtf")
  }

  assertString(project)
  assertFlag(quiet)


  # Prepare countData
  if (!inherits(countData, "dgCMatrix")) {
    if(!quiet) message("Converting countData to sparse matrix format...")
    countData <- as.matrix(countData, nrow = nrow(countData), byrow = FALSE)
    countData <- as(countData, "dgCMatrix")
  }


  # Prepare colData
  colData[["project"]] <- project
  colData[["nCount"]] <- colSums(countData)
  colData[["nTranscript"]] <- diff(countData@p)

  ## nGene
  gene_ids <- rowData[[active.gene.id]]
  n_genes_per_cell <- vapply(seq_len(ncol(countData)), function(j) {
    idx <- (countData@p[j] + 1) : countData@p[j + 1]
    if (length(idx) == 0) return(0L)
    transcripts <- countData@i[idx] + 1
    length(unique(gene_ids[transcripts]))
  },
  integer(1)
  )
  colData[["nGene"]] <- n_genes_per_cell

  ## column order
  colData <- colData[, c("project", "nCount", "nTranscript", "nGene",
                                               setdiff(names(colData), c("project", "nCount", "nTranscript", "nGene")))]


  # Prepare rowData
  rowData$nCell <- rowSums(countData > 0)
  rowData <- rowData[, c("nCell", setdiff(names(rowData), "nCell"))]


  # Create the SCE object
  if(!quiet) message("Creating a SingleCellExperiment object of ", ncol(countData), " cells and ", nrow(countData), " transcripts...")
  object <- SingleCellExperiment(
    assays = list(counts = countData),
    colData = colData,
    rowData = rowData
  )


  # Metadata slots
  metadata(object)$active.gene.id <- active.gene.id
  if (is.null(active.group.id)) {
    metadata(object)$active.group.id <- ""
  } else {
    metadata(object)$active.group.id <- active.group.id
  }
  if (is.null(active.transcript.id)) {
    metadata(object)$active.transcript.id <- ""
  } else {
    metadata(object)$active.transcript.id <- active.transcript.id
  }


  # Prepare rowRanges
  if (!is.null(gtf)) {
    if(!quiet) message("Adding rowRanges to object...")

    ## store provided rowRanges in metadata
    metadata(object)$GTF <- gtf

    ## assign names of GRange object
    names(gtf) <- mcols(gtf)[[gtf.transcript.id]]

    ## filter and match with countData rownames
    gtf <- gtf[rownames(countData)]
    if (!identical(rownames(object), names(gtf))) {
      stop("Please check that gtf mcol '", gtf.transcript.id, "' contains all transcripts found in countData.")
    }

    ## add rowRanges to the object - do not copy mcols to object rowData
    mcols(gtf) <- NULL
    rowRanges(object) <- gtf

    ## add back original rowData
    rowData(object) <- rowData

  }

  if(!quiet) message("Done.")

  return(object)
}

