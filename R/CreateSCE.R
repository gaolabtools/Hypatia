#' Create a SingleCellExperiment object
#'
#' Creates a `SingleCellExperiment` object from isoform counts, cell-level metadata, and transcript-level metadata.
#'
#' @param countData Numeric matrix, sparse matrix, or data frame of raw isoform counts. Row names are transcript IDs and column names are cell IDs.
#' @param colData A data frame of cell-level metadata. Row names must match the column names of `countData`.
#' @param rowData A data frame of transcript-level metadata. Row names must match the row names of `countData`. Must include the column named by `active.gene.id`.
#' @param active.gene.id Name of the `rowData` column containing gene IDs for each transcript.
#' @param active.group.id Optional `colData` column name to use as the default cell grouping in downstream functions.
#' @param active.transcript.id Optional `rowData` column name containing unique transcript IDs to report instead of row names.
#' @param gtf A `GRanges` object containing structural annotations for all transcripts in `countData`.
#' @param gtf.transcript.id Name of metadata column in `gtf` that corresponds to the transcript IDs in `countData`.
#' @param project Project name to be stored in the object's `colData` slot.
#' @param quiet Logical; if `TRUE`, suppresses messages.
#'
#' @returns A `SingleCellExperiment` object with a `"counts"` assay, QC columns in `colData`/`rowData`, and Hypatia settings in `metadata(object)`.
#'
#' @export
#' @import checkmate
#' @import Matrix
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#' @importFrom stats na.omit
#' @importFrom methods as

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
  if (!quiet) message("\nChecking inputs... ")

  ## countData
  assertMultiClass(countData, c("matrix", "dgCMatrix", "data.frame"))
  assertCharacter(rownames(countData), any.missing = FALSE, unique = TRUE)
  assertCharacter(colnames(countData), any.missing = FALSE, unique = TRUE)
  if (inherits(countData, "dgCMatrix")) {
    assertNumeric(countData@x, lower = 0, any.missing = FALSE)
  } else if (is.data.frame(countData)) {
    if (!all(vapply(countData, is.numeric, logical(1)))) {
      stop("countData must contain only numeric columns.")
    }
    assertNumeric(unlist(countData, use.names = FALSE), lower = 0, any.missing = FALSE)
  } else {
    assertNumeric(as.vector(countData), lower = 0, any.missing = FALSE)
  }
  if (!quiet) message("\u2714 countData")

  ## colData
  assertClass(colData, "data.frame")
  assertTRUE(identical(rownames(colData), colnames(countData)))
  if (!any(vapply(colData, function(x) length(unique(na.omit(x))) >= 2, logical(1)))) {
    stop("colData must have at least one column with at least 2 unique non-missing values.")
  }
  if (!is.null(active.group.id)) {
    assertChoice(active.group.id, colnames(colData), null.ok = TRUE)
    assertFALSE(anyMissing(colData[[active.group.id]]))
  }
  if (!quiet) message("\u2714 colData")

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
  if (!quiet) message("\u2714 rowData")

  ## gtf
  assertClass(gtf, "GRanges", null.ok = TRUE)
  if (!is.null(gtf)) {
    assertString(gtf.transcript.id, null.ok = FALSE)
    assertTRUE(gtf.transcript.id %in% names(mcols(gtf)))
    assertTRUE(all(rownames(countData) %in% mcols(gtf)[[gtf.transcript.id]]))
    if (!quiet) message("\u2714 gtf")
  }

  assertString(project)
  assertFlag(quiet)

  # Prepare countData
  if (!inherits(countData, "dgCMatrix")) {
    if (!quiet) message("Converting countData to sparse matrix format...")
    countData <- as.matrix(countData)
    countData <- as(countData, "dgCMatrix")
  }
  countData <- drop0(countData)

  # Prepare colData
  colData[["project"]] <- project
  colData[["nCount"]] <- colSums(countData)
  colData[["nTranscript"]] <- as.integer(colSums(countData > 0))

  ## nGene
  gene_ids <- rowData[[active.gene.id]]
  colData[["nGene"]] <- .NGenesPerCell(countData, gene_ids)

  ## column order
  qc_cols <- c("project", "nCount", "nTranscript", "nGene")
  colData <- colData[, c(qc_cols, setdiff(names(colData), qc_cols))]

  # Prepare rowData
  rowData$nCell <- rowSums(countData > 0)
  rowData <- rowData[, c("nCell", setdiff(names(rowData), "nCell"))]

  # Create the SCE object
  if (!quiet) {
    message(
      "Creating a SingleCellExperiment object of ",
      ncol(countData),
      " cells and ",
      nrow(countData),
      " transcripts..."
    )
  }
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
    if (!quiet) message("Adding rowRanges to object...")

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

  if (!quiet) message("Done.")

  return(object)
}
