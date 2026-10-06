# Prepare the GBM02PR, RCC12P, and HR182 objects

suppressPackageStartupMessages(library(Hypatia))
require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing input file(s): ", paste(missing, collapse = ", "))
}

validate_counts <- function(counts) {
  checkmate::assertMultiClass(counts, c("matrix", "dgCMatrix"))
  checkmate::assertCharacter(rownames(counts), any.missing = FALSE, unique = TRUE, min.len = 1)
  checkmate::assertCharacter(colnames(counts), any.missing = FALSE, unique = TRUE, min.len = 1)
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  checkmate::assertNumeric(values, lower = 0, any.missing = FALSE, finite = TRUE)
  if (any(values != floor(values))) stop("Raw counts must be integers.")
  invisible(counts)
}

load_application_dataset <- function(prefix, input_dir) {
  path <- file.path(input_dir, paste0("input_", prefix, ".rda"))
  require_files(path)
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  required <- paste0(prefix, c("_countData", "_rowData", "_colData"))
  if (!all(required %in% loaded)) {
    stop(path, " is missing: ", paste(setdiff(required, loaded), collapse = ", "))
  }
  counts <- env[[required[1]]]
  row_data <- as.data.frame(env[[required[2]]])
  col_data <- as.data.frame(env[[required[3]]])
  validate_counts(counts)
  if (!identical(rownames(counts), rownames(row_data)) ||
      !identical(colnames(counts), rownames(col_data))) {
    stop(path, ": metadata row names must match the count matrix in order.")
  }
  columns <- c("gene_name", "structural_category", "RTS_stage",
               "within_CAGE_peak", "within_polyA_site")
  if (!all(columns %in% names(row_data))) {
    stop(path, ": rowData is missing ", paste(setdiff(columns, names(row_data)), collapse = ", "))
  }
  checkmate::assertCharacter(row_data$gene_name, any.missing = FALSE, min.chars = 1)
  for (name in c("RTS_stage", "within_CAGE_peak", "within_polyA_site")) {
    checkmate::assertLogical(row_data[[name]])
  }
  list(counts = counts, row_data = row_data, col_data = col_data)
}

input_dir <- "manuscript/input"
output_dir <- "manuscript/figure2/data"
datasets <- c("gbm", "rcc", "heart")
inputs <- file.path(input_dir, paste0("input_", datasets, ".rda"))
require_files(inputs)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prepare_object <- function(prefix) {
  input <- load_application_dataset(prefix, input_dir)
  object <- Hypatia::CreateSCE(
    countData = input$counts, colData = input$col_data, rowData = input$row_data
  )
  object <- Hypatia::SetGenes(object, "gene_name")
  object <- Hypatia::SetTranscripts(object, "transcript_name")
  object <- Hypatia::SetGroups(object, "cell_type")
  object <- Hypatia::SubsetTranscripts(object, subset = nCell >= 3)
  object <- Hypatia::SubsetTranscripts(object, subset =
    (structural_category == "full-splice_match" & !RTS_stage) |
      (within_CAGE_peak & within_polyA_site & !RTS_stage))
  object <- Hypatia::SubsetCells(object, subset = nTranscript >= 50)
  Hypatia::NormalizeCounts(object, scale.factor = 10000)
}

objects <- setNames(lapply(datasets, prepare_object), datasets)
save(list = names(objects), envir = list2env(objects),
     file = file.path(output_dir, "hypatia_objects.rda"))
saveRDS(utils::sessionInfo(), file.path(output_dir, "prepare_session.rds"))
