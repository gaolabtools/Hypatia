# Run Hypatia and satuRn on the simulated counts

suppressPackageStartupMessages({
  library(Hypatia)
  library(satuRn)
})

validate_counts <- function(counts) {
  checkmate::assertMultiClass(counts, c("matrix", "dgCMatrix"))
  checkmate::assertCharacter(rownames(counts), any.missing = FALSE, unique = TRUE, min.len = 1)
  checkmate::assertCharacter(colnames(counts), any.missing = FALSE, unique = TRUE, min.len = 1)
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  checkmate::assertNumeric(values, lower = 0, any.missing = FALSE, finite = TRUE)
  if (any(values != floor(values))) stop("Raw counts must be integers.")
  invisible(counts)
}

# Check transcript IDs and cell group labels before running either tool.
validate_simulation <- function(sim_mat) {
  validate_counts(sim_mat)
  if (any(!grepl("^test[0-9]+[.]cat[0-9]+$", rownames(sim_mat)))) {
    stop("Expected transcript IDs test<integer>.cat<integer>.")
  }
  if (!setequal(sub("_.*$", "", colnames(sim_mat)), c("G1", "G2"))) {
    stop("Expected two simulated groups, with cell IDs G1_* and G2_*.")
  }
  invisible(sim_mat)
}

run_hypatia <- function(sim_mat) {
  set.seed(1029)
  object <- Hypatia::CreateSCE(
    countData = sim_mat,
    colData = data.frame(group = sub("_.*$", "", colnames(sim_mat)),
                         row.names = colnames(sim_mat)),
    rowData = data.frame(sim.gene.id = sub("[.]cat[0-9]+$", "", rownames(sim_mat)),
                         row.names = rownames(sim_mat)),
    active.gene.id = "sim.gene.id",
    active.group.id = "group"
  )
  result <- Hypatia::RunDIU(
    object, min.gene.pct = 0, min.gene.cts = 0, min.tx.cts = 1,
    bootstrap = TRUE, boot.iter = 5000,
    permutation = TRUE, perm.iter = 10000
  )
  # Keep gene-level bootstrap summaries and drop the larger per-transcript draws.
  result$data <- result$data[, setdiff(names(result$data),
    c("boot.prop.1", "boot.prop.2", "boot.prop.diff")), drop = FALSE]
  result
}

# Fit and test isoform usage with satuRn.
run_saturn <- function(sim_mat, run_id, input_file) {
  cores <- parallel::detectCores(logical = FALSE)
  workers <- if (is.na(cores)) 1L else min(20L, cores)
  set.seed(1029)
  transcript_id <- rownames(sim_mat)
  gene_id <- sub("[.]cat[0-9]+$", "", transcript_id)
  if (any(gene_id == transcript_id)) {
    stop("Could not derive every gene ID from the simulation transcript IDs.")
  }

  transcript_total <- as.numeric(Matrix::rowSums(sim_mat))
  nonzero_isoforms_per_gene <- table(gene_id[transcript_total > 0])
  genes_with_multiple_isoforms <- names(nonzero_isoforms_per_gene)[
    nonzero_isoforms_per_gene >= 2L
  ]

  # Keep detected transcripts from genes with at least two detected isoforms.
  keep <- transcript_total >= 1 &
    gene_id %in% genes_with_multiple_isoforms

  feature_filter <- data.frame(
    transcript_id = transcript_id,
    gene_id = gene_id,
    total_count = transcript_total,
    selected_by_limit = TRUE, # All simulated genes are included.
    kept = keep,
    filter_reason = ifelse(
      transcript_total < 1, "zero_total_count",
      ifelse(!gene_id %in% genes_with_multiple_isoforms,
             "fewer_than_two_detected_isoforms", "kept")
    ),
    stringsAsFactors = FALSE
  )

  sim_mat <- sim_mat[keep, , drop = FALSE]
  gene_id <- gene_id[keep]

  if (nrow(sim_mat) == 0L) {
    stop("No transcripts remain after filtering.")
  }

  group_label <- sub("_.*$", "", colnames(sim_mat))
  group_levels <- unique(group_label)
  if (length(group_levels) != 2L) {
    stop("satuRn benchmark expects exactly two simulation groups.")
  }
  group <- factor(group_label, levels = group_levels)

  row_data <- S4Vectors::DataFrame(
    isoform_id = rownames(sim_mat),
    gene_id = gene_id,
    row.names = rownames(sim_mat)
  )
  col_data <- S4Vectors::DataFrame(
    group = group,
    row.names = colnames(sim_mat)
  )

  saturn_object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = sim_mat),
    rowData = row_data,
    colData = col_data
  )
  S4Vectors::metadata(saturn_object)$formula <- ~ 0 + group

  design <- stats::model.matrix(~ 0 + group, data = as.data.frame(col_data))
  contrast_name <- paste0(group_levels[2], "_vs_", group_levels[1])
  contrast <- matrix(
    c(-1, 1),
    ncol = 1,
    dimnames = list(colnames(design), contrast_name)
  )

  if (workers > 1L && .Platform$OS.type != "windows") {
    bpparam <- BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  } else if (workers > 1L) {
    bpparam <- BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = TRUE)
  }

  message(
    "Fitting satuRn to ", nrow(sim_mat), " transcripts from ",
    length(unique(gene_id)), " genes and ", ncol(sim_mat), " cells using ",
    workers, " worker(s)."
  )

  fit_time <- system.time({
    saturn_object <- satuRn::fitDTU(
      object = saturn_object,
      formula = ~ 0 + group,
      parallel = workers > 1L,
      BPPARAM = bpparam,
      verbose = TRUE
    )
  })

  test_time <- system.time({
    saturn_object <- satuRn::testDTU(
      object = saturn_object,
      contrasts = contrast,
      diagplot1 = FALSE,
      diagplot2 = FALSE,
      sort = FALSE
    )
  })

  result_name <- paste0("fitDTUResult_", contrast_name)
  saturn_results <- as.data.frame(
    SummarizedExperiment::rowData(saturn_object)[[result_name]]
  )
  result_transcript_id <- rownames(saturn_results)
  result_order <- match(result_transcript_id, rownames(sim_mat))
  if (anyNA(result_order) || anyDuplicated(result_transcript_id)) {
    stop("Could not align satuRn test results to the filtered count matrix.")
  }
  saturn_results <- data.frame(
    transcript_id = result_transcript_id,
    gene_id = gene_id[result_order],
    total_count = transcript_total[match(
      result_transcript_id,
      transcript_id
    )],
    saturn_results,
    row.names = NULL,
    check.names = FALSE
  )

  benchmark_result <- list(
    tool = "satuRn",
    tool_version = as.character(utils::packageVersion("satuRn")),
    run_id = run_id,
    input_file = input_file,
    contrast = contrast_name,
    config = list(
      max_genes = 0L,
      workers = workers,
      minimum_transcript_total = 1L,
      minimum_detected_isoforms_per_gene = 2L,
      filter_strategy = "benchmark eligibility filter",
      formula = "~ 0 + group"
    ),
    dimensions = list(
      input_transcripts = nrow(feature_filter),
      input_genes = length(unique(feature_filter$gene_id)),
      tested_transcripts = nrow(sim_mat),
      tested_genes = length(unique(gene_id)),
      cells = ncol(sim_mat)
    ),
    timing = list(
      fit = fit_time,
      test = test_time
    ),
    feature_filter = feature_filter,
    results = saturn_results,
    session_info = utils::sessionInfo()
  )

  benchmark_result
}

run_id <- "run_3"
data_dir <- "manuscript/figure1/data"
input_file <- file.path(data_dir, paste0(run_id, "_sc_mat.rds"))
hypatia_file <- file.path(data_dir, paste0(run_id, "_hypatia_res.rds"))
saturn_file <- file.path(data_dir, paste0(run_id, "_saturn_res.rds"))
if (!file.exists(input_file)) stop("Missing simulation counts: ", input_file)
sim_mat <- readRDS(input_file)
validate_simulation(sim_mat)

hypatia_results <- run_hypatia(sim_mat)
saveRDS(hypatia_results, hypatia_file, compress = "gzip")
rm(hypatia_results)

saturn_results <- run_saturn(sim_mat, run_id, input_file)
saveRDS(saturn_results, saturn_file, compress = "gzip")
saveRDS(utils::sessionInfo(), file.path(data_dir, paste0(run_id, "_benchmark_run_session.rds")))
