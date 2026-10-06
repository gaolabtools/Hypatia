# Fits distributions for the simulation to filtered-gene LR-scRNA-seq data from GBM, heart, and RCC

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(Matrix)
})
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

transcript_quality_pass <- function(row_data) {
  pass <- with(row_data,
    (structural_category == "full-splice_match" & !RTS_stage) |
      (within_CAGE_peak & within_polyA_site & !RTS_stage))
  pass[is.na(pass)] <- FALSE
  pass
}

seed <- 1029
set.seed(seed)
application_data_dir <- "manuscript/input"
output_data_dir <- "manuscript/figure1/data"
max_fit_genes <- 5000L
alpha_grid <- 10^seq(-2, 2, length.out = 81)
k_bin_levels <- c("2", "3", "4", "5", "6+")
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
fit_output_file <- file.path(output_data_dir, "simulation_fit.rds")
parameter_output_file <- file.path(output_data_dir, "simulation_fit_parameters.tsv")
dirichlet_grid_output_file <- file.path(output_data_dir, "simulation_fit_dirichlet_grid.tsv")
real_datasets <- c(GBM = "gbm", Heart = "heart", RCC = "rcc")
require_files(file.path(application_data_dir, paste0("input_", real_datasets, ".rda")))
scope_levels <- c(names(real_datasets), "Combined")

# Data preparation ---------------------------------------------------------

prepare_fit_data <- function(dataset, prefix) {
  input <- load_application_dataset(prefix, application_data_dir)
  counts <- input$counts
  row_data <- input$row_data

  gene_ids <- as.character(row_data$gene_name)
  detected_cells <- Matrix::rowSums(counts > 0)
  quality_pass <- transcript_quality_pass(row_data)
  keep_transcripts <- detected_cells >= 3 & quality_pass

  counts <- counts[keep_transcripts, , drop = FALSE]
  gene_ids <- gene_ids[keep_transcripts]
  keep_cells <- Matrix::colSums(counts > 0) >= 50
  counts <- counts[, keep_cells, drop = FALSE]

  # Estimate cell depth from all transcripts that pass QC, before selecting genes with multiple isoforms.
  cell_total_count <- Matrix::colSums(counts)
  median_cell_total <- stats::median(cell_total_count)
  if (!is.finite(median_cell_total) || median_cell_total <= 0) {
    stop(dataset, " has a non-positive median cell total after filtering.")
  }
  cell_size_factors <- cell_total_count / median_cell_total

  # Remove transcripts with no remaining counts, then keep genes with at least two observed isoforms.
  keep_observed <- Matrix::rowSums(counts) > 0
  counts <- counts[keep_observed, , drop = FALSE]
  gene_ids <- gene_ids[keep_observed]
  observed_transcripts_per_gene <- table(gene_ids)
  multi_isoform_genes <- names(observed_transcripts_per_gene)[
    observed_transcripts_per_gene >= 2L
  ]
  keep_multi_isoform <- gene_ids %in% multi_isoform_genes
  counts <- counts[keep_multi_isoform, , drop = FALSE]
  gene_ids <- gene_ids[keep_multi_isoform]

  transcript_totals <- Matrix::rowSums(counts)
  transcript_split <- split(transcript_totals, gene_ids)
  gene_totals <- vapply(transcript_split, sum, numeric(1))
  n_cells <- ncol(counts)

  if (length(transcript_split) < 2L) {
    stop(dataset, " needs at least two observed multi-isoform genes after filtering.")
  }
  profiles <- tibble(
    dataset = dataset,
    gene = names(transcript_split),
    n_transcripts = lengths(transcript_split),
    gene_total_count = unname(gene_totals),
    gene_mean_count = unname(gene_totals) / n_cells,
    proportions = lapply(
      transcript_split,
      function(x) as.numeric(x / sum(x))
    )
  ) %>%
    mutate(
      dominant_proportion = map_dbl(proportions, max),
      shannon_entropy = map_dbl(
        proportions,
        function(p) -sum(p * log(p))
      ),
      k_bin = factor(
        if_else(n_transcripts >= 6L, "6+", as.character(n_transcripts)),
        levels = k_bin_levels
      )
    )

  cells <- tibble(
    dataset = dataset,
    cell_total_count = as.numeric(cell_total_count),
    cell_size_factor = as.numeric(cell_size_factors)
  )

  dataset_summary <- tibble(
    dataset = dataset,
    n_cells = n_cells,
    n_genes = nrow(profiles),
    n_transcripts = nrow(counts),
    median_cell_total = median_cell_total
  )

  list(
    profiles = profiles,
    cells = cells,
    dataset_summary = dataset_summary
  )
}

# Parametric fitting -------------------------------------------------------

fit_shifted_negative_binomial <- function(k, scope, offset = 2L) {
  x <- as.integer(k) - offset
  if (any(x < 0L)) {
    stop("Transcript counts must be at least the negative-binomial offset.")
  }

  sample_mean <- mean(x)
  sample_variance <- stats::var(x)
  start_size <- if (sample_variance > sample_mean && sample_mean > 0) {
    sample_mean^2 / (sample_variance - sample_mean)
  } else {
    100
  }
  start_mu <- max(sample_mean, 1e-4)

  fit <- stats::optim(
    par = log(c(size = start_size, mu = start_mu)),
    fn = function(log_parameters) {
      parameters <- exp(log_parameters)
      -sum(stats::dnbinom(
        x,
        size = parameters[["size"]],
        mu = parameters[["mu"]],
        log = TRUE
      ))
    },
    method = "BFGS"
  )
  estimates <- exp(fit$par)

  tibble(
    scope = scope,
    n = length(k),
    offset = offset,
    size = unname(estimates[["size"]]),
    mu = unname(estimates[["mu"]]),
    fitted_mean_transcripts = offset + unname(estimates[["mu"]]),
    negative_log_likelihood = fit$value,
    convergence = fit$convergence
  )
}

fit_lognormal <- function(x, scope, component, fixed_meanlog = NULL) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop(component, " values must be finite and positive.")
  }
  log_x <- log(x)
  meanlog <- if (is.null(fixed_meanlog)) mean(log_x) else fixed_meanlog
  sdlog <- sqrt(mean((log_x - meanlog)^2))

  tibble(
    component = component,
    scope = scope,
    n = length(x),
    meanlog = meanlog,
    sdlog = sdlog,
    fitted_median = exp(meanlog),
    fitted_mean = exp(meanlog + sdlog^2 / 2)
  )
}

simulate_dirichlet_metrics <- function(k, alpha) {
  dominant <- numeric(length(k))
  entropy <- numeric(length(k))

  for (current_k in sort(unique(k))) {
    index <- which(k == current_k)
    gamma_draws <- matrix(
      stats::rgamma(length(index) * current_k, shape = alpha),
      nrow = length(index),
      ncol = current_k
    )
    gamma_draws <- pmax(gamma_draws, .Machine$double.xmin)
    proportions <- gamma_draws / rowSums(gamma_draws)
    dominant[index] <- apply(proportions, 1, max)
    entropy[index] <- -rowSums(proportions * log(proportions))
  }

  list(dominant = dominant, entropy = entropy)
}

ks_distance <- function(x, y) {
  x <- sort(x[is.finite(x)])
  y <- sort(y[is.finite(y)])
  evaluation_points <- sort(unique(c(x, y)))
  max(abs(
    findInterval(evaluation_points, x) / length(x) -
      findInterval(evaluation_points, y) / length(y)
  ))
}

fit_dirichlet_alpha <- function(data, scope, k_bin, fit_index) {
  if (nrow(data) > max_fit_genes) {
    data <- data[sample.int(nrow(data), max_fit_genes), , drop = FALSE]
  }

  real_dominant <- data$dominant_proportion
  real_entropy <- data$shannon_entropy
  k <- data$n_transcripts

  grid_fit <- map_dfr(seq_along(alpha_grid), function(alpha_index) {
    alpha <- alpha_grid[[alpha_index]]
    set.seed(seed + fit_index * 1000L + alpha_index)
    simulated <- simulate_dirichlet_metrics(k, alpha)
    dominant_distance <- ks_distance(real_dominant, simulated$dominant)
    entropy_distance <- ks_distance(real_entropy, simulated$entropy)

    tibble(
      scope = scope,
      k_bin = k_bin,
      n = nrow(data),
      alpha = alpha,
      dominant_ks = dominant_distance,
      entropy_ks = entropy_distance,
      objective = mean(c(dominant_distance, entropy_distance))
    )
  })

  best <- grid_fit %>%
    arrange(objective, dominant_ks, entropy_ks, alpha) %>%
    slice(1)

  list(grid = grid_fit, best = best)
}

# Fit each dataset and the pooled data -------------------------------------

prepared <- imap(real_datasets, function(prefix, dataset) {
  message("Preparing ", dataset, ".")
  result <- prepare_fit_data(dataset, prefix)
  invisible(gc())
  result
})

profiles <- map_dfr(prepared, "profiles")
cells <- map_dfr(prepared, "cells")
dataset_summary <- map_dfr(prepared, "dataset_summary")

profile_scopes <- c(
  split(profiles, profiles$dataset),
  list(Combined = profiles)
)
cell_scopes <- c(
  split(cells, cells$dataset),
  list(Combined = cells)
)

transcript_number_fit <- imap_dfr(profile_scopes, function(data, scope) {
  fit_shifted_negative_binomial(data$n_transcripts, scope)
}) %>%
  mutate(scope = factor(scope, levels = scope_levels)) %>%
  arrange(scope)

gene_abundance_fit <- imap_dfr(profile_scopes, function(data, scope) {
  fit_lognormal(data$gene_mean_count, scope, "gene_abundance")
}) %>%
  mutate(scope = factor(scope, levels = scope_levels)) %>%
  arrange(scope)

cell_size_fit <- imap_dfr(cell_scopes, function(data, scope) {
  # Fix the fitted median at one because cell size factors are scaled to each dataset's median.
  fit_lognormal(
    data$cell_size_factor,
    scope,
    "cell_size_factor",
    fixed_meanlog = 0
  )
}) %>%
  mutate(scope = factor(scope, levels = scope_levels)) %>%
  arrange(scope)

dirichlet_results <- list()
fit_index <- 0L
for (scope in names(profile_scopes)) {
  scope_data <- profile_scopes[[scope]]
  for (k_bin in k_bin_levels) {
    fit_index <- fit_index + 1L
    bin_data <- scope_data %>%
      filter(as.character(.data$k_bin) == .env$k_bin)
    if (
      nrow(bin_data) > 0L &&
        !all(as.character(bin_data$k_bin) == k_bin)
    ) {
      stop("Dirichlet fitting data contain an unexpected K bin.")
    }
    if (nrow(bin_data) < 50L) {
      warning("Skipping ", scope, " K=", k_bin, ": fewer than 50 genes.")
      next
    }
    message("Fitting Dirichlet alpha for ", scope, ", K=", k_bin, ".")
    dirichlet_results[[paste(scope, k_bin, sep = "_")]] <-
      fit_dirichlet_alpha(bin_data, scope, k_bin, fit_index)
  }
}

if (!length(dirichlet_results)) stop("No K bin has at least 50 genes for Dirichlet fitting.")

dirichlet_grid <- map_dfr(dirichlet_results, "grid") %>%
  mutate(
    scope = factor(scope, levels = scope_levels),
    k_bin = factor(k_bin, levels = k_bin_levels)
  ) %>%
  arrange(scope, k_bin, alpha)

dirichlet_fit <- map_dfr(dirichlet_results, "best") %>%
  mutate(
    scope = factor(scope, levels = scope_levels),
    k_bin = factor(k_bin, levels = k_bin_levels)
  ) %>%
  arrange(scope, k_bin)

# Save the fitted parameters -----------------------------------------------

parameter_table <- bind_rows(
  transcript_number_fit %>%
    select(scope, n, offset, size, mu, fitted_mean_transcripts) %>%
    pivot_longer(
      cols = c(offset, size, mu, fitted_mean_transcripts),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(component = "observed_transcript_number", stratum = "K >= 2"),
  gene_abundance_fit %>%
    select(scope, n, meanlog, sdlog, fitted_median, fitted_mean) %>%
    pivot_longer(
      cols = c(meanlog, sdlog, fitted_median, fitted_mean),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(component = "gene_abundance", stratum = "All multi-isoform genes"),
  cell_size_fit %>%
    select(scope, n, meanlog, sdlog, fitted_median, fitted_mean) %>%
    pivot_longer(
      cols = c(meanlog, sdlog, fitted_median, fitted_mean),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(component = "cell_size_factor", stratum = "All filtered cells"),
  dirichlet_fit %>%
    select(scope, n, k_bin, alpha, dominant_ks, entropy_ks, objective) %>%
    pivot_longer(
      cols = c(alpha, dominant_ks, entropy_ks, objective),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(component = "isoform_proportions", stratum = paste0("K = ", k_bin)) %>%
    select(-k_bin)
) %>%
  select(component, scope, stratum, parameter, value, n) %>%
  arrange(component, scope, stratum, parameter)

fit_output <- list(
  seed = seed,
  filters = list(
    minimum_detected_cells_per_transcript = 3L,
    minimum_detected_transcripts_per_cell = 50L,
    minimum_observed_transcripts_per_gene = 2L,
    real_data_sqanti_filter = TRUE
  ),
  dataset_summary = dataset_summary,
  transcript_number = transcript_number_fit,
  gene_abundance = gene_abundance_fit,
  cell_size_factor = cell_size_fit,
  isoform_proportions = dirichlet_fit,
  dirichlet_grid = dirichlet_grid
)

saveRDS(fit_output, fit_output_file, compress = "gzip")
readr::write_tsv(parameter_table, parameter_output_file)
readr::write_tsv(dirichlet_grid, dirichlet_grid_output_file)

message("Wrote fitted objects to ", fit_output_file)
message("Wrote parameter table to ", parameter_output_file)
message("Wrote Dirichlet grid to ", dirichlet_grid_output_file)

saveRDS(utils::sessionInfo(), file.path(output_data_dir, "simulation_fit_session.rds"))
