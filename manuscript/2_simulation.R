# Simulate two groups from the fitted abundance, transcript-number, and isoform-proportion distributions

suppressPackageStartupMessages({
  library(tibble)
  library(purrr)
  library(Matrix)
  library(furrr)
  library(logger)
})
run_start <- Sys.time()
run_id <- "run_3"
n_cores <- 8L
data_dir <- "manuscript/figure1/data"
seed <- 1029
set.seed(seed)
n_tests <- 8000L
cell_number <- 500L
simulation_fit_file <- file.path(data_dir, "simulation_fit.rds")
simulation_fit_scope <- "Combined"
sc_disp <- 0.1 # negative-binomial size for cell-level gene counts
timeout <- 60  # time limit for generating each gene (seconds)
null_fraction <- 0.25
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(simulation_fit_file)) {
  stop(
    "Simulation fit not found: ", simulation_fit_file,
    ". Run manuscript/1_simulation-fit.R first."
  )
}

simulation_fit <- readRDS(simulation_fit_file)
required_fit_components <- c(
  "dataset_summary",
  "transcript_number",
  "gene_abundance",
  "cell_size_factor",
  "isoform_proportions"
)
missing_fit_components <- setdiff(required_fit_components, names(simulation_fit))
if (length(missing_fit_components) > 0L) {
  stop(
    "Simulation fit is missing: ",
    paste(missing_fit_components, collapse = ", ")
  )
}

select_fit_scope <- function(component, component_name, required_columns) {
  missing_columns <- setdiff(c("scope", required_columns), names(component))
  if (length(missing_columns) > 0L) {
    stop(
      component_name, " fit is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  selected <- component[
    as.character(component$scope) == simulation_fit_scope,
    ,
    drop = FALSE
  ]
  if (nrow(selected) != 1L) {
    stop(
      "Expected exactly one ", component_name, " fit for scope '",
      simulation_fit_scope, "'; found ", nrow(selected), "."
    )
  }
  selected
}

transcript_number_fit <- select_fit_scope(
  simulation_fit$transcript_number,
  "transcript-number",
  c("offset", "size", "mu")
)
gene_abundance_fit <- select_fit_scope(
  simulation_fit$gene_abundance,
  "gene-abundance",
  c("meanlog", "sdlog")
)
cell_size_fit <- select_fit_scope(
  simulation_fit$cell_size_factor,
  "cell-size",
  c("meanlog", "sdlog")
)

isoform_fit <- simulation_fit$isoform_proportions
missing_isoform_columns <- setdiff(
  c("scope", "k_bin", "alpha"),
  names(isoform_fit)
)
if (length(missing_isoform_columns) > 0L) {
  stop(
    "Isoform-proportion fit is missing: ",
    paste(missing_isoform_columns, collapse = ", ")
  )
}
isoform_fit <- isoform_fit[
  as.character(isoform_fit$scope) == simulation_fit_scope,
  ,
  drop = FALSE
]
required_k_bins <- c("2", "3", "4", "5", "6+")
if (
  nrow(isoform_fit) != length(required_k_bins) ||
    anyDuplicated(as.character(isoform_fit$k_bin)) ||
    !setequal(as.character(isoform_fit$k_bin), required_k_bins)
) {
  stop(
    "Isoform-proportion fit for scope '", simulation_fit_scope,
    "' must contain one row for each of: ",
    paste(required_k_bins, collapse = ", "), "."
  )
}
isoform_fit <- isoform_fit[
  match(required_k_bins, as.character(isoform_fit$k_bin)),
  ,
  drop = FALSE
]

transcript_offset <- as.integer(transcript_number_fit$offset[[1]])
transcript_size <- transcript_number_fit$size[[1]]
transcript_mu <- transcript_number_fit$mu[[1]]
fitted_gene_expr_meanlog <- gene_abundance_fit$meanlog[[1]]
gene_expr_sdlog <- gene_abundance_fit$sdlog[[1]]
cell_size_meanlog <- cell_size_fit$meanlog[[1]]
cell_size_sdlog <- cell_size_fit$sdlog[[1]]
isoform_alpha <- setNames(
  isoform_fit$alpha,
  as.character(isoform_fit$k_bin)
)

dataset_summary <- simulation_fit$dataset_summary
missing_summary_columns <- setdiff(
  c("dataset", "n_cells", "n_genes"),
  names(dataset_summary)
)
if (length(missing_summary_columns) > 0L) {
  stop(
    "Simulation-fit dataset summary is missing: ",
    paste(missing_summary_columns, collapse = ", ")
  )
}
if (
  any(!is.finite(dataset_summary$n_cells)) ||
    any(dataset_summary$n_cells <= 0) ||
    any(!is.finite(dataset_summary$n_genes)) ||
    any(dataset_summary$n_genes <= 0)
) {
  stop("Simulation-fit dataset cell and gene counts must be positive.")
}
# Use the gene-weighted geometric mean cell count to scale the pooled abundance fit.
reference_cell_number <- exp(stats::weighted.mean(
  log(dataset_summary$n_cells), dataset_summary$n_genes
))

simulated_cell_number <- 2L * cell_number
gene_count_scale <- reference_cell_number / simulated_cell_number
gene_expr_meanlog <- fitted_gene_expr_meanlog + log(gene_count_scale)

if (
  transcript_offset < 2L ||
    !isTRUE(all.equal(transcript_offset, transcript_number_fit$offset[[1]])) ||
    !is.finite(transcript_size) || transcript_size <= 0 ||
    !is.finite(transcript_mu) || transcript_mu < 0
) {
  stop("Invalid shifted negative-binomial transcript-number parameters.")
}
if (
  !is.finite(fitted_gene_expr_meanlog) ||
    !is.finite(gene_expr_meanlog) ||
    !is.finite(gene_expr_sdlog) || gene_expr_sdlog < 0
) {
  stop("Invalid gene-abundance log-normal parameters.")
}
if (
  !is.finite(cell_size_meanlog) ||
    !is.finite(cell_size_sdlog) || cell_size_sdlog < 0
) {
  stop("Invalid cell-size log-normal parameters.")
}
if (any(!is.finite(isoform_alpha)) || any(isoform_alpha <= 0)) {
  stop("All fitted Dirichlet alpha parameters must be finite and positive.")
}

get_isoform_alpha <- function(n_transcripts) {
  k_bin <- if (n_transcripts >= 6L) "6+" else as.character(n_transcripts)
  unname(isoform_alpha[[k_bin]])
}

simulate_isoform_proportions <- function(n_transcripts) {
  alpha <- get_isoform_alpha(n_transcripts)
  gamma_draws <- pmax(
    rgamma(n_transcripts, shape = alpha),
    .Machine$double.xmin
  )
  gamma_draws / sum(gamma_draws)
}

logfile <- file.path(data_dir, paste0(run_id, ".log"))
invisible(file.create(logfile))
log_appender(appender_tee(logfile))
log_info("--- Run details ---",
         "\nRun name: ", run_id,
         "\nData directory: ", data_dir,
         "\nCores: ", n_cores,
         "\nSeed: ", seed,
         "\nNumber of tests: ", n_tests,
         "\nCells per group: ", cell_number,
         "\nSimulation fit: ", simulation_fit_file,
         "\nSimulation fit scope: ", simulation_fit_scope,
         "\nTranscript-number offset: ", transcript_offset,
         "\nTranscript-number size: ", transcript_size,
         "\nTranscript-number mu: ", transcript_mu,
         "\nFitted gene expression meanlog: ", fitted_gene_expr_meanlog,
         "\nReference real cells: ", signif(reference_cell_number, 6),
         "\nTotal simulated cells: ", simulated_cell_number,
         "\nGene total-count scale: ", signif(gene_count_scale, 6),
         "\nScaled gene expression meanlog: ", gene_expr_meanlog,
         "\nGene expression sdlog: ", gene_expr_sdlog,
         "\nCell-size meanlog: ", cell_size_meanlog,
         "\nCell-size sdlog: ", cell_size_sdlog,
         "\nDirichlet alpha by K: ",
         paste(names(isoform_alpha), signif(isoform_alpha, 4), sep = "=", collapse = ", "),
         "\nDispersion parameter: ", sc_disp,
         "\nExact-null fraction: ", null_fraction,
         "\nTimeout (s): ", timeout)

# Use separate R sessions on Windows and forked workers on Unix.
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = n_cores)
} else {
  future::plan(future::multicore, workers = n_cores)
}


# Ground truth ------------------------------------------------------------
log_info("Simulating ground truth deltas and isoform proportion sets for ", n_tests, " genes...")

# Maximum usage change per gene
max_deltas <- rbeta(n = n_tests, shape1 = 0.5, shape2 = 10)

# Null genes share isoform proportions across groups, but their total expression can differ.
stopifnot(null_fraction >= 0, null_fraction <= 1)
n_exact_null <- round(n_tests * null_fraction)
exact_null <- rep(FALSE, n_tests)
if (n_exact_null > 0) {
  exact_null[sample.int(n_tests, n_exact_null)] <- TRUE
}
max_deltas[exact_null] <- 0

# Draw at least two isoforms per gene from the fitted shifted negative binomial.
cats <- transcript_offset + rnbinom(
  n = n_tests,
  size = transcript_size,
  mu = transcript_mu
)
log_info(
  "Simulated transcript categories per gene:",
  "\nMean: ", signif(mean(cats), 4),
  ", median: ", median(cats),
  ", maximum: ", max(cats)
)

# Usage changes and isoform proportions
gt_params <- tibble(
  testID = seq_along(max_deltas),
  max_delta = max_deltas,
  n_cats = cats,
  exact_null = exact_null
  )

gt_data <- future_pmap(
  .options = furrr_options(stdout = FALSE, seed = TRUE),
  .progress = TRUE,
  .l = list(
    gt_params$testID,
    gt_params$max_delta,
    gt_params$n_cats,
    gt_params$exact_null
  ),
  .f = function(testID, max_delta, n_cats, exact_null) {

    time_start <- Sys.time()

    repeat {

      # Generate isoform usage changes that sum to zero.
      if (n_cats == 2) {
        deltas <- c(max_delta, -max_delta)
      } else {

        repeat {
          remainder <- -max_delta
          deltas <- numeric(n_cats)
          for (i in 2:(n_cats - 1)) {
            delta_val <- (2 * rbeta(1, shape1 = 1, shape2 = 2) - 1) * abs(remainder)
            deltas[i] <- delta_val
            remainder <- remainder - delta_val
          }
          deltas[1] <- max_delta # Set the first isoform's change to the maximum.
          deltas[n_cats] <- remainder

          if (all(abs(deltas) <= max_delta) && sum(abs(deltas)) <= 2) {
            break
          }
        }
      }

      # Draw baseline isoform proportions from the Dirichlet fit for this isoform count.
      prop1 <- simulate_isoform_proportions(n_cats)

      # Apply the changes directly for two isoforms.
      if (n_cats == 2) {
        prop2 <- prop1 + deltas
      }
      # For more isoforms, try reassigning the changes to keep proportions between zero and one.
      else {
        shuffles <- 0
        repeat {
          prop2 <- prop1 + deltas
          if (all(prop2 > 0) && all(prop1 > 0) && all(prop2 <= 1)) {
            break
          } else {
            deltas <- sample(deltas)
            shuffles <- shuffles + 1
            if (shuffles >= n_cats * 2) {
              break
            }
          }
        }
      }

      if (all(prop2 > 0) && all(prop1 > 0) && all(prop2 <= 1)) {
        break
      }

      if (as.numeric(difftime(Sys.time(), time_start, units = "secs")) > timeout) {
        return(NULL)
      }
    }

    if (exact_null) {
      stopifnot(all(deltas == 0), all(prop1 == prop2))
    }

    gt_mat <- cbind(
      "testID" = testID,
      "catID" = seq_along(deltas),
      "gt.delta" = deltas,
      "gt.pct1" = prop1,
      "gt.pct2" = prop2,
      "gt.exact.null" = rep(exact_null, length(deltas))
    )
    gt_mat <- as(gt_mat, "dMatrix")
    rownames(gt_mat) <- paste0("test", testID, ".cat", seq_along(deltas))
    return(gt_mat)
  }
)

# Stop if any gene times out so the saved truth and counts describe the same set of genes.
failed <- which(vapply(gt_data, is.null, logical(1)))
if (length(failed)) {
  stop("Ground-truth generation timed out for test IDs: ", paste(failed, collapse = ", "),
       ". Increase timeout in this script and rerun; no new truth or counts were saved.")
}
log_info("Successful simulations: ", length(gt_data), "/", n_tests)

output_data <- as.data.frame(as(purrr::reduce(gt_data, rbind), "matrix"))
truth_output_file <- file.path(data_dir, paste0(run_id, "_results.tsv"))
write.table(
  output_data,
  file = truth_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_info("Done. Final ground truth table saved to ", truth_output_file)

# Observed counts ---------------------------------------------------------

# Keep the G1_C1, G2_C1, ... draw order so combining columns once per gene preserves seeded results.
simulate_transcript_counts <- function(counts1, counts2, prop1, prop2, transcript_ids) {
  columns <- vector("list", 2L * length(counts1))
  for (i in seq_along(counts1)) {
    columns[[2L * i - 1L]] <- methods::as(
      stats::rmultinom(1, counts1[i], prop1), "dgCMatrix")
    columns[[2L * i]] <- methods::as(
      stats::rmultinom(1, counts2[i], prop2), "dgCMatrix")
  }
  result <- do.call(cbind, columns)
  rownames(result) <- transcript_ids
  colnames(result) <- paste0(rep(c("G1_C", "G2_C"), length(counts1)),
                             rep(seq_along(counts1), each = 2L))
  result
}

log_info("Simulating observed counts...")

# Share cell depth factors across genes and scale each group's factors to a mean of one.
cell_size_grp1 <- rlnorm(
  n = cell_number,
  meanlog = cell_size_meanlog,
  sdlog = cell_size_sdlog
)
cell_size_grp2 <- rlnorm(
  n = cell_number,
  meanlog = cell_size_meanlog,
  sdlog = cell_size_sdlog
)
cell_size_grp1 <- cell_size_grp1 / mean(cell_size_grp1)
cell_size_grp2 <- cell_size_grp2 / mean(cell_size_grp2)

log_info(
  "Applied cell-size factors:",
  "\nGroup 1 median: ", signif(median(cell_size_grp1), 4),
  ", sdlog: ", signif(sd(log(cell_size_grp1)), 4),
  "\nGroup 2 median: ", signif(median(cell_size_grp2), 4),
  ", sdlog: ", signif(sd(log(cell_size_grp2)), 4)
)

cts_params <- tibble("gt_mat" = gt_data, "ncells" = cell_number)

# Generate counts for each gene.
obs_counts <- future_pmap(
  .options = furrr_options(stdout = FALSE, seed = TRUE),
  .progress = TRUE,
  .l = list(cts_params$gt_mat, cts_params$ncells),
  .f = function(gt_mat, ncells) {

    gt.pct1 <- gt_mat[ , "gt.pct1", drop = TRUE]
    gt.pct2 <- gt_mat[ , "gt.pct2", drop = TRUE]
    testID <- unique(gt_mat[ , "testID", drop = TRUE])
    catID <- gt_mat[ , "catID", drop = TRUE]

    # Draw each group's mean gene expression independently of its isoform count.
    mean_expr_grp1 <- rlnorm(
      n = 1,
      meanlog = gene_expr_meanlog,
      sdlog = gene_expr_sdlog
    )
    mean_expr_grp2 <- rlnorm(
      n = 1,
      meanlog = gene_expr_meanlog,
      sdlog = gene_expr_sdlog
    )

    # Draw gene counts for each cell, requiring at least one count in each group.
    repeat {
      sc_expr_grp1 <- rnbinom(
        n = ncells,
        size = sc_disp,
        mu = mean_expr_grp1 * cell_size_grp1
      )
      sc_expr_grp2 <- rnbinom(
        n = ncells,
        size = sc_disp,
        mu = mean_expr_grp2 * cell_size_grp2
      )
      if (sum(sc_expr_grp1) != 0 & sum(sc_expr_grp2) != 0) {
        break
      }
    }

    dgc_mat <- simulate_transcript_counts(
      sc_expr_grp1, sc_expr_grp2, gt.pct1, gt.pct2,
      paste0("test", testID, ".cat", catID)
    )

    return(dgc_mat)
  }
)

gene_total_counts <- vapply(
  obs_counts,
  function(x) if (is.null(x)) NA_real_ else sum(x),
  numeric(1)
)
count_bins <- cut(
  gene_total_counts,
  breaks = c(-Inf, 49, 100, 200, 499, Inf),
  labels = c("<50", "50-100", "101-200", "201-499", ">=500")
)
log_info(
  "Gene total-count distribution:\n",
  paste(capture.output(print(table(count_bins, useNA = "ifany"))), collapse = "\n")
)

# Combine the gene matrices and save the counts.
counts_output_file <- file.path(data_dir, paste0(run_id, "_sc_mat.rds"))
saveRDS(purrr::reduce(obs_counts, rbind), file = counts_output_file)

cat("\n")
log_info("Done. Saved single cell matrix to ", counts_output_file)


# Total runtime -----------------------------------------------------------
run_end <- Sys.time()
log_info("Total run time: ", as.numeric(difftime(run_end, run_start, units = "mins")), " mins")

future::plan(future::sequential)
saveRDS(utils::sessionInfo(), file.path(data_dir, paste0(run_id, "_simulation_session.rds")))
