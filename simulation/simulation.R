run_start <- Sys.time()

# Settings ------------------------------------------------------------

## Run
run_id <- "run_x"
n_cores <- 20
working_dir <- "simulation/"
seed <- 1029
set.seed(seed)

## Simulation
n_tests <- 30000 # number of 'genes' to simulate
cell_number <- 500 # cells per group (2 groups will be simulated)
avg_gene_expr <- 2 # gamma scale for mean gene expr
gene_expr_shape <- 0.2 # gamma shape for mean gene expr
sc_disp <- 0.1 # rnbinom size for dispersion param
timeout <- 60 # timeout for ground truth simulations (seconds)

## Testing
# run_id <- "test_run"
# n_cores <- 5
# n_tests <- 500
# cell_number <- 100

# Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(furrr)
  library(logger)
})

setwd(working_dir)
logfile <- paste0(run_id, ".log")
invisible(file.create(logfile))
log_appender(appender_tee(logfile))
log_info("--- Run details ---",
         "\nRun name: ", run_id,
         "\nWorking directory: ", working_dir,
         "\nCores: ", n_cores,
         "\nSeed: ", seed,
         "\nNumber of tests: ", n_tests,
         "\nCells per group: ", cell_number,
         "\nAverage gene expression: ", avg_gene_expr,
         "\nGene expression shape: ", gene_expr_shape,
         "\nDispersion parameter: ", sc_disp,
         "\nTimeout (s): ", timeout)

plan(multicore, workers = n_cores)


# Ground truth ------------------------------------------------------------
log_info("Simulating ground truth deltas and isoform proportion sets for ", n_tests, " genes...")

## Maximum delta for each test
max_deltas <- numeric(0)
while (length(max_deltas) < n_tests) {
  # gamma
  # max_delta_val <- rgamma(n_tests, shape = 1, scale = 0.1)
  # max_deltas <- c(max_deltas, max_delta_val[max_delta_val <= 1 & max_delta_val >= 0])

  # beta
  max_deltas <- rbeta(n = n_tests, shape1 = 0.5, shape2 = 10)
}
max_deltas <- max_deltas[1:n_tests]

## Distribution of categories
cats_range <- c(2, 10)
cats <- integer(0)
while (length(cats) < n_tests) {
  cat_vals <- rpois(n_tests, lambda = cats_range[1])
  cats <- c(cats, cat_vals[cat_vals >= cats_range[1] & cat_vals <= cats_range[2]])
}
cats <- cats[1:n_tests]

## Delta vectors and Proportions
gt_params <- tibble(
  testID = 1:length(max_deltas), # each test gets an ID
  max_delta = max_deltas,
  n_cats = cats
  )

# gt_params <- gt_params %>% sample_n(5) # for testing purposes

gt_data <- future_pmap(
  .options = furrr_options(stdout = FALSE, seed = TRUE),
  .progress = TRUE,
  .l = list(gt_params$testID, gt_params$max_delta, gt_params$n_cats),
  .f = function(testID, max_delta, n_cats) {

    time_start <- Sys.time()

    repeat {

      time_elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))

      ## Delta vector
      if (n_cats == 2) {
        deltas <- c(max_delta, -max_delta)
      } else {

        repeat {
          # uniform
          # deltas <- runif(n_cats - 2, min = -max_delta * (2 / n_cats), max = max_delta * (2 / n_cats)) # '2 /' for heterogeneity
          # last_delta <- -max_delta - sum(deltas)
          # deltas <- c(max_delta, deltas, last_delta)

          # beta
          remainder <- -max_delta
          deltas <- numeric(n_cats)
          for (i in 2:(n_cats - 1)) {
            delta_val <- (2 * rbeta(1, shape1 = 1, shape2 = 2) - 1) * abs(remainder)
            deltas[i] <- delta_val
            remainder <- remainder - delta_val
          }
          deltas[1] <- max_delta # first is max delta
          deltas[n_cats] <- remainder

          if (all(abs(deltas) <= max_delta) && sum(abs(deltas)) <= 2) {
            break
          }
        }
      }

      ## Proportions
      prop1 <- numeric(n_cats)
      prop1[1] <- rbeta(1, shape1 = 3, shape2 = 1) # control max prop
      remainder <- 1 - prop1[1]

      # case for 2 isoforms
      if (n_cats == 2) {
        prop1[2] <- remainder
        prop2 <- prop1 + deltas
      }
      # all other cases
      else {
        for (i in 2:(n_cats - 1)) {
          prop_val <- rbeta(1, shape1 = 1, shape2 = 2) * abs(remainder) # control remaining prop
          prop1[i] <- prop_val
          remainder <- remainder - prop_val
        }

        prop1[n_cats] <- remainder

        # try shuffling delta
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

      if (time_elapsed > timeout) {
        return(NULL)
        break
      }
    }

    gt_mat <- cbind("testID" = testID, "catID" = 1:length(deltas), "gt.delta" = deltas, "gt.pct1" = prop1, "gt.pct2" = prop2)
    gt_mat <- as(gt_mat, "dMatrix")
    rownames(gt_mat) <- paste0("test", testID, ".cat", 1:length(deltas))
    return(gt_mat)
  }
)

cat("\n")
log_info("Done. Successful simulations: ", length(gt_data), "/", n_tests, "(", round(length(gt_data) / n_tests * 100, 2), "%)")
# saveRDS(gt_data, file = paste0(run_id, "_gt_mat.rds"))
# quit(save = "no")

# Observed counts ---------------------------------------------------------
log_info("Simulating observed counts...")

cts_params <- tibble("gt_mat" = gt_data, "ncells" = cell_number)

## Single cell
# note: repeats ensure that each category will have at least 1 count across all cells
obs_counts <- future_pmap(
  # pmap(
  .options = furrr_options(stdout = FALSE, seed = TRUE),
  .progress = TRUE,
  .l = list(cts_params$gt_mat, cts_params$ncells),
  .f = function(gt_mat, ncells) {

    gt.pct1 <- gt_mat[ , "gt.pct1", drop = TRUE]
    gt.pct2 <- gt_mat[ , "gt.pct2", drop = TRUE]
    testID <- unique(gt_mat[ , "testID", drop = TRUE])
    catID <- gt_mat[ , "catID", drop = TRUE]

    ncats <- length(gt.pct1)

    ## mean cell expression of the test (gene)
    repeat {
      mean_expr_grp1 <- rgamma(n = 1, shape = gene_expr_shape, scale = avg_gene_expr)
      mean_expr_grp2 <- rgamma(n = 1, shape = gene_expr_shape, scale = avg_gene_expr)
      # make sure values are not too small
      if (mean_expr_grp1 >= 0.02 & mean_expr_grp2 >= 0.02) {
        break
      }
    }

    ## single cell expression of the test (gene)
    repeat {
      sc_expr_grp1 <- rnbinom(n = ncells, size = sc_disp, mu = mean_expr_grp1)
      sc_expr_grp2 <- rnbinom(n = ncells, size = sc_disp, mu = mean_expr_grp2)
      # must have at least 1 count in each cell group
      if (sum(sc_expr_grp1) != 0 & sum(sc_expr_grp2) != 0) {
        break
      }
    }

    ## single cell expression of the categories (isoforms)

    repeat {
      dgc_mat <- sparseMatrix(i = integer(0),
                              j = integer(0),
                              x = numeric(0),
                              dims = c(ncats, 0))

      for (i in 1:ncells) {
        # if (any(gt.pct1 == 0) || any(gt.pct2 == 0)) {
        #   stop("Categories have 0 probability")
        # }
        counts_grp1 <- rmultinom(n = 1, size = sc_expr_grp1[i], prob = gt.pct1) # same order as catID
        counts_grp1 <- as(counts_grp1, "dgCMatrix")
        colnames(counts_grp1) <- paste0("G1_C", i) # cell group ID is coded into column name
        rownames(counts_grp1) <- paste0("test", testID, ".cat", catID) # test and cat ID is coded into row name

        counts_grp2 <- rmultinom(n = 1, size = sc_expr_grp2[i], prob = gt.pct2)
        counts_grp2 <- as(counts_grp2, "dgCMatrix")
        colnames(counts_grp2) <- paste0("G2_C", i)
        rownames(counts_grp2) <- paste0("test", testID, ".cat", catID)

        dgc_mat <- cbind(dgc_mat, counts_grp1, counts_grp2)
      }

      # ensure that each category has at least 1 count across all cells
      # now set at >= 0 ... if set at >0, stuck when pct are near 0 with low expr
      if (all(rowSums(dgc_mat) >= 0)) {
        break
      }
    }

    return(dgc_mat)

  }
)

## save combined mat
saveRDS(purrr::reduce(obs_counts, rbind), file = paste0(run_id, "_sc_mat.rds"))

cat("\n")
log_info("Done. Saving single cell matrix to ", paste0(getwd(), "/", run_id, "_sc_mat.rds"), "...")

# Observed results --------------------------------------------------------
log_info("Done. Summarizing observed results...")

res_params <- tibble("sc_mat" = obs_counts)

obs_data <- future_pmap(
  .options = furrr_options(stdout = FALSE, seed = TRUE),
  .progress = TRUE,
  .l = list(res_params$sc_mat),
  .f = function(sc_mat) {

    ## pseudobulk counts
    cell_groups <- factor(sub("_.*", "", colnames(sc_mat)))
    bulk_mat <- sc_mat %*% model.matrix(~ 0 + cell_groups)
    colnames(bulk_mat) <- str_replace_all(levels(cell_groups), "G", "cts")

    ## proportions and delta
    grp_sums <- colSums(bulk_mat)
    grp_props <- bulk_mat %*% Diagonal(x = 1 / grp_sums)
    colnames(grp_props) <- c("prop1", "prop2")
    grp_props_delta <- cbind(grp_props, "delta" = grp_props[, "prop2"] - grp_props[, "prop1"])
    res_mat <- cbind(bulk_mat, grp_props_delta)

    ## Chi-sq test
    cont_table <- as.matrix(res_mat[, c("cts1", "cts2")])
    cont_table <- cont_table[!(cont_table[,1] == 0 & cont_table[,2] == 0), ] # removes rows that are both 0
    if (is.null(nrow(cont_table))) {
      res_mat <- cbind(res_mat, "X2.pval" = NA, "X2.stat" = NA, "X2.v" = NA, "FE.pval" = NA)
    } else {
      chisq_res <- suppressWarnings(chisq.test(cont_table))
      fisher_res <- suppressWarnings(fisher.test(cont_table, simulate.p.value = TRUE))
      chisq.v <- sqrt(chisq_res$statistic / (sum(cont_table) * 1)) # always 2 columns
      res_mat <- cbind(res_mat,
                       "X2.pval" = chisq_res$p.value,
                       "X2.stat" = chisq_res$statistic,
                       "X2.v" = chisq.v,
                       "FE.pval" = fisher_res$p.value)
    }

    return(res_mat)
  }
)

cat("\n")

# Output table ------------------------------------------------------------
log_info("Done. Generating table of results...")

output_data <- mapply(cbind, gt_data, obs_data, SIMPLIFY = FALSE)
output_data <- as(purrr::reduce(output_data, rbind), "matrix")

write.table(output_data, file = paste0(run_id, "_results.tsv"), sep = "\t", quote = FALSE)

# write.table(output_data, file = "test.tsv", sep = "\t", quote = FALSE)
# read.table("test.tsv", sep = "\t", header = TRUE, row.names = 1)

log_info("Done. Final table saved to ", paste0(getwd(), "/", run_id, "_results.tsv"))


# -------------------------------------------------------------------------
run_end <- Sys.time()
log_info("Total run time: ", as.numeric(difftime(run_end, run_start, units = "mins")), " mins")
