# DIU, DEI, DIV analysis on GBM02PR, RCC12P, and HR182

suppressPackageStartupMessages(library(Hypatia))
data_dir <- "manuscript/figure2/data"
input <- file.path(data_dir, "hypatia_objects.rda")
if (!file.exists(input)) stop("Missing Hypatia objects: ", input)
objects <- new.env(parent = emptyenv())
loaded <- load(input, envir = objects)
datasets <- c("gbm", "rcc", "heart")
if (!all(datasets %in% loaded)) stop("Expected gbm, rcc, and heart objects in ", input)

# Keep the seed and analysis order fixed to reproduce the resampling results.
set.seed(1029)
run_comparisons <- function(fun, suffix, filename, ...) {
  results <- setNames(lapply(datasets, function(dataset) fun(objects[[dataset]], ...)),
                      paste0(datasets, "_", suffix))
  save(list = names(results), envir = list2env(results),
       file = file.path(data_dir, filename))
}
run_comparisons(Hypatia::RunDIU, "diu", "diu.rda", permutation = TRUE)
run_comparisons(Hypatia::RunDEI, "dei", "dei.rda", only.pos = TRUE)
run_comparisons(Hypatia::RunDIV, "div", "div_new.rda")

contrasts <- list(
  neu_div = c("gbm", "GABA neuron", "Glut neuron"),
  cmfb_div = c("heart", "Cardiomyocyte", "Fibroblast"),
  bcell_div = c("rcc", "plasma cell", "follicular B cell")
)
results <- lapply(contrasts, function(contrast) {
  Hypatia::RunDIV(objects[[contrast[1]]], group.1 = contrast[2],
                  group.2 = contrast[3], include.single = FALSE)
})
save(list = names(results), envir = list2env(results),
     file = file.path(data_dir, "div_subset_new.rda"))
saveRDS(utils::sessionInfo(), file.path(data_dir, "comparisons_session.rds"))
