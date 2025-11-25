test_that("SetGenes updates the active.gene.id object metadata", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- SetGenes(gbm, id = "gene_name"))

  expect_true(metadata(gbm)$active.gene.id == "gene_name")

})
