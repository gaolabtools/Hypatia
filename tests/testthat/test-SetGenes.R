test_that("SetGenes updates the active.gene.id object metadata", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- SetGenes(gbm, id = "gene_name"))

  expect_equal(metadata(gbm)$active.gene.id, "gene_name")

})

test_that("SetGenes validates gene ID columns", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  rowData(gbm)$gene_missing <- rowData(gbm)$gene_id
  rowData(gbm)$gene_missing[1] <- NA_character_
  rowData(gbm)$gene_factor <- factor(rowData(gbm)$gene_id)

  expect_error(SetGenes(gbm, id = "not_a_column", quiet = TRUE))
  expect_error(SetGenes(gbm, id = "gene_missing", quiet = TRUE))
  expect_error(SetGenes(gbm, id = "gene_factor", quiet = TRUE))
  expect_message(SetGenes(gbm, id = "gene_name", quiet = FALSE), "different number")
  expect_message(SetGenes(gbm, id = "gene_name", quiet = TRUE), NA)
})
