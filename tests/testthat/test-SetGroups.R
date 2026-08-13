test_that("SetGroups updates the active.group.id object metadata", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error(gbm <- SetGroups(gbm, id = "cell_type"))

  expect_equal(metadata(gbm)$active.group.id, "cell_type")

})

test_that("SetGroups validates grouping columns", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  colData(gbm)$cell_factor <- factor(colData(gbm)$cell_type)
  colData(gbm)$`cell type` <- colData(gbm)$cell_type
  colData(gbm)$one_group <- "A"
  colData(gbm)$missing_group <- colData(gbm)$cell_type
  colData(gbm)$missing_group[1] <- NA_character_

  expect_no_error(gbm_factor <- SetGroups(gbm, id = "cell_factor"))
  expect_equal(metadata(gbm_factor)$active.group.id, "cell_factor")

  expect_no_error(gbm_space <- SetGroups(gbm, id = "cell type"))
  expect_equal(metadata(gbm_space)$active.group.id, "cell type")

  expect_error(SetGroups(gbm, id = "nCount"), "character or factor")
  expect_error(SetGroups(gbm, id = "one_group"), "at least 2")
  expect_error(SetGroups(gbm, id = "missing_group"))
  expect_error(SetGroups(gbm, id = "not_a_column"))
})

test_that("SetGroups factor columns work in downstream grouping", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )
  colData(gbm)$cell_factor <- factor(colData(gbm)$cell_type)
  gbm <- SetGenes(gbm, "gene_name", quiet = TRUE)
  gbm <- SetTranscripts(gbm, "transcript_name")
  gbm <- SetGroups(gbm, "cell_factor")

  expect_no_error(
    RunDIU(
      gbm,
      group.1 = "Glut neuron",
      group.2 = "GABA neuron",
      genes = "MEG3",
      quiet = TRUE
    )
  )
})
