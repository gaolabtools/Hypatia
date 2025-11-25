test_that("CreateSCE returns a SingleCellExperiment object", {

  expect_class(
    CreateSCE(
    countData = gbm_countData,
    colData = gbm_colData,
    rowData = gbm_rowData,
    ),
  "SingleCellExperiment")

  expect_class(
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.gene.id = "gene_name",
      active.group.id = "cell_type",
      active.transcript.id = "transcript_name",
      gtf = gbm_gtf,
      gtf.transcript.id = "transcript_id",
    ),
    "SingleCellExperiment")

})
