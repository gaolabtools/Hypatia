test_that("group.by supports column names with spaces", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )
  colData(gbm)[["cell type"]] <- colData(gbm)$cell_type
  gbm <- SetGenes(gbm, "gene_name", quiet = TRUE)
  gbm <- SetTranscripts(gbm, "transcript_name")
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  expect_no_error(
    RunDIU(
      gbm,
      group.by = "cell type",
      group.1 = "Glut neuron",
      group.2 = "GABA neuron",
      genes = "MEG3",
      quiet = TRUE
    )
  )

  expect_no_error(
    RunDIV(
      gbm,
      group.by = "cell type",
      group.1 = "Glut neuron",
      group.2 = "GABA neuron",
      genes = c("NAPEPLD", "NPEPPS", "GNG7", "ANKRD10", "PNISR"),
      boot.iter = 2,
      min.gene.pct = 0,
      quiet = TRUE
    )
  )

  expect_no_error(
    RunDEI(
      gbm,
      group.by = "cell type",
      group.1 = "Glut neuron",
      group.2 = "GABA neuron",
      transcripts = rowData(gbm)$transcript_name[1:2],
      min.pct = 0,
      quiet = TRUE
    )
  )

  expect_no_error(
    GetUsage(gbm, genes = "MEG3", group.by = "cell type", group.subset = "Glut neuron", quiet = TRUE)
  )

  expect_no_error(
    GetExpression(
      gbm,
      transcripts = rowData(gbm)$transcript_name[1:2],
      group.by = "cell type",
      group.subset = "Glut neuron",
      quiet = TRUE
    )
  )

  expect_no_error(
    PlotCellQC(gbm, group.by = "cell type", combine = FALSE)
  )
})
