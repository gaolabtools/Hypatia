test_that("analysis outputs match the regression baseline", {

  baseline <- readRDS(test_path("fixtures", "run-output-baseline.rds"))

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )
  gbm <- SetGenes(gbm, "gene_name", quiet = TRUE)
  gbm <- SetTranscripts(gbm, "transcript_name")

  diu <- RunDIU(
    gbm,
    group.1 = "Glut neuron",
    group.2 = "GABA neuron",
    genes = c("MEG3", "MALAT1"),
    min.gene.pct = 0,
    min.gene.cts = 0,
    min.tx.cts = 1,
    quiet = TRUE
  )

  set.seed(1024)
  div <- RunDIV(
    gbm,
    group.1 = "Glut neuron",
    group.2 = "GABA neuron",
    genes = c("NAPEPLD", "NPEPPS", "GNG7", "ANKRD10", "PNISR"),
    min.gene.pct = 0,
    min.gene.cts = 0,
    min.tx.cts = 1,
    boot.iter = 10,
    quiet = TRUE
  )

  dei_transcripts <- rowData(gbm)$transcript_name[rowData(gbm)$gene_name == "MEG3"][1:5]
  gbm_dei <- SubsetTranscripts(gbm, transcripts = dei_transcripts, quiet = TRUE)
  gbm_dei <- NormalizeCounts(gbm_dei, quiet = TRUE)
  dei <- RunDEI(
    gbm_dei,
    group.1 = "Glut neuron",
    group.2 = "GABA neuron",
    min.pct = 0,
    quiet = TRUE
  )

  expect_equal(diu, baseline$run_diu, tolerance = 1e-12)
  expect_equal(div, baseline$run_div, tolerance = 1e-12)
  expect_equal(dei, baseline$run_dei, tolerance = 1e-12)
})
