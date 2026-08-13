test_that("RunDEI outputs a data frame", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    gbm <- NormalizeCounts(gbm)
    res <- RunDEI(gbm)
    RunDEI(gbm, group.1 = "Tumor")
    RunDEI(gbm, group.1 = "Tumor", group.2 = "Oligodendrocyte")
    RunDEI(gbm, group.1 = "Tumor", group.2 = c("Oligodendrocyte", "Astrocyte"))
    })


  expect_class(res, "data.frame")
  expect_true(nrow(res) >= 3)

})

test_that("RunDEI respects transcript filters", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  transcripts <- rownames(gbm)[1:2]
  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = transcripts,
    min.pct = 0,
    quiet = TRUE
  )

  expect_setequal(res$transcript, transcripts)
})

test_that("RunDEI requires cells on both sides of a comparison", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  expect_error(
    RunDEI(gbm, group.1 = unique(colData(gbm)$cell_type), quiet = TRUE),
    "at least one cell in both groups"
  )
})

test_that("RunDEI works with active transcript IDs", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)
  gbm <- SetTranscripts(gbm, id = "transcript_name")

  transcripts <- rowData(gbm)$transcript_name[1:2]
  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = transcripts,
    min.pct = 0,
    quiet = TRUE
  )

  expect_setequal(res$transcript, transcripts)
})

test_that("RunDEI reports when no transcripts pass detection thresholds", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  expect_error(
    RunDEI(
      gbm,
      group.1 = "Tumor",
      group.2 = "Oligodendrocyte",
      min.pct = 1,
      quiet = TRUE
    ),
    "0 transcripts passed detection thresholds"
  )
})

test_that("RunDEI supports stats::p.adjust methods", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )
  gbm <- NormalizeCounts(gbm, quiet = TRUE)

  res <- RunDEI(
    gbm,
    group.1 = "Tumor",
    group.2 = "Oligodendrocyte",
    transcripts = rownames(gbm)[1:2],
    min.pct = 0,
    p.adj = "none",
    quiet = TRUE
  )

  expect_equal(res$padj, res$pval)
  expect_error(
    RunDEI(gbm, p.adj = "Bonferroni", quiet = TRUE),
    "element of set"
  )
})
