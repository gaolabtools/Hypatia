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
    permutation = FALSE,
    bootstrap = FALSE,
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
    cell.dispersion = FALSE,
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

  expect_equal(
    diu$data[names(baseline$run_diu$data)],
    baseline$run_diu$data,
    tolerance = 1e-12
  )
  unchanged_diu_stats <- setdiff(
    names(baseline$run_diu$stats),
    "effect.size"
  )
  expect_equal(
    diu$stats[unchanged_diu_stats],
    baseline$run_diu$stats[unchanged_diu_stats],
    tolerance = 1e-12
  )
  expect_equal(
    diu$stats$cramers.v,
    baseline$run_diu$stats$effect.size,
    tolerance = 1e-12
  )
  expect_true(all(lengths(diu$data$boot.prop.1) == 0))
  expect_true(all(lengths(diu$data$boot.prop.2) == 0))
  expect_true(all(lengths(diu$data$boot.prop.diff) == 0))
  expect_false("bootstrap" %in% names(diu$stats))
  expect_true(all(is.na(diu$stats$boot.prop.diff.mean)))
  expect_true(all(is.na(diu$stats$boot.prop.diff.lower)))
  expect_true(all(is.na(diu$stats$boot.prop.diff.upper)))
  expect_true(all(diu$stats$boot.valid.iter == 0L))

  baseline_div_data <- baseline$run_div$data[
    order(
      baseline$run_div$data$group.1,
      baseline$run_div$data$group.2,
      baseline$run_div$data$gene
    ),
    ,
    drop = FALSE
  ]
  rownames(baseline_div_data) <- NULL
  expect_equal(
    div$data[names(baseline$run_div$data)],
    baseline_div_data,
    tolerance = 1e-12
  )
  unchanged_div_stats <- setdiff(
    names(baseline$run_div$stats),
    c("log2FC", "div.class.1", "div.class.2")
  )
  baseline_div_stats <- baseline$run_div$stats[
    order(
      baseline$run_div$stats$group.1,
      baseline$run_div$stats$group.2,
      baseline$run_div$stats$padj,
      baseline$run_div$stats$gene
    ),
    ,
    drop = FALSE
  ]
  rownames(baseline_div_stats) <- NULL
  expect_equal(
    div$stats[unchanged_div_stats],
    baseline_div_stats[unchanged_div_stats],
    tolerance = 1e-12
  )
  expected_effective_1 <- c(
    ANKRD10 = 1L, GNG7 = 2L, NAPEPLD = 1L, NPEPPS = 2L, PNISR = 2L
  )
  expected_effective_2 <- c(
    ANKRD10 = 1L, GNG7 = 2L, NAPEPLD = 2L, NPEPPS = 2L, PNISR = 2L
  )
  expect_equal(
    div$stats$n.effective.1,
    unname(expected_effective_1[div$stats$gene])
  )
  expect_equal(
    div$stats$n.effective.2,
    unname(expected_effective_2[div$stats$gene])
  )
  expect_equal(
    div$stats$div.class.1,
    rep("polyform", 5)
  )
  expect_equal(
    div$stats$div.class.2,
    rep("polyform", 5)
  )
  expect_true(all(is.na(div$data$cell.n.1)))
  expect_true(all(is.na(div$data$cell.div.mean.2)))

  baseline_dei_data <- baseline$run_dei[
    order(
      baseline$run_dei$group.1,
      baseline$run_dei$group.2,
      baseline$run_dei$gene,
      baseline$run_dei$transcript
    ),
    ,
    drop = FALSE
  ]
  rownames(baseline_dei_data) <- NULL
  expect_equal(
    dei$data[c("group.1", "group.2", "gene", "transcript", "pct.1", "pct.2", "avgExpr.1", "avgExpr.2")],
    baseline_dei_data[c("group.1", "group.2", "gene", "transcript", "pct.1", "pct.2", "avgExpr.1", "avgExpr.2")],
    tolerance = 1e-12
  )
  baseline_dei_stats <- baseline$run_dei[
    order(
      baseline$run_dei$group.1,
      baseline$run_dei$group.2,
      baseline$run_dei$padj,
      baseline$run_dei$gene,
      baseline$run_dei$transcript
    ),
    ,
    drop = FALSE
  ]
  rownames(baseline_dei_stats) <- NULL
  expect_equal(
    dei$stats[c("group.1", "group.2", "gene", "transcript", "log2FC", "pval", "padj")],
    baseline_dei_stats[c("group.1", "group.2", "gene", "transcript", "log2FC", "pval", "padj")],
    tolerance = 1e-12
  )
})
