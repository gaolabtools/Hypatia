test_that("RunDIV works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE,
      active.group.id = "cell_type"
    )

  expect_no_error({
    res <- RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10)
    RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "Shannon", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 0, min.gene.pct = 0, boot.iter = 10, entropy.use = "NormalizedShannon", prop.thresh = 0.2, quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "Renyi", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "NormalizedRenyi", prop.thresh = 0.2, quiet = TRUE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "GiniSimpson", quiet = FALSE)
    # RunDIV(gbm, min.gene.cts = 10, min.gene.pct = 0.05, boot.iter = 10, entropy.use = "InverseSimpson", quiet = FALSE)
  })

  expect_class(res, "list")
  expect_class(res$data, "data.frame")
  expect_class(res$stats, "data.frame")
  expect_true("div.diff" %in% names(res$stats))
  expect_false("delta.div" %in% names(res$stats))

})

test_that("RunDIV requires cells on both sides of a comparison", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  expect_error(
    RunDIV(gbm, group.1 = unique(colData(gbm)$cell_type), quiet = TRUE),
    "at least one cell in both groups"
  )
})

test_that("RunDIV works with active transcript IDs", {

  countData <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 2, 3, 4),
    x = c(5, 3, 2, 6),
    dims = c(2, 4),
    dimnames = list(paste0("tx", 1:2), paste0("cell", 1:4))
  )
  colData <- data.frame(group = c("A", "A", "B", "B"), row.names = colnames(countData))
  rowData <- data.frame(
    gene_id = rep("gene1", 2),
    tx_name = c("tx_a", "tx_b"),
    row.names = rownames(countData)
  )
  object <- CreateSCE(countData, colData, rowData, active.group.id = "group", quiet = TRUE)
  object <- SetTranscripts(object, id = "tx_name")

  expect_no_error({
    res <- RunDIV(
      object,
      group.1 = "A",
      group.2 = "B",
      min.gene.cts = 0,
      min.gene.pct = 0,
      min.tx.cts = 0,
      boot.iter = 3,
      boot.fraction = 1,
      genes = "gene1",
      p.adj = "none",
      quiet = TRUE
    )
  })
  expect_equal(res$stats$gene, "gene1")
  expect_equal(res$stats$padj, res$stats$pval)
})

test_that("RunDIV rejects p.adjust method names not used by stats::p.adjust", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      active.group.id = "cell_type",
      quiet = TRUE
    )

  expect_error(
    RunDIV(gbm, p.adj = "Bonferroni", quiet = TRUE),
    "element of set"
  )
})

test_that("RunDIV calculates diversity using only transcripts that pass min.tx.cts", {

  countData <- Matrix::Matrix(
    matrix(
      c(
        10, 10, 8, 8,
         5,  5, 10, 10,
         1,  1, 1, 1
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  colData <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(countData)
  )
  rowData <- data.frame(
    gene_id = rep("gene1", 3),
    row.names = rownames(countData)
  )
  object <- CreateSCE(
    countData,
    colData,
    rowData,
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(1024)
  res <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 3,
    boot.iter = 4,
    boot.fraction = 1,
    include.single = FALSE,
    p.adj = "none",
    quiet = TRUE
  )

  expected_div_1 <- (1 - ((2 / 3)^3 + (1 / 3)^3)) / 2
  expected_div_2 <- (1 - ((4 / 9)^3 + (5 / 9)^3)) / 2

  expect_equal(res$data$n.transcripts, 2)
  expect_equal(res$data$div.1[[1]], rep(expected_div_1, 4), tolerance = 1e-12)
  expect_equal(res$data$div.2[[1]], rep(expected_div_2, 4), tolerance = 1e-12)
  expect_equal(res$stats$avgDiv.1, expected_div_1, tolerance = 1e-12)
  expect_equal(res$stats$avgDiv.2, expected_div_2, tolerance = 1e-12)
})

test_that("RunDIV supports a fixed number of bootstrap cell draws", {

  countData <- Matrix::Matrix(
    matrix(
      c(
        9, 3, 8, 2,
        1, 7, 2, 8
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:2), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  colData <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(countData)
  )
  rowData <- data.frame(
    gene_id = rep("gene1", 2),
    row.names = rownames(countData)
  )
  object <- CreateSCE(
    countData,
    colData,
    rowData,
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(1024)
  grp1_idx <- sample(seq_len(2), size = 1, replace = TRUE)
  grp2_idx <- sample(seq_len(2), size = 1, replace = TRUE)
  tsallis <- function(x) {
    prop <- x / sum(x)
    (1 - sum(prop^3)) / 2
  }
  expected_div_1 <- tsallis(countData[, grp1_idx])
  expected_div_2 <- tsallis(countData[, grp2_idx + 2])

  set.seed(1024)
  res <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 0,
    boot.iter = 1,
    boot.ncells = 1,
    include.single = FALSE,
    p.adj = "none",
    quiet = TRUE
  )

  expect_equal(res$data$div.1[[1]], expected_div_1, tolerance = 1e-12)
  expect_equal(res$data$div.2[[1]], expected_div_2, tolerance = 1e-12)
  expect_equal(res$stats$avgDiv.1, expected_div_1, tolerance = 1e-12)
  expect_equal(res$stats$avgDiv.2, expected_div_2, tolerance = 1e-12)

  expect_error(
    RunDIV(object, boot.fraction = 0.5, boot.ncells = 1, quiet = TRUE),
    "only one"
  )
  expect_error(
    RunDIV(object, boot.ncells = 1.5, quiet = TRUE),
    "type 'count'"
  )
})

test_that("RunDIV reports biological cell-to-cell diversity dispersion", {
  countData <- Matrix::Matrix(
    matrix(
      c(
        8, 0, 4, 0, 3, 1,
        2, 0, 0, 0, 1, 3,
        0, 0, 1, 0, 0, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tx", 1:3), paste0("cell", 1:6))
    ),
    sparse = TRUE
  )
  object <- CreateSCE(
    countData,
    data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(countData)),
    data.frame(gene_id = rep("gene1", 3), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(1024)
  without_dispersion <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 2,
    boot.iter = 2,
    boot.fraction = 1,
    include.single = FALSE,
    cell.dispersion = FALSE,
    p.adj = "none",
    quiet = TRUE
  )
  set.seed(1024)
  res_messages <- capture.output(
    res <- RunDIV(
      object,
      group.1 = "A",
      group.2 = "B",
      min.gene.cts = 0,
      min.gene.pct = 0,
      min.tx.cts = 2,
      cell.dispersion = TRUE,
      boot.iter = 2,
      boot.fraction = 1,
      include.single = FALSE,
      p.adj = "none",
      quiet = FALSE
    ),
    type = "message"
  )
  expect_true(any(grepl("Calculating cell-to-cell dispersion...", res_messages, fixed = TRUE)))

  div_a <- c(0.24, 0)
  div_b <- c(0.28125, 0.28125)
  expect_equal(res$data$n.transcripts, 2)
  expect_equal(res$data$cell.n.1, 2L)
  expect_equal(res$data$cell.n.2, 2L)
  expect_equal(res$data$cell.div.mean.1, mean(div_a), tolerance = 1e-12)
  expect_equal(res$data$cell.div.median.1, median(div_a), tolerance = 1e-12)
  expect_equal(res$data$cell.div.sd.1, sd(div_a), tolerance = 1e-12)
  expect_equal(res$data$cell.div.iqr.1, IQR(div_a), tolerance = 1e-12)
  expect_equal(res$data$cell.div.mean.2, mean(div_b), tolerance = 1e-12)
  expect_false("cell.monoform.frac.1" %in% names(res$data))
  expect_false("cell.monoform.frac.2" %in% names(res$data))
  legacy_cols <- setdiff(names(without_dispersion$data), grep("^cell\\.", names(without_dispersion$data), value = TRUE))
  expect_equal(res$data[legacy_cols], without_dispersion$data[legacy_cols], tolerance = 1e-12)
  expect_equal(res$stats, without_dispersion$stats, tolerance = 1e-12)

  normalized <- Hypatia:::.CellDiversityDispersion(
    Matrix::Matrix(
      matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE,
             dimnames = list(c("tx1", "tx2"), c("c1", "c2"))),
      sparse = TRUE
    ),
    c("gene1", "gene1"),
    Hypatia:::.DiversityFunction("NormalizedShannon")
  )
  expect_equal(normalized$cell.div.mean, 0)
  expect_false("cell.monoform.frac" %in% names(normalized))

  empty <- Hypatia:::.CellDiversityDispersion(
    Matrix::Matrix(matrix(0, nrow = 2, ncol = 1,
                          dimnames = list(c("tx1", "tx2"), "c1")), sparse = TRUE),
    c("gene1", "gene1"),
    Hypatia:::.DiversityFunction("Tsallis")
  )
  expect_equal(empty$cell.n, 0L)
  expect_true(is.na(empty$cell.div.mean))

  disabled <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 2,
    boot.iter = 1,
    boot.ncells = 1,
    cell.dispersion = FALSE,
    quiet = TRUE
  )
  expect_type(disabled$data$cell.n.1, "integer")
  expect_true(all(is.na(disabled$data$cell.n.1)))
  expect_true(all(is.na(disabled$data$cell.div.mean.2)))
})

test_that("diversity functions classify genes by effective isoform counts", {
  countData <- Matrix::Matrix(
    matrix(
      c(
        5, 5, 5, 5,
        1, 1, 2, 2
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("tx1", "tx2"), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  object <- CreateSCE(
    countData,
    data.frame(group = c("A", "A", "B", "B"), row.names = colnames(countData)),
    data.frame(gene_id = c("gene1", "gene1"), row.names = rownames(countData)),
    active.group.id = "group",
    quiet = TRUE
  )

  set.seed(1024)
  run_default <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 3,
    prop.thresh = 0.2,
    boot.iter = 3,
    boot.fraction = 1,
    p.adj = "none",
    quiet = TRUE
  )
  set.seed(1024)
  run_high_threshold <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    entropy.thresh = 0.5,
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 3,
    prop.thresh = 0.5,
    boot.iter = 3,
    boot.fraction = 1,
    p.adj = "none",
    quiet = TRUE
  )

  expect_equal(run_default$stats$n.effective.1, 1L)
  expect_equal(run_default$stats$n.effective.2, 2L)
  expect_equal(run_default$stats$div.class.1, "monoform")
  expect_equal(run_default$stats$div.class.2, "polyform")
  expect_equal(run_high_threshold$stats$n.effective.2, 1L)
  expect_equal(run_high_threshold$stats$div.class.2, "monoform")
  expect_equal(run_default$data, run_high_threshold$data, tolerance = 1e-12)
  numeric_stats <- c("avgDiv.1", "avgDiv.2", "div.diff", "pval", "padj")
  expect_equal(
    run_default$stats[numeric_stats],
    run_high_threshold$stats[numeric_stats],
    tolerance = 1e-12
  )

  set.seed(1024)
  run_shannon_top2 <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    entropy.use = "Shannon",
    top.n = 2,
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 3,
    boot.iter = 3,
    boot.fraction = 1,
    p.adj = "none",
    quiet = TRUE
  )
  expect_equal(run_shannon_top2$stats$div.class.1, "monoform")
  expect_equal(run_shannon_top2$stats$div.class.2, "polyform")

  set.seed(1024)
  run_tsallis_one <- RunDIV(
    object,
    group.1 = "A",
    group.2 = "B",
    entropy.use = "Tsallis",
    order = 1,
    min.gene.cts = 0,
    min.gene.pct = 0,
    min.tx.cts = 3,
    boot.iter = 3,
    boot.fraction = 1,
    p.adj = "none",
    quiet = TRUE
  )
  expect_equal(
    run_tsallis_one$stats[c("avgDiv.1", "avgDiv.2")],
    run_shannon_top2$stats[c("avgDiv.1", "avgDiv.2")],
    tolerance = 1e-12
  )
  expect_true(is.na(run_tsallis_one$stats$div.class.1))
  expect_true(is.na(run_tsallis_one$stats$div.class.2))

  get_result <- GetDiversity(
    object,
    genes = "gene1",
    min.tx.cts = 3,
    prop.thresh = 0.2,
    quiet = TRUE
  ) %>%
    distinct(group, gene, n.effective, div.class)
  expect_equal(get_result$n.effective, c(1L, 2L))
  expect_equal(get_result$div.class, c("monoform", "polyform"))

  get_shannon_top2 <- GetDiversity(
    object,
    genes = "gene1",
    entropy.use = "Shannon",
    top.n = 2,
    min.tx.cts = 3,
    quiet = TRUE
  ) %>%
    distinct(group, gene, div.class)
  expect_equal(get_shannon_top2$div.class, c("monoform", "polyform"))

  get_tsallis_one <- GetDiversity(
    object,
    genes = "gene1",
    entropy.use = "Tsallis",
    order = 1,
    min.tx.cts = 3,
    quiet = TRUE
  ) %>%
    distinct(group, gene, div.class)
  expect_true(all(is.na(get_tsallis_one$div.class)))

  get_tsallis_one_explicit <- GetDiversity(
    object,
    genes = "gene1",
    entropy.use = "Tsallis",
    order = 1,
    entropy.thresh = 0.5,
    min.tx.cts = 3,
    quiet = TRUE
  ) %>%
    distinct(group, gene, div.class)
  expect_equal(
    get_tsallis_one_explicit$div.class,
    c("monoform", "polyform")
  )

  plot_result <- PlotDiversity(
    object,
    genes = "gene1",
    min.tx.cts = 3,
    prop.thresh = 0.2,
    quiet = TRUE
  )$data %>%
    arrange(group_var)
  expect_equal(plot_result$n.effective, c(1L, 2L))
  expect_equal(plot_result$class, c("monoform", "polyform"))

  plot_shannon_top2 <- PlotDiversity(
    object,
    genes = "gene1",
    entropy.use = "Shannon",
    top.n = 2,
    min.tx.cts = 3,
    quiet = TRUE
  )$data %>%
    arrange(group_var)
  expect_equal(plot_shannon_top2$class, c("monoform", "polyform"))

  plot_tsallis_one <- PlotDiversity(
    object,
    genes = "gene1",
    entropy.use = "Tsallis",
    order = 1,
    min.tx.cts = 3,
    quiet = TRUE
  )$data
  expect_true(all(is.na(plot_tsallis_one$class)))

  expect_equal(Hypatia:::.EffectiveIsoformCount(rep(0.2, 5), 0.2), 5L)
  expect_equal(Hypatia:::.EffectiveIsoformCount(rep(1 / 6, 6), 0.2), 0L)
  expect_equal(Hypatia:::.EffectiveIsoformClass(0L), "polyform")
  expect_equal(Hypatia:::.EffectiveIsoformClass(1L), "monoform")
  expect_equal(Hypatia:::.EffectiveIsoformClass(2L), "polyform")

  expect_error(
    RunDIV(object, prop.thresh = 0, quiet = TRUE),
    "greater than 0"
  )
  expect_error(
    GetDiversity(object, genes = "gene1", prop.thresh = 1.1, quiet = TRUE),
    "not <= 1"
  )
  expect_error(
    GetDiversity(object, genes = "gene1", top.n = 1, quiet = TRUE),
    "at least 2"
  )
  expect_error(
    PlotDiversity(object, genes = "gene1", top.n = 2.5, quiet = TRUE),
    "type 'count'"
  )
  expect_error(
    RunDIV(object, renormalize = NA, quiet = TRUE),
    "May not be NA"
  )
  expect_error(
    GetDiversity(object, genes = "gene1", entropy.use = "Shannon2", quiet = TRUE),
    "element of set"
  )
  expect_error(
    GetDiversity(object, genes = "gene1", entropy.use = "Renyi2", quiet = TRUE),
    "element of set"
  )

  get_messages <- capture.output(
    invisible(GetDiversity(
      object,
      genes = "gene1",
      entropy.use = "Shannon",
      min.tx.cts = 3,
      quiet = FALSE
    )),
    type = "message"
  )
  expect_equal(sum(grepl("Using entropy.thresh", get_messages)), 1)
  expect_true(any(grepl(
    "Using entropy.thresh = 0.5 for monoform/polyform classification.",
    get_messages,
    fixed = TRUE
  )))

  unsupported_messages <- capture.output(
    invisible(GetDiversity(
      object,
      genes = "gene1",
      entropy.use = "Tsallis",
      order = 1,
      min.tx.cts = 3,
      quiet = FALSE
    )),
    type = "message"
  )
  expect_equal(sum(grepl("No default entropy.thresh", unsupported_messages)), 1)
  expect_true(any(grepl(
    "monoform/polyform classifications will be NA unless entropy.thresh is supplied",
    unsupported_messages,
    fixed = TRUE
  )))
  expect_length(
    capture.output(
      invisible(GetDiversity(object, genes = "gene1", min.tx.cts = 3, quiet = TRUE)),
      type = "message"
    ),
    0
  )

  plot_messages <- capture.output(
    invisible(PlotDiversity(
      object,
      genes = "gene1",
      entropy.thresh = 0.4,
      min.tx.cts = 3,
      quiet = FALSE
    )),
    type = "message"
  )
  expect_equal(sum(grepl("Using entropy.thresh", plot_messages)), 1)
  expect_true(any(grepl("entropy.thresh = 0.4", plot_messages, fixed = TRUE)))

  set.seed(1024)
  run_messages <- capture.output(
    invisible(RunDIV(
      object,
      group.1 = "A",
      group.2 = "B",
      min.gene.cts = 0,
      min.gene.pct = 0,
      min.tx.cts = 3,
      boot.iter = 1,
      boot.fraction = 1,
      p.adj = "none",
      quiet = FALSE
    )),
    type = "message"
  )
  expect_equal(sum(grepl("Using entropy.thresh", run_messages)), 1)
  expect_true(any(grepl("entropy.thresh = 0.243", run_messages, fixed = TRUE)))
})

test_that("entropy helper supports top-isoform selection for every metric", {
  prop <- c(0.5, 0.3, 0.2)
  selected <- c(0.5, 0.3)
  normalized <- selected / sum(selected)

  expect_equal(
    Hypatia:::.DiversityFunction("Shannon")(prop),
    -sum(prop * log(prop)),
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("Shannon", top.n = 2)(prop),
    -sum(selected * log(selected)),
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("Renyi")(prop),
    (1 / (1 - 2)) * log(sum(prop^2)),
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("Renyi", top.n = 2)(prop),
    -log(sum(selected^2)),
    tolerance = 1e-12
  )

  expected <- c(
    Shannon = -sum(selected * log(selected)),
    NormalizedShannon = -sum(selected * log(selected)) / log(2),
    Renyi = -log(sum(selected^2)),
    NormalizedRenyi = -log(sum(selected^2)) / log(2),
    GiniSimpson = 1 - sum(selected^2),
    Tsallis = (1 - sum(selected^3)) / 2,
    InverseSimpson = 1 / sum(selected^2)
  )
  expected_normalized <- c(
    Shannon = -sum(normalized * log(normalized)),
    NormalizedShannon = -sum(normalized * log(normalized)) / log(2),
    Renyi = -log(sum(normalized^2)),
    NormalizedRenyi = -log(sum(normalized^2)) / log(2),
    GiniSimpson = 1 - sum(normalized^2),
    Tsallis = (1 - sum(normalized^3)) / 2,
    InverseSimpson = 1 / sum(normalized^2)
  )

  for (metric in names(expected)) {
    expect_equal(
      Hypatia:::.DiversityFunction(metric, top.n = 2)(prop),
      unname(expected[[metric]]),
      tolerance = 1e-12
    )
    expect_equal(
      Hypatia:::.DiversityFunction(metric, top.n = 2, renormalize = TRUE)(prop),
      unname(expected_normalized[[metric]]),
      tolerance = 1e-12
    )
  }

  expect_equal(
    Hypatia:::.DiversityFunction("Shannon", top.n = 10)(prop),
    Hypatia:::.DiversityFunction("Shannon")(prop),
    tolerance = 1e-12
  )

  edge_prop <- c(0.1, 0.1, rep(0.08, 10))
  expect_equal(
    Hypatia:::.DiversityFunction("Shannon", top.n = 2)(edge_prop),
    -2 * 0.1 * log(0.1),
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("Shannon", top.n = 2, renormalize = TRUE)(edge_prop),
    log(2),
    tolerance = 1e-12
  )
})

test_that("normalized entropy returns NA when only one isoform is present", {
  prop <- c(1, 0, 0)

  expect_true(is.na(Hypatia:::.DiversityFunction("NormalizedShannon")(prop)))
  expect_true(is.na(Hypatia:::.DiversityFunction("NormalizedRenyi")(prop)))
  expect_true(is.na(Hypatia:::.DiversityFunction("NormalizedRenyi", order = 1)(prop)))
})

test_that("Tsallis and Renyi use their Shannon limits at order one", {
  prop <- c(0.5, 0.3, 0.2)
  shannon <- -sum(prop * log(prop))
  normalized_shannon <- shannon / log(length(prop))

  expect_equal(
    Hypatia:::.DiversityFunction("Tsallis", order = 1)(prop),
    shannon,
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("Renyi", order = 1)(prop),
    shannon,
    tolerance = 1e-12
  )
  expect_equal(
    Hypatia:::.DiversityFunction("NormalizedRenyi", order = 1)(prop),
    normalized_shannon,
    tolerance = 1e-12
  )
  expect_true(is.na(Hypatia:::.DiversityThreshold("Tsallis", order = 1)))
  expect_true(is.na(Hypatia:::.DiversityThreshold("Tsallis", order = 0.5)))
  expect_equal(Hypatia:::.DiversityThreshold("Tsallis", order = 3), 0.243)
  expect_true(is.na(Hypatia:::.DiversityThreshold("Renyi", order = 1)))
  expect_true(is.na(Hypatia:::.DiversityThreshold("Renyi", order = 3)))
  expect_equal(Hypatia:::.DiversityThreshold("Renyi", order = 2), 0.435)
  expect_equal(
    Hypatia:::.DiversityThreshold("Tsallis", entropy.thresh = 0.4, order = 1),
    0.4
  )
  expect_equal(Hypatia:::.DiversityThreshold("NormalizedRenyi", order = 1), 0)
})
