make_output_contract_sce <- function() {
  count_data <- matrix(
    c(
      10, 9, 8, 7, 2, 2, 1, 1, 5, 5, 4, 4,
      1, 1, 2, 2, 8, 7, 9, 8, 5, 4, 5, 4,
      6, 5, 6, 5, 3, 4, 3, 4, 9, 8, 9, 8,
      2, 3, 2, 3, 7, 6, 7, 6, 1, 2, 1, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("tx", seq_len(4)), paste0("cell", seq_len(12)))
  )
  col_data <- data.frame(
    group = rep(c("A", "B", "C"), each = 4),
    row.names = colnames(count_data)
  )
  row_data <- data.frame(
    gene_id = rep(c("gene1", "gene2"), each = 2),
    row.names = rownames(count_data)
  )

  object <- CreateSCE(
    count_data,
    col_data,
    row_data,
    active.group.id = "group",
    quiet = TRUE
  )
  assay(object, "logcounts") <- log1p(count_data)
  object
}

expect_output_schema <- function(object, schema) {
  expect_s3_class(object, "data.frame")
  expect_identical(names(object), names(schema))
  expect_identical(
    vapply(object, typeof, character(1)),
    schema
  )
}

expect_rows_ordered_by <- function(object, columns) {
  expected_order <- do.call(
    order,
    c(as.list(object[columns]), list(na.last = TRUE))
  )
  expect_identical(expected_order, seq_len(nrow(object)))
}

usage_schema <- c(
  group = "character", gene = "character", gene.pct = "double",
  transcript = "character", cts = "double", prop = "double",
  cell.n = "integer", cell.prop.mean = "double",
  cell.prop.median = "double", cell.prop.sd = "double",
  cell.prop.iqr = "double", cell.prop.zero.frac = "double"
)

diversity_schema <- c(
  group = "character", gene = "character", gene.pct = "double",
  n.transcripts = "integer", n.effective = "integer",
  transcript = "character", cts = "double", prop = "double",
  div = "double", div.class = "character", cell.n = "integer",
  cell.div.mean = "double", cell.div.median = "double",
  cell.div.sd = "double", cell.div.iqr = "double"
)

expression_schema <- c(
  group = "character", gene = "character", transcript = "character",
  pct = "double", avgExpr = "double", cell.n = "integer",
  cell.expr.median = "double", cell.expr.sd = "double",
  cell.expr.iqr = "double"
)

diu_data_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  gene.pct.1 = "double", gene.pct.2 = "double",
  transcript = "character", cts.1 = "double", cts.2 = "double",
  prop.1 = "double", prop.2 = "double", prop.diff = "double",
  boot.prop.1 = "list", boot.prop.2 = "list", boot.prop.diff = "list",
  cell.n.1 = "integer", cell.prop.mean.1 = "double",
  cell.prop.median.1 = "double", cell.prop.sd.1 = "double",
  cell.prop.iqr.1 = "double", cell.prop.zero.frac.1 = "double",
  cell.n.2 = "integer", cell.prop.mean.2 = "double",
  cell.prop.median.2 = "double", cell.prop.sd.2 = "double",
  cell.prop.iqr.2 = "double", cell.prop.zero.frac.2 = "double"
)

diu_stats_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  max.prop.diff = "double", transcript = "character", pval = "double",
  padj = "double", pval.perm = "double", padj.perm = "double",
  pval.perm.delta = "double", padj.perm.delta = "double",
  cramers.v = "double", approx = "character",
  boot.prop.diff.mean = "double", boot.prop.diff.lower = "double",
  boot.prop.diff.upper = "double", boot.valid.iter = "integer"
)

div_data_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  gene.pct.1 = "double", gene.pct.2 = "double",
  n.transcripts = "integer", div.1 = "list", div.2 = "list",
  cell.n.1 = "integer", cell.div.mean.1 = "double",
  cell.div.median.1 = "double", cell.div.sd.1 = "double",
  cell.div.iqr.1 = "double", cell.n.2 = "integer",
  cell.div.mean.2 = "double", cell.div.median.2 = "double",
  cell.div.sd.2 = "double", cell.div.iqr.2 = "double"
)

div_stats_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  avgDiv.1 = "double", avgDiv.2 = "double", div.diff = "double",
  n.effective.1 = "integer",
  n.effective.2 = "integer", div.class.1 = "character",
  div.class.2 = "character", div.diff.lower = "double",
  div.diff.upper = "double", support.positive = "double",
  support.negative = "double", boot.valid.iter = "integer",
  supported = "logical"
)

dei_data_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  transcript = "character", pct.1 = "double", pct.2 = "double",
  avgExpr.1 = "double", avgExpr.2 = "double", cell.n.1 = "integer",
  cell.expr.median.1 = "double", cell.expr.sd.1 = "double",
  cell.expr.iqr.1 = "double", cell.n.2 = "integer",
  cell.expr.median.2 = "double", cell.expr.sd.2 = "double",
  cell.expr.iqr.2 = "double"
)

dei_stats_schema <- c(
  group.1 = "character", group.2 = "character", gene = "character",
  transcript = "character", log2FC = "double", pval = "double",
  padj = "double"
)

test_that("Get function schemas and empty-result types are stable", {
  object <- make_output_contract_sce()

  usage <- GetUsage(
    object, c("gene1", "gene2"), min.tx.cts = 0, quiet = TRUE
  )
  diversity <- GetDiversity(
    object, c("gene1", "gene2"), min.tx.cts = 0, quiet = TRUE
  )
  expression <- GetExpression(
    object, transcripts = rownames(object), quiet = TRUE
  )
  empty_usage <- GetUsage(
    object, c("gene1", "gene2"), min.tx.cts = 1e9, quiet = TRUE
  )
  empty_diversity <- GetDiversity(
    object, c("gene1", "gene2"), min.tx.cts = 1e9, quiet = TRUE
  )

  expect_output_schema(usage, usage_schema)
  expect_output_schema(diversity, diversity_schema)
  expect_output_schema(expression, expression_schema)
  expect_output_schema(empty_usage, usage_schema)
  expect_output_schema(empty_diversity, diversity_schema)
  expect_equal(nrow(empty_usage), 0L)
  expect_equal(nrow(empty_diversity), 0L)
})

test_that("GetExpression declares transcripts as optional in gene mode", {
  object <- make_output_contract_sce()

  expect_identical(formals(GetExpression)$transcripts, NULL)
  result <- GetExpression(object, genes = "gene1", quiet = TRUE)

  expect_setequal(result$transcript, c("tx1", "tx2"))
  expect_output_schema(result, expression_schema)
})

test_that("Run function schemas, list columns, and ordering are stable", {
  object <- make_output_contract_sce()

  diu <- RunDIU(
    object,
    min.gene.pct = 0,
    min.gene.cts = 0,
    min.tx.cts = 0,
    bootstrap = FALSE,
    quiet = TRUE
  )
  set.seed(1024)
  div <- RunDIV(
    object,
    min.gene.pct = 0,
    min.gene.cts = 0,
    min.tx.cts = 0,
    boot.iter = 3,
    boot.fraction = 1,
    quiet = TRUE
  )
  dei <- suppressWarnings(RunDEI(object, min.pct = 0, quiet = TRUE))

  expect_named(diu, c("data", "stats"))
  expect_named(div, c("data", "stats"))
  expect_named(dei, c("data", "stats"))
  expect_output_schema(diu$data, diu_data_schema)
  expect_output_schema(diu$stats, diu_stats_schema)
  expect_output_schema(div$data, div_data_schema)
  expect_output_schema(div$stats, div_stats_schema)
  expect_output_schema(dei$data, dei_data_schema)
  expect_output_schema(dei$stats, dei_stats_schema)

  expect_s3_class(div$data$div.1, "AsIs")
  expect_s3_class(div$data$div.2, "AsIs")
  expect_true(all(lengths(div$data$div.1) == 3L))
  expect_true(all(lengths(div$data$div.2) == 3L))

  expect_rows_ordered_by(
    diu$data, c("group.1", "group.2", "gene", "transcript")
  )
  expect_rows_ordered_by(
    diu$stats, c("group.1", "group.2", "padj", "gene", "transcript")
  )
  expect_rows_ordered_by(div$data, c("group.1", "group.2", "gene"))
  expect_rows_ordered_by(
    div$stats, c("group.1", "group.2", "gene")
  )
  expect_rows_ordered_by(
    dei$data, c("group.1", "group.2", "gene", "transcript")
  )
  expect_rows_ordered_by(
    dei$stats, c("group.1", "group.2", "padj", "gene", "transcript")
  )
})

test_that("RunDIU keeps output types across methods and reports signed maxima", {
  object <- make_output_contract_sce()

  set.seed(1024)
  fisher <- RunDIU(
    object,
    method.use = "Fisher",
    min.gene.pct = 0,
    min.gene.cts = 0,
    min.tx.cts = 0,
    bootstrap = FALSE,
    quiet = TRUE
  )

  expect_output_schema(fisher$stats, diu_stats_schema)
  expect_true(all(is.na(fisher$stats$approx)))

  for (i in seq_len(nrow(fisher$stats))) {
    stat_row <- fisher$stats[i, , drop = FALSE]
    data_rows <- fisher$data[
      fisher$data$group.1 == stat_row$group.1 &
        fisher$data$group.2 == stat_row$group.2 &
        fisher$data$gene == stat_row$gene,
      ,
      drop = FALSE
    ]
    selected_difference <- data_rows$prop.diff[
      match(stat_row$transcript, data_rows$transcript)
    ]

    expect_equal(stat_row$max.prop.diff, selected_difference)
    expect_equal(abs(stat_row$max.prop.diff), max(abs(data_rows$prop.diff)))
  }
})
