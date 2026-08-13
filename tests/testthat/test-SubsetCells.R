test_that("SubsetCells works", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_no_error({
    gbm_nCount <- SubsetCells(gbm, subset = nCount >= 1)
    gbm_nCount_inv <- SubsetCells(gbm, subset = nCount >= 1, invert = TRUE)
    gbm_celltype <- SubsetCells(gbm, subset = cell_type == "Oligodendrocyte")
    gbm_cells <- SubsetCells(gbm, cells = colnames(gbm)[1:5])
  })

  expect_true({
    all(colData(gbm_nCount)$nCount >= 1)
    all(colData(gbm_nCount_inv)$nCount < 1)
    all(colData(gbm_celltype)$cell_type == "Oligodendrocyte")
    all(colnames(gbm_cells) == colnames(gbm)[1:5])
  })

})

test_that("SubsetCells requires exactly one selector", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_error(SubsetCells(gbm, quiet = TRUE), "exactly one")
  expect_error(
    SubsetCells(gbm, subset = nCount >= 1, cells = colnames(gbm)[1:5], quiet = TRUE),
    "exactly one"
  )
})

test_that("SubsetCells supports non-syntactic colData names", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  colData(gbm)[["cell type"]] <- colData(gbm)$cell_type

  expect_no_error({
    gbm_celltype <- SubsetCells(gbm, subset = `cell type` == "Oligodendrocyte", quiet = TRUE)
  })
  expect_true(all(colData(gbm_celltype)[["cell type"]] == "Oligodendrocyte"))
})

test_that("SubsetCells reports invalid subset expressions", {

  gbm <-
    CreateSCE(
      countData = gbm_countData,
      colData = gbm_colData,
      rowData = gbm_rowData,
      quiet = TRUE
    )

  expect_error(
    SubsetCells(gbm, subset = missing_column == "x", quiet = TRUE),
    "subset term was not valid"
  )
})
