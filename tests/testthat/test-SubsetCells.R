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
