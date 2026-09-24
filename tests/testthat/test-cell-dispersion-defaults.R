test_that("cell dispersion is opt-in for all public functions", {
  functions <- list(
    GetUsage,
    GetDiversity,
    GetExpression,
    RunDIU,
    RunDIV,
    RunDEI
  )

  expect_true(all(vapply(
    functions,
    function(fun) identical(formals(fun)$cell.dispersion, FALSE),
    logical(1)
  )))
})
