test_that("reference fit updates the seasonal cycle until convergence", {
  calls <- 0L
  local_mocked_bindings(
    LSmultipleRed = function(...) {
      calls <<- calls + 1L
      list(coef = c(0, 0, 2))
    },
    .package = "RHtests"
  )

  result <- RHtests:::fit_reference_cycle(
    base = c(1, 2, 1, 2, 3, 4, 3, 4), seasonal = rep(c(1, 2), 4),
    previous_cycle = c(1, 2), cycle_index = rep(1:2, 4), cycle_values = 1:2,
    time = 1:8, points = c(4, 8), continuous = rep(1, 8), reference_segment = 1, n_cycles = 2
  )
  expect_equal(calls, 2L)
  expect_equal(result$cycle, c(2, 4))
  expect_equal(result$seasonal, rep(c(2, 4), 4))
  expect_equal(result$base, c(0, 0, 0, 0, 2, 2, 2, 2))
})

test_that("reference fit preserves the seasonal values without changepoints", {
  local_mocked_bindings(LSmultipleRed = function(...) list(coef = c(0, 0)), .package = "RHtests")
  result <- RHtests:::fit_reference_cycle(
    base = 1:4, seasonal = c(1, 2, 1, 2), previous_cycle = c(1, 2),
    cycle_index = c(1, 2, 1, 2), cycle_values = 1:2,
    time = 1:4, points = 4, continuous = rep(1, 4), reference_segment = 1, n_cycles = 2
  )
  expect_equal(result$base, 1:4)
  expect_equal(result$cycle, c(1, 2, 1, 2))
})
