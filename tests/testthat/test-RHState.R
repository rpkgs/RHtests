test_that("readers return validated S7 states", {
  base <- system.file("extdata/Example1.dat", package = "RHtests")
  state <- Read(base)
  expect_true(S7::S7_inherits(state, RHtests:::RHState))
  expect_false(state@reference)
  expect_equal(length(state@series$Y0), length(state@series$IY0))
  expect_equal(length(state@season$Icy), state@season$Nt)

  base <- system.file("extdata/Example2.dat", package = "RHtests")
  reference <- system.file("extdata/Example2_Ref.dat", package = "RHtests")
  state_ref <- Read.wRef(base, reference)
  expect_true(S7::S7_inherits(state_ref, RHtests:::RHState))
  expect_true(state_ref@reference)

  points <- data.table::data.table(date = 19740201, kind = 1)
  expect_error(StepSize(InCs = points, output = tempfile(), .state = state_ref), "without a reference series")
  expect_error(StepSize.wRef(InCs = points, output = tempfile(), .state = state), "with a reference series")

  state_90 <- Read(system.file("extdata/Example1.dat", package = "RHtests"), plev = 0.9)
  expect_error(FindU(output = tempfile(), .state = state_90), "confidence level")
  output <- tempfile()
  on.exit(unlink(Sys.glob(paste0(output, "*"))))
  result <- RHtests_stepsize(TP2 = points, prefix = output, .state = state_90, verbose = FALSE)
  expect_equal(result$TP$plev, 0.9)
})

test_that("state fields are unpacked by name, not position", {
  base <- system.file("extdata/Example1.dat", package = "RHtests")
  state <- Read(base)
  reordered <- RHtests:::RHState(
    series = state@series[rev(names(state@series))],
    season = state@season[rev(names(state@season))],
    critical_values = state@critical_values[rev(names(state@critical_values))],
    reference = FALSE, plev = state@plev
  )
  points <- data.table::data.table(date = 19740201, kind = 1)
  output <- tempfile()
  on.exit(unlink(Sys.glob(paste0(output, "*"))))
  original <- StepSize(InCs = points, output = output, is_plot = FALSE, .state = state)
  actual <- StepSize(InCs = points, output = output, is_plot = FALSE, .state = reordered)
  expect_equal(actual, original)
})

test_that("S7 state rejects inconsistent inputs", {
  base <- system.file("extdata/Example1.dat", package = "RHtests")
  state <- Read(base)
  series <- state@series
  season <- state@season
  values <- state@critical_values

  series$IY0 <- series$IY0[-1]
  expect_error(RHtests:::RHState(series, season, values, FALSE, 0.95), "same nonzero length")
  series <- state@series
  season$Nt <- season$Nt + 1
  expect_error(RHtests:::RHState(series, season, values, FALSE, 0.95), "Season length")
  season <- state@season
  values$PFmax <- matrix(0, 1, 2)
  expect_error(RHtests:::RHState(series, season, values, FALSE, 0.95), "Critical value table")
})
