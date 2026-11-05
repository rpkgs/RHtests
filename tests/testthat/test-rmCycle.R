test_that("RemoveSeasonal removes the mean for each month-day", {
  x <- c(5, 2, NA, 4, 8)
  month_day <- c(201, 101, 201, 101, 201)
  result <- RHtests:::RemoveSeasonal(y = x, t_mmdd = month_day)

  expect_equal(result$mu, c(3, 6.5))
  expect_equal(result$anomaly, c(-1.5, -1, NA, 1, 1.5))
  expect_equal(RHtests:::RemoveSeasonal(integer(), integer()), list(mu = numeric(), anomaly = numeric()))
})

test_that("seasonal regression returns a named seasonal mean", {
  state <- Read(system.file("extdata/Example1.dat", package = "RHtests"))
  series <- state@series
  season <- state@season
  initial <- RHtests:::RemoveSeasonal(series$Y0, series$Imd)
  fit <- RHtests:::LSmultiRedCycle(
    initial$anomaly, series$Ti, length(series$Y0), 1, series$IY0flg,
    seasonal_mean = initial$mu, series$Imd, season$Icy, season$Nt
  )
  expect_length(fit$seasonal_mean, season$Nt)
  expect_null(fit$EB)
})
