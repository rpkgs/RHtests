test_that("date columns convert to YYYYMMDD", {
  dates <- tibble::tibble(year = c(2000, 2001), month = c(2, 12), day = c(29, 1))
  expect_equal(RHtests:::ymd2num(dates), c(20000229, 20011201))
})

test_that("changepoint dates and adjustment segments are shared", {
  dates <- c(20000101, 20000103, 20000105, 20000107, 20000109)
  points <- RHtests:::locate_changepoints(
    data.frame(date = as.Date(c("2000-01-02", "2000-01-05")), kind = c(0, 1)), dates
  )
  expect_equal(points$Ips, c(1, 3, 5))
  expect_equal(points$Ids, c(0, 1))
  expect_equal(RHtests:::segment_bounds(points$Ips), list(start = c(1, 2, 4), end = c(1, 3, 5)))
  expect_equal(RHtests:::segment_bounds(5), list(start = 1, end = 5))
  expect_equal(RHtests:::select_base_segment(points$Ips, 0), 2)
  expect_equal(RHtests:::select_base_segment(points$Ips, 10000), 3)
  expect_equal(RHtests:::select_base_segment(points$Ips, 1), 1)
  expect_equal(RHtests:::select_base_segment(5, 0), 1)
})

test_that("segment adjustments apply to each observation", {
  segments <- RHtests:::segment_bounds(c(2, 5, 6))
  expect_equal(RHtests:::add_segment_adjustment(1:6, segments, c(1, -2, 0)), c(2, 3, 1, 2, 3, 6))
  expect_equal(RHtests:::add_segment_adjustment(1:3, RHtests:::segment_bounds(3), 0), 1:3)
})

test_that("PMF probability compares reduced and full residuals", {
  time <- 1:12
  series <- c(1, 3, 2, 4, 3, 5, 8, 7, 9, 8, 10, 9)
  full_residuals <- RHtests:::LS_MultiBreaks(series, time, c(6, 12))$resi
  reduced_residuals <- RHtests:::LS_MultiBreaks(series, time, 12)$resi
  full_sse <- sum(full_residuals^2)
  statistic <- (sum(reduced_residuals^2) - full_sse) * 9 / full_sse
  expect_equal(
    RHtests:::pmf_probability(series, full_residuals, time, 12, time),
    pf(statistic, 1, 9)
  )
  expect_equal(RHtests:::pmf_probability(series, rep(100, 12), time, 12, time), 0)
})

test_that("changepoint classification respects confidence boundaries", {
  classify <- RHtests:::classify_changepoint
  expect_identical(classify(0, 0.8, 0.9, 0.95, 2, 3, 4), "No  ")
  expect_identical(classify(0, 0.9, 0.95, 0.95, 2, 3, 4), "?   ")
  expect_identical(classify(0, 0.95, 1, 0.95, 2, 3, 4), "YifD")
  expect_identical(classify(0, 0.8, 0.9, 0.95, 4, 3, 4), "Yes ")
  expect_identical(classify(1, 0, 0, 0.95, 2, 3, 4), "No  ")
  expect_identical(classify(1, 0, 0, 0.95, 3, 3, 4), "?   ")
  expect_identical(classify(1, 0, 0, 0.95, 4, 3, 4), "Yes ")
})
