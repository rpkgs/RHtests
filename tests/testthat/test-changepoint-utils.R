test_that("adjustment segment follows the original selection rules", {
  Ips <- c(10, 20, 25)
  expect_equal(adjustment_segment(Ips, 2, 0), 1)
  expect_equal(adjustment_segment(Ips, 2, 2), 2)
  expect_equal(adjustment_segment(Ips, 2, 4), 3)
  expect_equal(adjustment_segment(25, 0, 0), 1)
  expect_equal(adjustment_segment(25, 0, 1), 1)
})

test_that("changepoint dates match exact or preceding observations", {
  dates <- c(20000101, 20000103, 20000103, 20000107)
  expect_equal(match_changepoint_dates(c(20000103, 20000105, 20000109), dates),
    c(2, 3, 4, 4))
  expect_equal(match_changepoint_dates(numeric(), dates), 4)
  expect_warning(result <- match_changepoint_dates(19991231, dates), "no non-missing")
  expect_equal(result, c(-Inf, 4))
  expect_true(is.na(match_changepoint_dates(NA_real_, dates)[1]))
})
