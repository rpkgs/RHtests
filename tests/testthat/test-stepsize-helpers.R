test_that("date columns convert to YYYYMMDD", {
  dates <- tibble::tibble(year = c(2000, 2001), month = c(2, 12), day = c(29, 1))
  expect_equal(RHtests:::.make_date(dates), c(20000229, 20011201))
})

test_that("changepoint dates and adjustment segments are shared", {
  dates <- c(20000101, 20000103, 20000105, 20000107, 20000109)
  points <- RHtests:::stepsize_points(
    data.frame(date = as.Date(c("2000-01-02", "2000-01-05")), kind = c(0, 1)), dates
  )
  expect_equal(points$Ips, c(1, 3, 5))
  expect_equal(points$Ids, c(0, 1))
  expect_equal(RHtests:::stepsize_bounds(points$Ips), list(start = c(1, 2, 4), end = c(1, 3, 5)))
  expect_equal(RHtests:::stepsize_bounds(5), list(start = 1, end = 5))
  expect_equal(RHtests:::stepsize_segment(points$Ips, 0), 2)
  expect_equal(RHtests:::stepsize_segment(points$Ips, 10000), 3)
  expect_equal(RHtests:::stepsize_segment(points$Ips, 1), 1)
  expect_equal(RHtests:::stepsize_segment(5, 0), 1)
})

test_that("changepoint classification respects confidence boundaries", {
  classify <- RHtests:::stepsize_classify
  expect_identical(classify(0, 0.8, 0.9, 0.95, 2, 3, 4), "No  ")
  expect_identical(classify(0, 0.9, 0.95, 0.95, 2, 3, 4), "?   ")
  expect_identical(classify(0, 0.95, 1, 0.95, 2, 3, 4), "YifD")
  expect_identical(classify(0, 0.8, 0.9, 0.95, 4, 3, 4), "Yes ")
  expect_identical(classify(1, 0, 0, 0.95, 2, 3, 4), "No  ")
  expect_identical(classify(1, 0, 0, 0.95, 3, 3, 4), "?   ")
  expect_identical(classify(1, 0, 0, 0.95, 4, 3, 4), "Yes ")
})
