test_that("fitted data uses named tibble columns and preserves missing rows", {
  original <- tibble::tibble(year = c(2001, 2001, 2001), month = 1, day = 1:3, data = c(10, NA, 12))
  fit <- list(Y0 = c(1, 2), meanhat = c(0.5, 1.5), trend = 0)
  result <- RHtests:::stepsize_data(
    original, c(TRUE, FALSE, TRUE), fit, 1, c(11, 13), c(1, 1), c(2, 3),
    0, NULL, c(11, 13), c(1, 2), c(1, 2), 2, c(20010101, 20010103), 1, 1, 0, 1
  )
  expect_s3_class(result, "tbl_df")
  expect_identical(names(result), RHtests:::fitdata_varnames_noref)
  expect_equal(result$base_mean_adj, c(11, NA, 13))
  expect_true(all(is.na(result$QM_adjusted)))
})

test_that("leap-day restoration preserves fitted values and adjustments", {
  original <- tibble::tibble(year = 2000, month = c(2, 2, 3), day = c(28, 29, 1), data = c(10, 12, 14))
  odata <- matrix(NA_real_, 3, 6)
  odata[, 2:3] <- cbind(c(20000228, 20000229, 20000301), original$data)
  odata[c(1, 3), 4:6] <- rbind(c(5, 8, 9), c(7, 11, 13))
  colnames(odata) <- c("id", "date", "base", "fitted", "adjusted", "QM_adjusted")
  odata <- tibble::as_tibble(odata)

  result <- RHtests:::stepsize_leap(
    odata, original, 20000301, NULL, NULL, NULL, NULL, 1, 1, 0, 1,
    fitted_cols = "fitted", adjusted_cols = "adjusted"
  )
  expect_equal(unlist(result[2, 4:6], use.names = FALSE), c(5, 10, 12))

  result <- RHtests:::stepsize_leap(
    odata, original, c(20000228, 20000301), c(10, 14), c(1, 2), 0,
    list(Mq = 1), 1, 1, 0, 1, fitted_cols = "fitted", adjusted_cols = "adjusted"
  )
  expect_equal(unlist(result[2, 4:5], use.names = FALSE), c(7, 9))
  expect_true(is.na(result$QM_adjusted[2]))
})
