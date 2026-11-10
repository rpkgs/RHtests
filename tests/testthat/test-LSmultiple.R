test_that("LS_MultiBreaks fits the supplied predictor", {
  x <- seq_len(8)
  y <- c(2, 1, 4, 3, 5, 3, 8, 4)
  result <- RHtests:::LS_MultiBreaks(Y = y, X = x, Ips = c(4, 8))
  expected <- lm(y ~ x + I(x > 4))

  expect_equal(as.vector(result$coef), unname(coef(lm(y ~ I(x - mean(x)) + I(x > 4)))))
  expect_null(result$sig)
  expect_equal(as.vector(result$fitted), as.vector(fitted(expected)))
  expect_equal(as.vector(result$resi), as.vector(residuals(expected)))

  single <- LS_OneBreak(y, x - mean(x), 4)
  expect_equal(as.vector(single$coef), as.vector(result$coef))
  expect_equal(as.vector(single$fitted), as.vector(fitted(expected)))
  expect_equal(LS_OneBreak(y, x - mean(x), 4, only.SSE = TRUE)$SSE, sum(residuals(expected)^2))
})
