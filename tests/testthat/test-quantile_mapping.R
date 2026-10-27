test_that("quantile_mapping matches QMadjGaussian", {
  set.seed(1)
  breaks <- c(60, 120)
  y <- rnorm(180) + rep(c(2, 0, -1), each = 60)

  for (ref in 1:3) for (window in c(0, 24)) {
    expected <- QMadjGaussian(y, c(breaks, 180), 3, ref, window, Nt = 12)$PA
    expect_equal(quantile_mapping(y, breaks, ref, 3, window = window), expected)
  }
})

test_that("quantile_mapping keeps the reference segment and aligns means", {
  y <- rep(seq_len(30), 2) + rep(c(2, 0), each = 30)
  r <- quantile_mapping(y, 30, ref = 2, Mq = 1)
  expect_equal(r, rep(seq_len(30), 2))
  expect_identical(r[31:60], y[31:60])
})
