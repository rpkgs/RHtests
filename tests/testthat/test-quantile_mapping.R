test_that("quantile_mapping matches QMadjGaussian", {
  set.seed(1)
  breaks <- c(60, 120)
  y <- rnorm(180) + rep(c(2, 0, -1), each = 60)

  for (ref in 1:3) for (window in c(0, 24)) {
    expected <- QMadjGaussian(y, c(breaks, 180), 3, ref, window, Nt = 12)$PA
    expect_equal(quantile_mapping(y, breaks, ref, 3, window = window), expected)
  }
})

test_that("shared core uses reference differences and preserves tie ranks", {
  y <- rep(rep(1:10, each = 2), 2) + rep(c(2, 0), each = 20)
  dif <- rep(c(2, 0), each = 20)
  r <- RHtests:::.quantile_mapping(y, 20, 2, 2, 5, 0, dif)
  expect_equal(r$y, rep(rep(1:10, each = 2), 2))
  expect_equal(r$AdjM, c(-2, 0))
  expect_equal(r$p[1:2], c(2, 1) / 20)
  expect_true(all(is.na(r$W[21:40])))
})

test_that("quantile_mapping handles no breaks and short windows", {
  y <- c(seq_len(30) + 2, seq_len(30))
  expect_identical(quantile_mapping(y, numeric()), y)
  expect_identical(QMadjGaussian(y, length(y), 10, 1, 0, 12)$PA, y)

  r <- RHtests:::.quantile_mapping(y, 30, 2, 10, 5, 2)
  expect_equal(r$Mq, 1)
  expect_true(all(is.finite(r$y)))
  expect_equal(r$y[31:60], y[31:60])
})

test_that("quantile_mapping keeps the reference segment and aligns means", {
  y <- rep(seq_len(30), 2) + rep(c(2, 0), each = 30)
  r <- quantile_mapping(y, 30, ref = 2, Mq = 1)
  expect_equal(r, rep(seq_len(30), 2))
  expect_identical(r[31:60], y[31:60])
})
