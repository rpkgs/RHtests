test_that("QM results match the saved legacy results", {
  # 旧版固定基准：三段序列、中间基准段、局部窗口、同值及参考缺测
  expected <- readRDS(test_path("fixtures", "QMadjGaussian.rds"))
  actual <- QMadjGaussian(expected$y, expected$Ips, 3, 2, 24, Nt = 12)
  expect_equal(actual, expected$noref, tolerance = 1e-12)

  actual <- QMadjGaussian.wRef(
    expected$y, expected$Dif, expected$Ips, 3, 2, 0, 12, 2
  )
  expect_equal(actual, expected$ref, tolerance = 1e-12)
})

test_that("QMadjGaussian adjusts to the chosen segment", {
  y <- rep(seq_len(30), 2) + rep(c(2, 0), each = 30)

  for (Mq in c(-1, 0, 1, 3)) {
    for (ref in 1:2) {
      r <- QMadjGaussian(y, c(30, 60), Mq, ref, 0, Nt = 12)
      expected <- rep(seq_len(30), 2) + if (ref == 1) 2 else 0
      expect_equal(r$PA, expected)
      expect_equal(r$PA, y + replace(r$W, is.na(r$W), 0))
      expect_null(dimnames(r$osp))
      expect_equal(dim(r$osmean), c(r$Mq + 2, 5))
    }
  }
})

test_that("QMadjGaussian.wRef preserves the reference segment", {
  ref <- rep(seq_len(30), 2)
  y <- ref + rep(c(2, 0), each = 30)

  for (Mq in c(-1, 0, 1, 3)) {
    r <- QMadjGaussian.wRef(y, y - ref, c(30, 60), Mq, 2, 999, 12, 0)
    expect_equal(r$PA, ref)
    expect_identical(r$PA[31:60], y[31:60])
    expect_equal(r$AdjM, c(-2, 0))
    expect_true(all(is.na(r$DB)))
    expect_null(dimnames(r$osp))
  }
})

test_that("QM does not depend on global variables", {
  y <- rep(seq_len(20), 2) + rep(c(2, 0), each = 20)
  difference <- rep(NA_real_, 40)
  difference[c(1:4, 21:24)] <- rep(c(2, 0), each = 4)
  globals <- mget(c("Nt", "ErrorMSG"), envir = .GlobalEnv,
    ifnotfound = list(NULL), inherits = FALSE)

  r <- QMadjGaussian(y, c(20, 40), 3, 2, 0, Nt = 12)
  expect_equal(r$Mq, 3)
  expect_warning(
    QMadjGaussian.wRef(y, difference, c(20, 40), 1, 2, 0, 12, 0),
    "参考序列有效数据不足"
  )
  expect_identical(mget(c("Nt", "ErrorMSG"), envir = .GlobalEnv,
    ifnotfound = list(NULL), inherits = FALSE), globals)
})

test_that("QMadjGaussian.wRef reduces categories for missing references", {
  ref <- rep(seq_len(60), 2)
  y <- ref + rep(c(2, 0), each = 60)
  Dif <- y - ref
  Dif[1:40] <- NA

  r <- QMadjGaussian.wRef(y, Dif, c(60, 120), 3, 2, 0, 12, 0)
  expect_equal(r$Mq, 1)
  expect_equal(r$Nseg.mn, 20)
  expect_equal(r$PA, ref)
})
