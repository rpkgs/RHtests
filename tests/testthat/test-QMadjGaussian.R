test_that("QM results match the saved legacy results", {
  # 旧版固定基准：三段序列、中间基准段、局部窗口、同值及参考缺测
  expected <- readRDS(test_path("fixtures", "QMadjGaussian.rds"))
  fun <- QMadjGaussian
  environment(fun) <- list2env(list(Nt = 12), parent = environment(fun))

  actual <- fun(expected$y, expected$Ips, 3, 2, 24)
  expect_equal(actual, expected$noref, tolerance = 1e-12)

  actual <- QMadjGaussian.wRef(
    expected$y, expected$Dif, expected$Ips, 3, 2, 0, 12, 2
  )
  expect_equal(actual, expected$ref, tolerance = 1e-12)
})

test_that("QMadjGaussian adjusts to the chosen segment", {
  # 隔离 Nt，不修改全局环境
  fun <- QMadjGaussian
  environment(fun) <- list2env(list(Nt = 12), parent = environment(fun))
  y <- rep(seq_len(30), 2) + rep(c(2, 0), each = 30)

  for (Mq in c(-1, 0, 1, 3)) {
    for (ref in 1:2) {
      r <- fun(y, c(30, 60), Mq, ref, 0)
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
