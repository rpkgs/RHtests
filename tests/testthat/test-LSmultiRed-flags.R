test_that("regression uses explicit continuity flags and cycle length", {
  n <- 40
  y <- sin(seq_len(n) / 3) + seq_len(n) / 100
  flags <- c(rep(1, n - 1), 0)
  first <- RHtests:::LSmultipleRed(y, seq_len(n), n, flags, Nt = 1)
  flags[10] <- 0
  second <- RHtests:::LSmultipleRed(y, seq_len(n), n, flags, Nt = 1)
  expect_false(isTRUE(all.equal(first, second)))
  expect_error(RHtests:::LSmultipleRed(y, seq_len(n), n, Nt = 1), "IY0flg")
})
