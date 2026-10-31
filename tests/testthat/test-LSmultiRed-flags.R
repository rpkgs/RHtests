test_that("regression uses the supplied continuity flags", {
  saved <- mget(c("Nt", "IY0flg"), envir = .GlobalEnv, ifnotfound = list(NULL))
  present <- vapply(c("Nt", "IY0flg"), exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    for (name in names(saved)) {
      if (present[name]) assign(name, saved[[name]], envir = .GlobalEnv)
      else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) rm(list = name, envir = .GlobalEnv)
    }
  })
  assign("Nt", 1, envir = .GlobalEnv)
  n <- 40
  y <- sin(seq_len(n) / 3) + seq_len(n) / 100
  flags <- c(rep(1, n - 1), 0)
  assign("IY0flg", rep(0, n), envir = .GlobalEnv)
  first <- RHtests:::LSmultipleRed(y, seq_len(n), n, flags)
  assign("IY0flg", rep(1, n), envir = .GlobalEnv)
  second <- RHtests:::LSmultipleRed(y, seq_len(n), n, flags)
  expect_equal(first, second)
  expect_error(RHtests:::LSmultipleRed(y, seq_len(n), n), "IY0flg")
})
