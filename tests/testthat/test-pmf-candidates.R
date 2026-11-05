test_that("initial candidate statistics match the original calculation", {
  set.seed(42)
  N <- 120
  Ti <- seq_len(N)
  Y1 <- rnorm(N) + (Ti > 60)
  IY0flg <- rep(1, N)
  indices <- c(0, 40, 75)

  expected <- lapply(indices, function(I) {
    if (I <= 0) return(list(PFx = 0, prob = 0, corL = 0))
    Ips <- c(I, N)
    fit <- LS_MultiBreaks(Y1, Ti, Ips)
    red <- Rphi(fit$resi, Ips, 1, IY0flg)
    raw <- PMFxKc(Y1, Ti, 0, N, I)
    adjusted <- PMFxKc(red$W + fit$fitted, Ti, 0, N, I)
    list(PFx = raw$PFc, prob = adjusted$prob, corL = red$corl)
  })

  expect_equal(pmf_candidates(indices, Y1, Ti, N, IY0flg), expected)
})
