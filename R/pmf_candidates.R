# 初始候选突变点的统计量，供 FindU() 和 FindUD() 共用。
pmf_candidates <- function(indices, Y1, Ti, N, IY0flg) {
  lapply(indices, function(I) {
    if (I <= 0) return(list(PFx = 0, prob = 0, corL = 0))

    Ips <- c(I, N)
    fit <- LSmultiple(Y1, Ti, Ips)
    red <- Rphi(fit$resi, Ips, 1, IY0flg)
    W <- red$W + fit$fitted
    list(
      PFx = PMFxKc(Y1, Ti, 0, N, I)$PFc,
      prob = PMFxKc(W, Ti, 0, N, I)$prob,
      corL = red$corl
    )
  })
}
