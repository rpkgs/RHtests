# 比较保留与移除候选突变点后的残差，计算 PMF 检验概率。
pmf_probability <- function(series, full_residuals, Ti, other_points, segment) {
  reduced_residuals <- LS_MultiBreaks(series, Ti, other_points)$resi
  full_sse <- sum(full_residuals[segment]^2)
  reduced_sse <- sum(reduced_residuals[segment]^2)
  df <- length(segment) - 3
  statistic <- (reduced_sse - full_sse) * df / full_sse
  if (statistic < 0) 0 else pf(statistic, 1, df)
}

# 逐个检验候选突变点，返回统计量与判定结果。
test_changepoints <- function(Y1, Ti, Ips, Ids, IY0, fit, plev, stepsizes, phi, PFmax) {
  Ns <- length(Ips) - 1
  RW <- LS_MultiBreaks(fit$W, Ti, Ips)$resi
  RWL <- LS_MultiBreaks(fit$WL, Ti, Ips)$resi
  RWU <- LS_MultiBreaks(fit$WU, Ti, Ips)$resi

  d_TP <- list()
  if (Ns > 0) {
    for (i in 1:Ns) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I3 <- Ips[i + 1]
      Ic <- Ips[i]
      Id <- Ids[i]
      Nseg <- I3 - I1 + 1

      PFx95 <- getPFx95(fit$cor, Nseg, phi, PFmax)
      PFx95l <- getPFx95(fit$corl, Nseg, phi, PFmax)
      PFx95h <- getPFx95(fit$corh, Nseg, phi, PFmax)
      SSEf.Iseg <- sum(fit$resi[I1:I3]^2)
      Ips0 <- Ips[-i]
      otmp <- LS_MultiBreaks(Y1, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      Pk0 <- Pk.PMFT(Nseg)
      PFx <- Fx * Pk0[Ic - I1 + 1]

      otmp <- LS_MultiBreaks(fit$W, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RW[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      if (Fx < 0) {
        PFx <- 0
        Fx <- 0
        prob <- 0
      } else {
        prob <- pf(Fx, 1, Nseg - 3)
      }

      # 在自相关系数的上下界分别检验，得到概率范围。
      probL0 <- pmf_probability(fit$WL, RWL, Ti, Ips0, I1:I3)
      probU0 <- pmf_probability(fit$WU, RWU, Ti, Ips0, I1:I3)
      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      Idc <- classify_changepoint(Id, probL, probU, plev, PFx, PFx95l, PFx95h)

      d_TP[[i]] <- data.table(
        kind = Id, Idc,
        date = IY0[Ic], Ic, Nseg, stepsize = stepsizes[i],
        probL, probU, plev, PFx, PFx95l, PFx95h, prob, PFx95
      )
    }
  }
  d_TP
}
# 在基准与参考序列的差值上逐个检验突变点。
test_changepoints_ref <- function(Y0, IY0, Ips, Ids, fit, plev, phi, PTmax) {
  Ns <- length(Ips) - 1
  d_TP <- list()
  if (Ns > 0) {
    for (i in 1:Ns) {
      Ic <- Ips[i]
      Id <- Ids[i]
      I0 <- if (i == 1) 0 else Ips[i - 1]
      I3 <- Ips[i + 1] # from t-1:t+1
      Nseg <- I3 - I0
      PTx95 <- getPTx95(fit$cor, Nseg, phi, PTmax)
      PTx95L <- getPTx95(fit$corl, Nseg, phi, PTmax)
      PTx95U <- getPTx95(fit$corh, Nseg, phi, PTmax)

      Pk0 <- Pk.PMT(Nseg)
      otmp <- PTKIc(fit$W[(I0 + 1):I3], Pk0, Ic - I0)
      prob <- otmp$prob
      otmp <- PTKIc(fit$WL[(I0 + 1):I3], Pk0, Ic - I0)
      probL0 <- otmp$prob
      otmp <- PTKIc(fit$WU[(I0 + 1):I3], Pk0, Ic - I0)
      probU0 <- otmp$prob
      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      otmp <- PTKIc(Y0[(I0 + 1):I3], Pk0, Ic - I0)
      PTx0 <- otmp$PTk
      Idc <- classify_changepoint(Id, probL, probU, plev, PTx0, PTx95L, PTx95U)

      d_TP[[i]] <- data.table(
        kind = Id, Idc,
        date = IY0[Ic], Ic, Nseg, stepsize = NA,
        probL, probU, plev, PFx = PTx0, PFx95l = PTx95L, PFx95h = PTx95U, prob, PFx95 = PTx95
      )
    }
  }
  d_TP
}
