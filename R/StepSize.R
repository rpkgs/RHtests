#' StepSize
#'
#' @inheritParams FindU
#' @param InCs 待检验的突变点，data.frame 或文件路径，含列 `date` 与 `kind`。
#' @param debug
#' - `true`  : `_Fstat.txt` will be shown
#' - `false` : no `_Fstat.txt`
#'
#' @example R/examples/run_ex01.R
#' @export
StepSize <- function(
  InSeries = NULL, InCs, output,
  MissingValueCode = "-999.99",
  GUI = FALSE, plev = 0.95, Iadj = 10000, Mq = 10, Ny4a = 0,
  is_plot = TRUE, ..., debug = TRUE, .state = NULL
) {
  # 读取计算状态，或使用上游传入的状态。
  state <- if (is.null(.state)) {
    Read(InSeries, MissingValueCode, plev)
  } else .state
  check_rh_state(state, reference = FALSE, plev = plev)
  c(Y0, Ti, IY0, IY0flg, Imd, itable, ori.itable, ooflg) %<-%
    state@series[c("Y0", "Ti", "IY0", "IY0flg", "Imd", "itable", "ori.itable", "ooflg")]
  c(Icy, Nt) %<-% state@season[c("Icy", "Nt")]
  c(phi, PFmax) %<-% state@critical_values[c("phi", "PFmax")]

  # 将候选突变日期映射到观测位置。
  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Ny4a * Nt

  c(Ips, Ids, Ns) %<-% locate_changepoints(InCs, IY0)[c("Ips", "Ids", "Ns")]
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    message <- paste("StepSize: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
    stop(message, call. = FALSE)
  }
  # 去除季节周期，分别拟合无突变点和含突变点的趋势。
  cycle <- rmCycle(itable)
  Y1 <- cycle$Base
  EB <- cycle$EB
  if (length(EB) != length(Icy)) {
    stop("Annual cycle length differs from original dataset", call. = FALSE)
  }

  Iseg.adj <- select_base_segment(Ips, Iadj)
  segments <- segment_bounds(Ips)

  Ip0 <- N
  fit_base <- LSmultiRedCycle(Y1, Ti, Ip0, 1, IY0flg, EB, Imd, Icy, Nt)
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  fit <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj, IY0flg, EB, Imd, Icy, Nt)
  Y1 <- fit$Y0
  cor <- fit$cor
  corl <- fit$corl
  corh <- fit$corh
  pcor <- pt(abs(cor) * sqrt((N - 2) / (1 - cor^2)), N - 2)
  EB1 <- fit$EB

  EBfull <- EB1[match(Imd, Icy)]
  EEB <- mean(EB1, na.rm = T)

  # 先去除趋势，再按基准段进行分位数匹配，最后加回趋势。
  Rb <- Y1 - fit$trend * Ti + EBfull
  QMout <- NULL
  if (Ns > 0) {
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj, Nt)
    B <- QMout$PA
  } else {
    B <- Rb
  }

  adj <- fit$Y0 + EBfull
  B <- B + fit$trend * Ti

  ofilePdf <- paste0(output, "_F.pdf")
  ofileSout <- paste0(output, "_Fstat.txt")
  ofileIout <- paste0(output, "_fCs.txt")
  file.create(ofileSout)
  file.create(ofileIout)
  file.create(ofilePdf)

  if (debug) {
    log_fit_summary(
      ofileSout, plev, N, fit_base, Ns, QMout, Ny4a,
      Iseg.adj, Ips, fit, cor, corl, corh, pcor
    )
  }

  # 相邻段均值之差为突变幅度，各段均值向基准段订正。
  stepsizes <- diff(fit$mu)
  if (Ns > 0) {
    adj <- add_segment_adjustment(adj, segments, fit$mu[Iseg.adj] - fit$mu)
    if (debug) {
      for (i in seq_len(Ns)) {
        fprintf("%d %d stepsize= %.4f \n", Ips[i], IY0[Ips[i]], stepsizes[i],
          file = ofileSout, append = TRUE
        )
      }
    }
  }
  Ehat <- mean(fit$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  if (is_plot) {
    plot_stepsize(
      cycle, ofilePdf, EBfull, EEB, B, QMout, Mq, Ns, adj,
      Ips, Iseg.adj, fit, ori.itable, ooflg
    )
  }

  # 恢复原始日期和缺测位置，检验各候选突变点。
  odata <- build_fitted_table(
    ori.itable, ooflg, fit, EEB, adj, EBfull, meanhat0,
    Ns, QMout, B, Y1, Ti, Ips, IY0, Mq, Iseg.adj, Nadj, Nt
  )

  d_TP <- test_changepoints(Y1, Ti, Ips, Ids, IY0, fit, plev, stepsizes, phi, PFmax)
  if (debug) write_test_log(d_TP, "PMF", ofileSout)
  d_TP <- write_points_table(d_TP, ofileIout)

  odata <- dt_num2date(odata)
  # write.table(file = paste0(output, "_F.dat"), odata, na = MissingValueCode, col.names = TRUE, row.names = F)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from=ofileIout, to=ofileMout,overwrite=TRUE)
    # cat("StepSize finished successfully...\n")
    list(data = odata, TP = d_TP) # fit = odata,
  }
}

# 写出趋势拟合与订正结果摘要。
log_fit_summary <- function(ofileSout, plev, N, fit_base, Ns, QMout, Ny4a,
                         Iseg.adj, Ips, otmp, cor, corl, corh, pcor) {
  fprintf("The nominal level of confidence (1-alpha)= %.4f \n", plev, file = ofileSout)
  fprintf("Input data filename: N= %d \n", N, file = ofileSout, append = TRUE)
  fprintf(
    " Ignore changepoints -> trend0 = %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f )\n",
    fit_base$trend, fit_base$betaL, fit_base$betaU, fit_base$p.tr,
    fit_base$cor, fit_base$corl, fit_base$corh, file = ofileSout, append = TRUE
  )
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
    file = ofileSout, append = T
  )

  if (Ns > 0) {
    segments <- segment_bounds(Ips)
    adjustment <- function(i) mean(QMout$W[segments$start[i]:segments$end[i]])
    # Fcat 为分位组边界，DP 为各组均值差。
    write_qm_log(QMout, segments, Iseg.adj, Ny4a, adjustment, ofileSout)
  }

  fprintf(
    "#steps=  %d ; trend= %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f ) %.4f \n",
    Ns, otmp$trend, otmp$betaL, otmp$betaU, otmp$p.tr, cor, corl, corh,
    pcor, file = ofileSout, append = TRUE
  )
}

# 组装拟合序列，并恢复闰日数据。
build_fitted_table <- function(ori.itable, ooflg, otmp, EEB, adj, EBfull, meanhat0,
                          Ns, QMout, B, Y1, Ti, Ips, IY0, Mq, Iseg.adj, Nadj, Nt) {
  ori.itable <- as_tibble(ori.itable)
  n <- nrow(ori.itable)
  expand <- function(x) expand_values(x, ooflg, n)
  odata <- tibble(
    id = as.double(seq_len(n)),
    date = ymd2num(ori.itable),
    base = as.double(ori.itable$data),
    `base_trend+meanShift` = expand(otmp$meanhat + EEB),
    base_mean_adj = expand(adj),
    base_anomaly = expand(otmp$Y0),
    fit_of_base_anomaly = expand(otmp$meanhat),
    `annualC+meanShifts` = expand(otmp$meanhat + EBfull),
    QM_adjusted = if (Ns > 0 && QMout$Mq > 1) expand(B) else rep(NA_real_, n),
    fit_of_deseason_base = expand(meanhat0)
  )
  # odata[ooflg,1] <- Ti
  # odata[ooflg,3] <- round(otmp$Y0+EBfull,4)
  # odata[ooflg,9] <- round(oR,4)

  odata <- restore_leap_days(
    odata, ori.itable, IY0[Ips], Y1 - otmp$trend * Ti + EBfull,
    Ti, otmp$trend, QMout, Mq, Iseg.adj, Nadj, Nt,
    fitted_cols = c(
      "base_trend+meanShift", "fit_of_base_anomaly",
      "annualC+meanShifts", "fit_of_deseason_base"
    ),
    adjusted_cols = c("base_mean_adj", "base_anomaly")
  )
  odata
}

# 比较保留与移除候选突变点后的残差，计算 PMF 检验概率。
pmf_probability <- function(series, full_residuals, Ti, other_points, segment) {
  reduced_residuals <- LSmultiple(series, Ti, other_points)$resi
  full_sse <- sum(full_residuals[segment]^2)
  reduced_sse <- sum(reduced_residuals[segment]^2)
  df <- length(segment) - 3
  statistic <- (reduced_sse - full_sse) * df / full_sse
  if (statistic < 0) 0 else pf(statistic, 1, df)
}

# 逐个检验候选突变点，返回统计量与判定结果。
test_changepoints <- function(Y1, Ti, Ips, Ids, IY0, fit, plev, stepsizes, phi, PFmax) {
  Ns <- length(Ips) - 1
  RW <- LSmultiple(fit$W, Ti, Ips)$resi
  RWL <- LSmultiple(fit$WL, Ti, Ips)$resi
  RWU <- LSmultiple(fit$WU, Ti, Ips)$resi

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
      otmp <- LSmultiple(Y1, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      Pk0 <- Pk.PMFT(Nseg)
      PFx <- Fx * Pk0[Ic - I1 + 1]

      otmp <- LSmultiple(fit$W, Ti, Ips0)
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
