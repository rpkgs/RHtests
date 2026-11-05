#' StepSize.wRef
#'
#' @param Bseries base series，文件路径或数据框。
#' @param Rseries reference series，文件路径或数据框。
#' @inheritParams FindU.wRef
#' @inheritParams FindU
#' @param InCs 待检验的突变点，data.frame 或文件路径，含列 `date` 与 `kind`。
#' @param debug
#' - `true`  : `_Fstat.txt` will be shown
#' - `false` : no `_Fstat.txt`
#' @example R/examples/run_ex02.R
#' @export
StepSize.wRef <- function(
  Bseries = NULL, Rseries = NULL, InCs, output,
  MissingValueCode = "-999.99",
  plev = 0.95, Iadj = 10000, Mq = 10, GUI = FALSE, Ny4a = 0,
  is_plot = TRUE, ..., debug = TRUE, .state = NULL
) {
  # 1. 读取基准与参考序列；已有状态时复用。
  state <- if (is.null(.state)) {
    Read.wRef(Bseries, Rseries, MissingValueCode, plev)
  } else {
    .state
  }
  check_rh_state(state, reference = TRUE, plev = plev)
  c(Y0, IY0, IY0flg, TiB, IYBflg, bdata, ori.bdata, owflg) %<-%
    state@series[c("Y0", "IY0", "IY0flg", "TiB", "IYBflg", "bdata", "ori.bdata", "owflg")]
  c(Icy, Nt) %<-% state@season[c("Icy", "Nt")]
  c(phi, PTmax) %<-% state@critical_values[c("phi", "PTmax")]

  # 2. 定位候选突变点，并选择订正基准段。
  if (Ny4a > 0 && Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Nt * Ny4a
  c(Ips, Ids, Ns) %<-% locate_changepoints(InCs, IY0)[c("Ips", "Ids", "Ns")]
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    message <- paste("StepSize.wRef: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
    stop(message, call. = FALSE)
  }

  Iseg.adj <- select_base_segment(Ips, Iadj)

  # 3. 在基准与参考序列的差值上复检突变点。
  fit_diff <- Rphi(Y0, Ips, Ns, IY0flg)
  df <- N - 2 - Ns
  p.cor.diff <- pt(abs(fit_diff$cor) * sqrt(df / (1 - fit_diff$cor^2)), df)
  d_TP <- test_changepoints_ref(Y0, IY0, Ips, Ids, fit_diff, plev, phi, PTmax)

  # 4. 计算差值序列各段均值，估计相对参考序列的偏移。
  muDif <- rep(0, Ns + 1)
  segments <- segment_bounds(Ips)
  omuDif <- Y0
  oY0 <- Y0
  for (i in 1:(Ns + 1)) {
    I0 <- segments$start[i]
    I2 <- segments$end[i]
    muDif[i] <- mean(Y0[I0:I2])
    omuDif[I0:I2] <- muDif[i]
  }

  # 将差值序列的突变位置映射回基准序列。
  Ips0 <- Ips
  IY1 <- ymd2num(bdata)
  IYM <- bdata[, 2] * 100 + bdata[, 3]
  inds <- sort(unique(IYM))
  rtmp <- cbind(1:length(IY1), IY1)

  Ips <- merge(IY0[Ips0], rtmp, by.x = 1, by.y = "IY1")[, 2]
  Ips[Ns + 1] <- length(IY1)
  segments <- segment_bounds(Ips)

  Ti <- TiB
  N <- length(TiB)
  IY0flg <- IYBflg

  # 5. 按差值均值订正基准序列，并拟合订正前后的趋势。
  cycle_base <- RemoveSeasonal(bdata$data, IYM)
  adjBase <- cycle_base$anomaly
  seasonal_values <- cycle_base$mu[match(IYM, inds)]

  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I0 <- segments$start[i]
      I2 <- segments$end[i]
      DeltaD <- muDif[i] - muDif[Iseg.adj]
      adjBase[I0:I2] <- adjBase[I0:I2] + seasonal_values[I0:I2] - DeltaD # adjBase is base series adj to Iseg.adj
    }
  }

  cycle_adjusted <- RemoveSeasonal(adjBase, IYM) # adjusted and de-seasonalized Base series
  adjusted_seasonal_mean <- cycle_adjusted$mu
  mean_adjusted_cycle <- mean(adjusted_seasonal_mean)
  Aadj <- cycle_adjusted$anomaly # de-seasonalize Badj
  Ipd <- c(N)
  fit_adjusted <- LSmultiRedCycle(Aadj, Ti, Ipd, Iseg.adj, IY0flg, adjusted_seasonal_mean, IYM, Icy, Nt)

  # de-seasonalized Base series
  Ipd <- length(cycle_base$anomaly)
  fit_base <- LSmultiRedCycle(
    cycle_base$anomaly, Ti, Ipd, Iseg.adj, IY0flg, cycle_base$mu, IYM, Icy, Nt
  )
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  # 6. 拟合基准序列的分段趋势，迭代更新季节周期。
  EPB <- mean(bdata[, "data"], na.rm = TRUE)
  fit_result <- fit_reference_cycle(
    cycle_base$anomaly, seasonal_values, adjusted_seasonal_mean, IYM, inds, Ti, Ips, IY0flg, Iseg.adj, Nt
  )
  c(Base, seasonal_values, fitted_seasonal_mean, otmp) %<-%
    fit_result[c("base", "seasonal", "cycle", "fit")]
  sig <- otmp$sig
  df <- N - 2 - Nt - Ns
  p.cor <- pt(abs(otmp$cor) * sqrt(df / (1 - otmp$cor^2)), df)
  Ehat <- mean(otmp$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  # 7. 按分位数、基准序列均值和差值均值分别订正。
  QMout <- NULL
  if (Ns > 0) {
    QMout <- QMadjGaussian.wRef(bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a)
    B <- QMout$PA
  } else {
    B <- Base
  }

  adj <- Base + seasonal_values
  adjB <- Base + seasonal_values
  meanhatD <- rep(0, N)
  if (Ns > 0) {
    delta <- sig[Iseg.adj + 1] - sig[-1]
    delta_diff <- muDif[Iseg.adj] - muDif
    adj <- add_segment_adjustment(adj, segments, delta)
    adjB <- add_segment_adjustment(adjB, segments, delta_diff)
    meanhatD <- add_segment_adjustment(
      mean_adjusted_cycle + fit_adjusted$mu[1] + fit_adjusted$trend * Ti, segments, -delta_diff
    )
  } else {
    meanhatD <- EPB + fit_adjusted$mu[1] + fit_adjusted$trend * Ti
  }

  # 各段系数相对首段；后续突变幅度取相邻系数之差。
  stepsizes <- if (Ns > 0) c(sig[3], diff(sig[-(1:2)])) else numeric()

  # 8. 恢复日期、缺测值和闰日，组装输出序列。
  odata <- build_ref_fitted_table(
    ori.bdata, owflg, Ti, IY0[Ips0], Base, seasonal_values, fitted_seasonal_mean, otmp, meanhatD,
    adj, adjB, meanhat0, B, QMout, Mq, Iseg.adj, Nadj, Nt
  )

  # 9. 写出统计结果，按需绘制诊断图。
  ofileSout <- paste0(output, "_Fstat.txt")
  ofileIout <- paste0(output, "_fCs.txt")
  if (debug) {
    log_stepsize_ref(
      ofileSout, Bseries, Rseries, fit_diff, p.cor.diff, d_TP,
      fit_base, fit_adjusted, muDif, QMout, segments, Iseg.adj, Ny4a,
      otmp, p.cor, stepsizes
    )
  } else {
    cat("", file = ofileSout, append = FALSE)
  }
  fprintf("%s changepoints in Series %s \n", Ns, Bseries, file = ofileIout, append = FALSE)
  d_TP <- write_points_table(
    d_TP, ofileIout, reference = TRUE, append = TRUE, stepsizes = stepsizes
  )
  if (is_plot) {
    plot_StepSize.ref(
      output, Base, seasonal_values, fitted_seasonal_mean, B,
      oY0, omuDif, otmp, meanhatD,
      QMout, Mq, Ns, adj, adjB, Ips, Iseg.adj, IY0, ori.bdata, owflg
    )
  }

  if (GUI) {
    return(0)
  } else {
    list(data = odata, TP = d_TP)
  }
}
