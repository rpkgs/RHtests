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
  # 1. 读取序列与检验参数；已有状态时复用。
  state <- if (is.null(.state)) {
    Read(InSeries, MissingValueCode, plev)
  } else {
    .state
  }
  check_rh_state(state, reference = FALSE, plev = plev)
  c(Y0, Ti, IY0, IY0flg, Imd, itable, ori.itable, ooflg) %<-%
    state@series[c("Y0", "Ti", "IY0", "IY0flg", "Imd", "itable", "ori.itable", "ooflg")]
  c(Icy, Nt) %<-% state@season[c("Icy", "Nt")]
  c(phi, PFmax) %<-% state@critical_values[c("phi", "PFmax")]

  # 2. 定位候选突变点，并选择订正基准段。
  if (Ny4a > 0 && Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Ny4a * Nt

  c(Ips, Ids, Ns) %<-% locate_changepoints(InCs, IY0)[c("Ips", "Ids", "Ns")]
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    message <- paste("StepSize: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
    stop(message, call. = FALSE)
  }
  Iseg.adj <- select_base_segment(Ips, Iadj)
  segments <- segment_bounds(Ips)

  # 3. 去除季节周期，分别拟合无突变和分段趋势。
  cycle <- RemoveSeasonal(Y0, Imd)
  Y1 <- cycle$anomaly
  seasonal_mean <- cycle$mu
  if (length(seasonal_mean) != length(Icy)) {
    stop("Annual cycle length differs from original dataset", call. = FALSE)
  }

  Ip0 <- N
  fit_base <- LSmultiRedCycle(Y1, Ti, Ip0, 1, IY0flg, seasonal_mean, Imd, Icy, Nt)
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  fit <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj, IY0flg, seasonal_mean, Imd, Icy, Nt)
  Y1 <- fit$Y0
  cor <- fit$cor
  corl <- fit$corl
  corh <- fit$corh
  pcor <- pt(abs(cor) * sqrt((N - 2) / (1 - cor^2)), N - 2)
  fitted_seasonal_mean <- fit$seasonal_mean
  seasonal_values <- fitted_seasonal_mean[match(Imd, Icy)]
  mean_seasonal_cycle <- mean(fitted_seasonal_mean, na.rm = TRUE)

  # 4. 去除趋势后进行分位数匹配，再加回趋势。
  Rb <- Y1 - fit$trend * Ti + seasonal_values
  QMout <- NULL
  if (Ns > 0) {
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj, Nt)
    B <- QMout$PA
  } else {
    B <- Rb
  }

  adj <- fit$Y0 + seasonal_values
  B <- B + fit$trend * Ti

  # 5. 计算相邻段的突变幅度，并将各段均值订正到基准段。
  stepsizes <- diff(fit$mu)
  if (Ns > 0) {
    adj <- add_segment_adjustment(adj, segments, fit$mu[Iseg.adj] - fit$mu)
  }
  Ehat <- mean(fit$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  # 6. 复检各候选突变点，恢复日期、缺测值和闰日。
  odata <- build_fitted_table(
    ori.itable, ooflg, fit, mean_seasonal_cycle, adj, seasonal_values, meanhat0,
    Ns, QMout, B, Y1, Ti, Ips, IY0, Mq, Iseg.adj, Nadj, Nt
  )

  d_TP <- test_changepoints(Y1, Ti, Ips, Ids, IY0, fit, plev, stepsizes, phi, PFmax)
  TP <- write_points_table(d_TP, paste0(output, "_fCs.txt"))

  odata <- dt_num2date(odata)

  # 7. 写出统计结果，按需绘制诊断图。
  ofileSout <- paste0(output, "_Fstat.txt")
  ofilePdf <- paste0(output, "_F.pdf")

  if (debug) {
    log_stepsize(
      ofileSout, plev, N, fit_base, QMout, Ny4a, Iseg.adj, Ips,
      fit, cor, corl, corh, pcor, IY0, stepsizes, d_TP
    )
  } else {
    cat("", file = ofileSout, append = FALSE)
  }
  if (is_plot) {
    plot_stepsize(
      cycle, ofilePdf, seasonal_values, mean_seasonal_cycle, B, QMout, Mq, Ns, adj,
      Ips, Iseg.adj, fit, ori.itable, ooflg
    )
  }

  if (GUI) {
    return(0)
  } else {
    list(data = odata, TP = TP)
  }
}
