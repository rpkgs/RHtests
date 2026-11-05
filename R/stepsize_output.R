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

# 组装参考序列流程的拟合结果，并恢复闰日数据。
build_ref_fitted_table <- function(ori.bdata, owflg, Ti, break_dates, Base, EB, EB1,
                                   otmp, meanhatD, adj, adjB, meanhat0, B, QMout,
                                   Mq, Iseg.adj, Nadj, Nt) {
  ori.bdata <- as_tibble(ori.bdata)
  n <- nrow(ori.bdata)
  expand <- function(x) expand_values(x, owflg, n)
  id <- rep(NA_real_, n)
  id[owflg] <- Ti
  odata <- tibble(
    id = id,
    date = ymd2num(ori.bdata),
    base = as.double(ori.bdata$data),
    `base_trend+meanShift.ref` = expand(meanhatD),
    base_mean_adj.ref = expand(adjB),
    `base_trend+meanShift.base` = expand(otmp$meanhat + mean(EB1)),
    base_mean_adj.base = expand(adj),
    base_anomaly = expand(Base),
    fit_of_base_anomaly = expand(otmp$meanhat),
    `annualC+meanShifts` = expand(otmp$meanhat + EB),
    # if(QMout$Mq>1) # vital change here, kongdd, 2020-07-24
    QM_adjusted = if (length(break_dates) > 1) expand(B) else rep(NA_real_, n),
    fit_of_deseason_base = expand(meanhat0)
  )


  odata <- restore_leap_days(
    odata, ori.bdata, break_dates, Base - otmp$trend * Ti + EB,
    Ti, otmp$trend, QMout, Mq, Iseg.adj, Nadj, Nt,
    fitted_cols = c(
      "base_trend+meanShift.ref", "base_trend+meanShift.base",
      "fit_of_base_anomaly", "annualC+meanShifts", "fit_of_deseason_base"
    ),
    adjusted_cols = c("base_mean_adj.ref", "base_mean_adj.base", "base_anomaly")
  )
  dt_num2date(odata)
}
