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
  plev = 0.95, Iadj = 10000, Mq = 10, GUI = F, Ny4a = 0,
  is_plot = TRUE, ..., debug = TRUE, .state = NULL
) {
  # 读取计算状态，或使用上游传入的状态。
  state <- if (is.null(.state)) {
    Read.wRef(Bseries, Rseries, MissingValueCode, plev)
  } else .state
  check_rh_state(state, reference = TRUE, plev = plev)
  c(Y0, IY0, IY0flg, TiB, IYBflg, bdata, ori.bdata, owflg) %<-%
    state@series[c("Y0", "IY0", "IY0flg", "TiB", "IYBflg", "bdata", "ori.bdata", "owflg")]
  c(Icy, Nt) %<-% state@season[c("Icy", "Nt")]
  c(phi, PTmax) %<-% state@critical_values[c("phi", "PTmax")]

  # 定位候选突变点。
  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Nt * Ny4a
  c(Ips, Ids, Ns) %<-% locate_changepoints(InCs, IY0)[c("Ips", "Ids", "Ns")]
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    message <- paste("StepSize.wRef: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
    stop(message, call. = FALSE)
  }

  # if(Iadj==1|Iseg.longest==0) Iseg.adj<-Ns+1 else Iseg.adj<-Iseg.longest
  Iseg.adj <- select_base_segment(Ips, Iadj)

  # 在基准与参考序列的差值上检验突变点。
  fit_diff <- Rphi(Y0, Ips, Ns, IY0flg)
  df <- N - 2 - Ns
  p.cor <- pt(abs(fit_diff$cor) * sqrt(df / (1 - fit_diff$cor^2)), df)
  ofileIout <- paste0(output, "_fCs.txt")
  ofilePdf <- paste0(output, "_F.pdf")
  ofileSout <- paste0(output, "_Fstat.txt")
  # ofileRout<-paste(output,"_Base_Ref.fitUDfinal",sep="")
  file.create(ofileSout)
  file.create(ofileIout)
  # file.create(ofileRout)

  if (debug) log_ref_inputs(ofileSout, Bseries, Rseries, fit_diff, p.cor)

  fprintf("%s changepoints in Series %s \n", Ns, Bseries, file = ofileIout)
  d_TP <- test_changepoints_ref(Y0, IY0, Ips, Ids, fit_diff, plev, phi, PTmax)
  if (debug) write_test_log(d_TP, "PMT", ofileSout)

  # 计算差值序列各段均值，用于估计相对参考序列的偏移。
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
  # write.table(cbind(IY0,round(oY0,4),round(omuDif,4),round(Wo,4)),
  #             file=ofileRout,col.names=F,row.names=F)

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

  cycle_base <- rmCycle(bdata)
  adjBase <- cycle_base$Base
  EB <- cycle_base$EB[match(IYM, inds)]

  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I0 <- segments$start[i]
      I2 <- segments$end[i]
      DeltaD <- muDif[i] - muDif[Iseg.adj]
      adjBase[I0:I2] <- adjBase[I0:I2] + EB[I0:I2] - DeltaD # adjBase is base series adj to Iseg.adj
    }
  }

  Ydata <- cbind(bdata[, 1:3], adjBase)
  cycle_adjusted <- rmCycle(Ydata) # adjusted and de-seasonalized Base series
  EB0 <- cycle_adjusted$EB
  EEBd <- mean(EB0)
  Aadj <- cycle_adjusted$Base # de-seasonalize Badj
  Ipd <- c(N)
  # dtmp<-LSmultipleRed(Aadj,Ti,Ipd)
  fit_adjusted <- LSmultiRedCycle(Aadj, Ti, Ipd, Iseg.adj, IY0flg, EB0, IYM, Icy, Nt)

  # de-seasonalized Base series
  Ipd <- length(cycle_base$Base)
  # dtmp   <- LSmultipleRed(tbase,Ti,Ipd)
  fit_base <- LSmultiRedCycle(cycle_base$Base, Ti, Ipd, Iseg.adj, IY0flg, cycle_base$EB, IYM, Icy, Nt)
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  if (debug) log_ref_fits(ofileSout, fit_base, fit_adjusted, muDif)

  # 拟合基准序列，并迭代更新季节周期直至收敛。
  EPB <- mean(bdata[, "data"], na.rm = TRUE)
  fit_result <- fit_reference_cycle(
    cycle_base$Base, EB, EB0, IYM, inds, Ti, Ips, IY0flg, Iseg.adj, Nt
  )
  Base <- fit_result$base
  EB <- fit_result$seasonal
  EB1 <- fit_result$cycle
  otmp <- fit_result$fit
  sig <- otmp$sig
  df <- N - 2 - Nt - Ns
  p.cor <- pt(abs(otmp$cor) * sqrt(df / (1 - otmp$cor^2)), df)
  Ehat <- mean(otmp$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat
  # 分别按分位数匹配和各段均值订正基准序列。
  QMout <- NULL
  if (Ns > 0) {
    #   Rb<-Base-otmp$trend*Ti+EB
    #   QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    QMout <- QMadjGaussian.wRef(bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a)
    B <- QMout$PA

    if (debug) {
      adjustment <- function(i) QMout$AdjM[i]
      # Fcat 为分位组边界，DP 为各组均值差。
      write_qm_log(QMout, segments, Iseg.adj, Ny4a, adjustment, ofileSout)
    }
  } else {
    B <- Base
  }
  # else B<-Base-otmp$trend*Ti+EB
  # B<-B+otmp$trend*Ti

  adj <- Base + EB
  adjB <- Base + EB
  meanhatD <- rep(0, N)
  if (Ns > 0) {
    delta <- sig[Iseg.adj + 1] - sig[-1]
    delta_diff <- muDif[Iseg.adj] - muDif
    adj <- add_segment_adjustment(adj, segments, delta)
    adjB <- add_segment_adjustment(adjB, segments, delta_diff)
    meanhatD <- add_segment_adjustment(
      EEBd + fit_adjusted$mu[1] + fit_adjusted$trend * Ti, segments, -delta_diff
    )
  } else {
    meanhatD <- EPB + fit_adjusted$mu[1] + fit_adjusted$trend * Ti
  }

  if (is_plot) {
    plot_StepSize.ref(
      output, Base, EB, EB1, B,
      oY0, omuDif, otmp, meanhatD,
      QMout, Mq, Ns, adj, adjB, Ips, Iseg.adj, IY0, ori.bdata, owflg
    )
  }

  # 各段系数相对首段；后续突变幅度取相邻系数之差。
  stepsizes <- if (Ns > 0) c(sig[3], diff(sig[-(1:2)])) else numeric()
  if (debug) log_ref_result(ofileSout, otmp, p.cor, Ns, stepsizes)

  # 恢复原始日期和缺测位置，组装输出序列。
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
    QM_adjusted = if (Ns > 0) expand(B) else rep(NA_real_, n),
    fit_of_deseason_base = expand(meanhat0)
  )
  # odata[owflg,3] <- round(Base+EB,4)
  # odata[owflg,11]<-round(Ro,4)

  odata <- restore_leap_days(
    odata, ori.bdata, IY0[Ips0], Base - otmp$trend * Ti + EB,
    Ti, otmp$trend, QMout, Mq, Iseg.adj, Nadj, Nt,
    fitted_cols = c("base_trend+meanShift.ref", "base_trend+meanShift.base",
                    "fit_of_base_anomaly", "annualC+meanShifts", "fit_of_deseason_base"),
    adjusted_cols = c("base_mean_adj.ref", "base_mean_adj.base", "base_anomaly")
  )

  odata <- dt_num2date(odata)
  # write.table(file=paste0(output,"_F.dat"), odata, col.names=TRUE, row.names=F, na=MissingValueCode)

  d_TP <- write_points_table(d_TP, ofileIout, reference = TRUE, stepsizes = stepsizes)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    # cat("StepSize.wRef finished successfully...\n")
    list(data = odata, TP = d_TP)
  }
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
